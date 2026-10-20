"""One running copy: a second start hands its files to the window already open.

Double-clicking a sequence file starts ``Vettore <file>``. When a copy is already running, the
new process gives the file to it and exits, so the file opens as a tab in the window the user
already has. How:

- The running copy records its local server (port, session token, process id) in
  ``instance.json`` under the settings home (``$VETTORE_HOME`` or ``~/.vettore``), readable by
  the user only, since the token guards the local API.
- A new start reads the record. It is ignored, and later replaced, when the process is gone or
  its server does not answer with that token (a stale record left by a crash or a restart).
- Otherwise the new start posts its files to ``/api/hand-over``. The running copy queues them
  with the files from its own command line (:class:`Inbox`); its page takes them within a
  second. The post returns once the page has taken them, or reports that it did not (a page
  that is frozen or was closed): the new start then opens its own window instead, so a
  double-click never ends in nothing.
- Several files selected and opened together start several processes at once: an exclusive
  lock on ``instance.lock`` makes them start one after the other, so the first becomes the
  running copy and the others hand their files to it.

This module imports nothing heavy, so a start that only hands over its files exits quickly.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import sys
import threading
import time
import urllib.request
from collections.abc import Iterable, Iterator
from dataclasses import asdict, dataclass
from pathlib import Path

from .. import settings

log = logging.getLogger(__name__)

TOKEN_HEADER = "X-Vettore-Token"
RECORD_FILE = "instance.json"
LOCK_FILE = "instance.lock"

#: How long the running copy waits for its page to take handed-over files (seconds). A page
#: polls every half second; a window that is still loading needs a few seconds more.
HAND_OVER_WAIT = 12.0

#: Dispatched on the page when files were handed over, so it takes them at once.
FILES_EVENT = "vettore:files-handed-over"


@dataclass(frozen=True)
class Record:
    """Where a running copy can be reached."""

    port: int
    token: str
    pid: int
    executable: str = ""

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.port}"


def record_path() -> Path:
    return settings.home() / RECORD_FILE


def read_record() -> Record | None:
    try:
        data = json.loads(record_path().read_text(encoding="utf-8"))
        return Record(port=int(data["port"]), token=str(data["token"]), pid=int(data["pid"]),
                      executable=str(data.get("executable", "")))
    except (OSError, ValueError, KeyError, TypeError):
        return None


def register(port: int, token: str) -> Record:
    """Record this process as the running copy (replacing any earlier record)."""
    record = Record(port=port, token=token, pid=os.getpid(), executable=sys.executable)
    path = record_path()
    tmp = path.with_name(f".{path.name}.tmp")
    try:
        # Readable by the user only: the token gives access to the local API. (A leftover
        # temporary file would keep its own permissions, so it goes first.)
        tmp.unlink(missing_ok=True)
        fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(asdict(record), stream)
        for attempt in range(5):
            try:
                os.replace(tmp, path)
                break
            except PermissionError:   # Windows: another start is reading the record
                if attempt == 4:
                    raise
                time.sleep(0.05)
    except OSError:
        log.exception("Could not record the running copy; files opened later start a new window")
        with contextlib.suppress(OSError):
            tmp.unlink()
    return record


def unregister(record: Record) -> None:
    """Remove the record if it is still this process's (a later copy may have replaced it)."""
    current = read_record()
    if current is not None and (current.pid, current.token) == (record.pid, record.token):
        with contextlib.suppress(OSError):
            record_path().unlink()


def pid_alive(pid: int) -> bool:
    """Whether a process with this id exists (it may belong to another program: the server
    check in :func:`find_live` tells)."""
    if pid <= 0:
        return False
    if sys.platform == "win32":
        return _windows_pid_alive(pid)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:   # exists, owned by someone else
        return True
    return True


def _windows_pid_alive(pid: int) -> bool:
    # os.kill(pid, 0) would send Ctrl+C on Windows: ask the kernel instead.
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    process_query_limited_information, still_active, access_denied = 0x1000, 259, 5
    handle = kernel32.OpenProcess(process_query_limited_information, False, pid)
    if not handle:
        return ctypes.get_last_error() == access_denied
    try:
        code = wintypes.DWORD()
        if not kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
            return True
        return code.value == still_active
    finally:
        kernel32.CloseHandle(handle)


def _request(record: Record, path: str, body: object | None, timeout: float) -> dict:
    request = urllib.request.Request(
        f"{record.url}/api/{path}",
        data=None if body is None else json.dumps(body).encode("utf-8"),
        headers={TOKEN_HEADER: record.token, "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def answers(record: Record, timeout: float = 2.0) -> bool:
    """The recorded server is up and accepts the recorded token."""
    try:
        return _request(record, "health", None, timeout) == {"status": "ok"}
    except (OSError, ValueError):
        return False


def find_live() -> Record | None:
    """The running copy, or None when there is none (no record, or a stale one)."""
    record = read_record()
    if record is None:
        return None
    if not pid_alive(record.pid):
        log.info("Ignoring the record of a Vettore that is no longer running (pid %d)",
                 record.pid)
        return None
    if not answers(record):
        log.info("Ignoring the record of pid %d: its server does not answer", record.pid)
        return None
    return record


def hand_over(record: Record, paths: Iterable[Path], wait: float = HAND_OVER_WAIT) -> bool:
    """Give files to the running copy; True once its page has taken them. With no files, the
    running copy only brings its window to the front."""
    allow_foreground(record.pid)
    body = {"paths": [str(p) for p in paths]}
    try:
        return _request(record, "hand-over", body, wait + 5).get("delivered") is True
    except (OSError, ValueError):
        log.exception("Could not hand the files over to the running Vettore (pid %d)",
                      record.pid)
        return False


def allow_foreground(pid: int) -> None:
    """Windows lets only the program the user just started take the foreground: pass that
    right on to the running copy, so it can bring its window to the front."""
    if sys.platform != "win32":
        return
    try:
        import ctypes

        ctypes.windll.user32.AllowSetForegroundWindow(pid)
    except Exception:  # cosmetic: the files still open
        log.debug("AllowSetForegroundWindow failed", exc_info=True)


@contextlib.contextmanager
def startup_lock(timeout: float = 30.0) -> Iterator[bool]:
    """Hold an exclusive lock while deciding whether to hand over or to become the running
    copy. Yields whether the lock was obtained; after ``timeout`` the start goes ahead
    without it."""
    path = settings.home() / LOCK_FILE
    try:
        fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o600)
    except OSError:
        log.exception("Could not open %s; starting without it", path)
        yield False
        return
    acquired = False
    try:
        deadline = time.monotonic() + timeout
        while True:
            try:
                _lock(fd)
                acquired = True
                break
            except OSError:
                if time.monotonic() > deadline:
                    log.warning("Another Vettore kept %s for %.0f s; starting anyway", path,
                                timeout)
                    break
                time.sleep(0.1)
        yield acquired
    finally:
        if acquired:
            with contextlib.suppress(OSError):
                _unlock(fd)
        os.close(fd)


def _lock(fd: int) -> None:
    if sys.platform == "win32":
        import msvcrt

        os.lseek(fd, 0, os.SEEK_SET)
        msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
    else:
        import fcntl

        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)


def _unlock(fd: int) -> None:
    if sys.platform == "win32":
        import msvcrt

        os.lseek(fd, 0, os.SEEK_SET)
        msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
    else:
        import fcntl

        fcntl.flock(fd, fcntl.LOCK_UN)


def files_script() -> str:
    """JavaScript telling the page that files are waiting."""
    return f"window.dispatchEvent(new CustomEvent({json.dumps(FILES_EVENT)}))"


class Inbox:
    """Files waiting for the page: those from the command line, then those handed over by later
    starts. The page takes everything waiting at once (``GET /api/startup-files``)."""

    def __init__(self, paths: Iterable[Path] = ()) -> None:
        self._cond = threading.Condition()
        self._items: list[tuple[int, Path]] = [(0, p) for p in paths]
        self._next = 1
        self._taken = 0   # every ticket up to this one has been taken by the page

    def put(self, paths: Iterable[Path]) -> int:
        """Queue files; returns the ticket to wait on. An empty hand-over gets a ticket too:
        it is taken at the page's next look."""
        with self._cond:
            ticket = self._next
            self._next += 1
            self._items += [(ticket, p) for p in paths]
            return ticket

    def take(self) -> list[Path]:
        with self._cond:
            paths = [p for _, p in self._items]
            self._items.clear()
            self._taken = self._next - 1
            self._cond.notify_all()
            return paths

    def wait_or_withdraw(self, ticket: int, timeout: float) -> bool:
        """True once the page took the ticket's files. Otherwise they are withdrawn, in the
        same step, so they can never open both here and in the window the caller starts."""
        with self._cond:
            if self._cond.wait_for(lambda: self._taken >= ticket, timeout):
                return True
            self._items = [(t, p) for t, p in self._items if t != ticket]
            return False
