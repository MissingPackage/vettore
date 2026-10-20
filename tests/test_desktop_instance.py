"""One running copy: a second start hands its files to the window already open.

The last tests start real processes: a running copy (server only, no window) and later starts
with files, the way a double-click starts ``Vettore <file>``. The test plays the page, which
asks ``/api/startup-files`` every half second.
"""

import json
import os
import socket
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from vettore.desktop import instance
from vettore.desktop.instance import Inbox, Record
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"
H = {TOKEN_HEADER: TOKEN}


@pytest.fixture(autouse=True)
def private_home(tmp_path, monkeypatch):
    home = tmp_path / "vettore-home"
    monkeypatch.setenv("VETTORE_HOME", str(home))
    return home


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _dead_pid() -> int:
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    return child.pid


# ---- the queue of files waiting for the page ---------------------------------------------

def test_the_page_takes_command_line_files_then_handed_over_ones():
    inbox = Inbox([Path("a.gb")])
    ticket = inbox.put([Path("b.dna"), Path("c.fasta")])
    assert inbox.take() == [Path("a.gb"), Path("b.dna"), Path("c.fasta")]
    assert inbox.wait_or_withdraw(ticket, 0) is True
    assert inbox.take() == []


def test_files_the_page_does_not_take_in_time_are_withdrawn():
    inbox = Inbox()
    ticket = inbox.put([Path("late.gb")])
    assert inbox.wait_or_withdraw(ticket, 0.1) is False
    assert inbox.take() == []   # never opened here: the second start opens them itself


def test_a_hand_over_without_files_is_delivered_at_the_page_next_look():
    inbox = Inbox()
    ticket = inbox.put([])
    threading.Timer(0.1, inbox.take).start()
    assert inbox.wait_or_withdraw(ticket, 5) is True


def test_hand_over_endpoint_waits_for_the_page():
    shown = []
    app = create_app(TOKEN, None, [], on_hand_over=lambda: shown.append(True),
                     hand_over_wait=10)
    sample = DATA / "sample-d.dna"
    taken: list[dict] = []

    def page():
        with TestClient(app, base_url="http://127.0.0.1") as c:
            deadline = time.monotonic() + 10
            while not taken and time.monotonic() < deadline:
                taken.extend(c.get("/api/startup-files", headers=H).json()["results"])
                time.sleep(0.05)

    reader = threading.Thread(target=page)
    reader.start()
    with TestClient(app, base_url="http://127.0.0.1") as c:
        res = c.post("/api/hand-over", headers=H, json={"paths": [str(sample)]})
    reader.join()
    assert res.json() == {"delivered": True}
    assert shown == [True]
    assert [r["path"] for r in taken] == [str(sample)]
    assert len(taken[0]["documents"]) == 1


def test_hand_over_to_a_page_that_never_looks_is_not_delivered():
    app = create_app(TOKEN, None, [], hand_over_wait=0.2)
    c = TestClient(app, base_url="http://127.0.0.1")
    res = c.post("/api/hand-over", headers=H, json={"paths": [str(DATA / "sample-d.dna")]})
    assert res.json() == {"delivered": False}
    assert c.get("/api/startup-files", headers=H).json() == {"results": [], "errors": []}
    assert c.post("/api/hand-over", json={"paths": []}).status_code == 403   # token required


# ---- the record of the running copy ------------------------------------------------------

def test_the_record_is_private_and_removed_only_by_its_owner():
    record = instance.register(4567, "secret")
    assert instance.read_record() == record
    assert record.pid == os.getpid()
    if sys.platform != "win32":
        assert (instance.record_path().stat().st_mode & 0o777) == 0o600
    instance.unregister(Record(port=4567, token="another", pid=record.pid))
    assert instance.read_record() == record          # someone else's: kept
    instance.unregister(record)
    assert instance.read_record() is None


def test_a_damaged_record_is_ignored(private_home):
    private_home.mkdir(parents=True, exist_ok=True)
    instance.record_path().write_text("{not json", encoding="utf-8")
    assert instance.read_record() is None
    assert instance.find_live() is None


def test_process_liveness():
    assert instance.pid_alive(os.getpid())
    assert not instance.pid_alive(_dead_pid())
    assert not instance.pid_alive(0)


def _write_record(record: Record) -> None:
    instance.record_path().write_text(json.dumps(record.__dict__), encoding="utf-8")


def test_stale_records_are_ignored():
    instance.register(1, "x")   # creates the home folder
    _write_record(Record(port=_free_port(), token="x", pid=_dead_pid()))
    assert instance.find_live() is None               # the process is gone
    # Alive (this very process) but nothing answers on the port: also stale.
    _write_record(Record(port=_free_port(), token="x", pid=os.getpid()))
    assert instance.find_live() is None


def test_overlapping_starts_take_turns():
    with instance.startup_lock() as first, instance.startup_lock(timeout=0.3) as second:
        assert first is True and second is False
    with instance.startup_lock(timeout=0.3) as again:
        assert again is True


def test_the_page_is_told_when_files_arrive():
    # frontend/src/App.tsx listens for this event.
    assert instance.FILES_EVENT == "vettore:files-handed-over"
    assert instance.files_script() == (
        'window.dispatchEvent(new CustomEvent("vettore:files-handed-over"))')


# ---- real processes ---------------------------------------------------------------------

def _get(port: int, path: str, token: str = TOKEN) -> dict:
    req = urllib.request.Request(f"http://127.0.0.1:{port}/api/{path}",
                                 headers={TOKEN_HEADER: token})
    with urllib.request.urlopen(req, timeout=5) as res:
        return json.loads(res.read())


def _wait_up(proc: subprocess.Popen, port: int, token: str = TOKEN) -> None:
    """Until the server answers and the process has recorded itself as the running copy
    (right after its server is up)."""
    deadline = time.monotonic() + 90
    while time.monotonic() < deadline:
        assert proc.poll() is None, "the running copy exited early"
        try:
            if _get(port, "health", token) == {"status": "ok"}:
                record = instance.read_record()
                if record is not None and record.pid == proc.pid:
                    return
        except OSError:
            pass
        time.sleep(0.2)
    raise AssertionError("the running copy did not start")


class Page:
    """Plays the page of a running copy: asks for waiting files every 0.1 s."""

    def __init__(self, port: int, token: str = TOKEN) -> None:
        self.port, self.token = port, token
        self.results: list[dict] = []
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def _run(self) -> None:
        while not self._stop.is_set():
            try:
                self.results += _get(self.port, "startup-files", self.token)["results"]
            except OSError:
                pass
            time.sleep(0.1)

    def close(self) -> None:
        self._stop.set()
        self._thread.join(5)


class Starts:
    """Vettore processes started by a test, each logging to its own file."""

    def __init__(self, home: Path, logs: Path) -> None:
        self.env = {**os.environ, "VETTORE_HOME": str(home)}
        self.logs = logs
        self.procs: list[tuple[subprocess.Popen, Path]] = []

    def start(self, *args: str) -> subprocess.Popen:
        # --no-window keeps the test off the display; --single-instance gives it the decision
        # a start with a window makes (hand over, or become the running copy). A start that
        # wrongly did not hand over would serve forever: the timeouts below catch it.
        path = self.logs / f"start-{len(self.procs)}.log"
        with open(path, "wb") as log:
            proc = subprocess.Popen([sys.executable, "-m", "vettore.app", "--no-window",
                                     "--single-instance", *args],
                                    env=self.env, stdout=subprocess.DEVNULL, stderr=log)
        self.procs.append((proc, path))
        return proc

    def log(self, proc: subprocess.Popen) -> str:
        path = next(p for q, p in self.procs if q is proc)
        return path.read_text(encoding="utf-8", errors="replace")

    @staticmethod
    def stop(proc: subprocess.Popen) -> None:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)

    def stop_all(self) -> None:
        for proc, _ in self.procs:
            self.stop(proc)


@pytest.fixture
def starts(private_home, tmp_path):
    logs = tmp_path / "logs"
    logs.mkdir()
    s = Starts(private_home, logs)
    yield s
    s.stop_all()


def test_a_second_start_hands_its_file_to_the_running_copy_and_exits(starts, tmp_path):
    port = _free_port()
    first = starts.start("--port", str(port), "--token", TOKEN)
    _wait_up(first, port)
    assert instance.read_record().pid == first.pid
    page = Page(port)
    try:
        sample = tmp_path / "colleague.dna"
        sample.write_bytes((DATA / "sample-d.dna").read_bytes())
        second = starts.start(str(sample))
        assert second.wait(timeout=60) == 0, starts.log(second)
        assert "Handed 1 file(s) to the running Vettore" in starts.log(second)
        assert [r["path"] for r in page.results] == [str(sample)]
        assert page.results[0]["documents"][0]["sequence"]
        assert instance.read_record().pid == first.pid     # still the running copy

        # Started again with no file: the running copy is shown, no second copy starts.
        third = starts.start()
        assert third.wait(timeout=60) == 0, starts.log(third)
        assert "Handed 0 file(s)" in starts.log(third)
        assert len(page.results) == 1
    finally:
        page.close()


def test_a_start_whose_files_the_page_does_not_take_runs_on_its_own(starts, tmp_path):
    """A running copy whose page never looks (frozen, or a browser tab that was closed): the
    new start does not end in nothing, it runs on its own and becomes the running copy."""
    port = _free_port()
    first = starts.start("--port", str(port), "--token", TOKEN)
    _wait_up(first, port)
    sample = tmp_path / "plasmid.gb"
    sample.write_bytes((DATA / "MT612431.1.gb").read_bytes())
    port2 = _free_port()
    started = time.monotonic()
    second = starts.start("--port", str(port2), "--token", "other", str(sample))
    _wait_up(second, port2, "other")
    assert time.monotonic() - started >= instance.HAND_OVER_WAIT
    assert "did not open the files" in starts.log(second)
    assert instance.read_record().pid == second.pid
    assert [r["path"] for r in _get(port2, "startup-files", "other")["results"]] == [str(sample)]
    assert _get(port, "startup-files")["results"] == []   # withdrawn from the first


def test_a_stale_record_is_replaced_by_the_next_start(starts, tmp_path):
    port = _free_port()
    first = starts.start("--port", str(port), "--token", TOKEN)
    _wait_up(first, port)
    starts.stop(first)                 # killed: its record stays behind, pointing at nothing
    assert instance.read_record().pid == first.pid

    sample = tmp_path / "plasmid.gb"
    sample.write_bytes((DATA / "MT612431.1.gb").read_bytes())
    port2 = _free_port()
    successor = starts.start("--port", str(port2), "--token", "other", str(sample))
    _wait_up(successor, port2, "other")
    record = instance.read_record()
    assert (record.pid, record.port) == (successor.pid, port2)
    results = _get(port2, "startup-files", "other")["results"]
    assert [r["path"] for r in results] == [str(sample)]
