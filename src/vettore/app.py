"""Desktop entry point: starts the local API server and opens the UI in a native window.

A start with a copy of Vettore already running hands its files to that copy and exits
(``desktop.instance``). That decision comes first and needs none of the biology engine, so the
engine (``server``, ``fileio``, uvicorn) is imported only once this process is going to run.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import logging
import os
import secrets
import socket
import sys
import threading
import time
import urllib.request
import webbrowser
from pathlib import Path
from typing import Any

from . import settings
from .desktop import close_guard, instance
from .desktop.window import WindowPresence

log = logging.getLogger("vettore")

TITLE = "Vettore"

# Save dialog filters per kind of file; the engine picks the format from the extension.
SAVE_TYPES = {
    "genbank": "GenBank (*.gb;*.gbk)",
    "fasta": "FASTA (*.fasta;*.fa)",
    "dna": "DNA file (*.dna)",
    "svg": "SVG image (*.svg)",
    "png": "PNG image (*.png)",
}
SEQUENCE_KINDS = ("genbank", "fasta", "dna")

# Browser event carrying the full paths of files dropped on the native window.
DROP_EVENT = "vettore:dropped-paths"


def file_types() -> tuple[str, ...]:
    """Open dialog filters: every extension the engine reads (it also recognises files by their
    content, so "All files" opens files with a wrong or missing extension). "*" matches names
    without a dot."""
    from . import fileio

    return (
        "Sequence files (" + ";".join(f"*{ext}" for ext in fileio.FORMATS) + ")",
        "All files (*.*;*)",
    )


def save_file_types(kind: str) -> tuple[str, ...]:
    """Filters for the save dialog, the requested kind first (Windows appends its extension)."""
    if kind not in SEQUENCE_KINDS:
        return (SAVE_TYPES.get(kind, "All files (*.*;*)"),)
    return (SAVE_TYPES[kind], *(SAVE_TYPES[k] for k in SEQUENCE_KINDS if k != kind))


def dropped_paths(event: dict[str, Any]) -> dict[str, list[str]]:
    """Split the files of a drop event into full paths and names pywebview could not resolve.

    pywebview adds ``pywebviewFullPath`` to the files of its own copy of the event (the
    browser's File objects never carry a path)."""
    files = (event.get("dataTransfer") or {}).get("files") or []
    paths: list[str] = []
    unresolved: list[str] = []
    for item in files:
        if full := item.get("pywebviewFullPath"):
            paths.append(full)
        elif name := item.get("name"):
            unresolved.append(name)
    return {"paths": paths, "unresolved": unresolved}


def drop_script(payload: dict[str, list[str]]) -> str:
    detail = json.dumps(payload)
    return f"window.dispatchEvent(new CustomEvent({json.dumps(DROP_EVENT)}, {{detail: {detail}}}))"


def _bind_drop(window) -> None:
    """Forward dropped files' full paths to the page, so a dropped document remembers where it
    came from and Save writes back in place. Called on every page load."""
    from webview.dom import DOMEventHandler

    def on_drop(event: dict[str, Any]) -> None:
        payload = dropped_paths(event)
        if payload["paths"] or payload["unresolved"]:
            window.run_js(drop_script(payload))

    try:
        window.dom.document.events.drop += DOMEventHandler(on_drop)
        window.run_js("window.__vettoreNativeDrop = true")
    except Exception:  # without it, dropped files are uploaded by content instead
        log.exception("Could not listen for dropped files")


def static_dir() -> Path | None:
    candidates = []
    if env := os.environ.get("VETTORE_STATIC"):
        candidates.append(Path(env))
    if bundle := getattr(sys, "_MEIPASS", None):  # PyInstaller
        candidates.append(Path(bundle) / "vettore" / "static")
    # A source checkout serves the live frontend build before any stale packaged copy.
    candidates.append(Path(__file__).resolve().parents[2] / "frontend" / "dist")
    candidates.append(Path(__file__).parent / "static")
    return next((c for c in candidates if (c / "index.html").is_file()), None)


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_server(port: int, token: str, startup_files: list[Path] | None = None,
                 on_unsaved=None, on_hand_over=None):
    import uvicorn

    from .server import TOKEN_HEADER, create_app

    app = create_app(token, static_dir(), startup_files or [], on_unsaved, on_hand_over)
    # log_config=None: uvicorn's default formatter calls isatty() on streams that do not
    # exist in windowed builds; our own logging setup is used instead.
    config = uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning",
                            log_config=None)
    server = uvicorn.Server(config)
    threading.Thread(target=server.run, daemon=True).start()
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline:
        try:
            req = urllib.request.Request(f"http://127.0.0.1:{port}/api/health",
                                         headers={TOKEN_HEADER: token})
            urllib.request.urlopen(req, timeout=1).read()
            return server
        except OSError:
            time.sleep(0.1)
    raise RuntimeError("The local server did not start")


class DesktopApi:
    """Methods exposed to JavaScript as ``window.pywebview.api``."""

    def __init__(self) -> None:
        # Private on purpose: pywebview recursively exposes every public attribute of the
        # js_api object to JavaScript; walking the native window object never terminates on
        # Windows (Bounds.Empty.Empty...) and blocks the UI thread.
        self._window = None
        self._guard: close_guard.CloseGuard | None = None
        self._presence: WindowPresence | None = None
        self._unsaved = 0

    def _attach(self, window) -> None:
        """Called once the native window exists: guard its closing, follow its state."""
        self._window = window
        self._guard = close_guard.CloseGuard(window)
        self._guard.set_unsaved(self._unsaved)
        window.events.closing += self._guard.on_closing
        self._presence = WindowPresence(window, TITLE)

    def set_unsaved(self, count: int) -> None:
        """How many documents have unsaved changes (the page reports it through the API)."""
        self._unsaved = count
        if self._guard is not None:
            self._guard.set_unsaved(count)

    def confirm_quit(self) -> None:
        """The page saved or discarded every unsaved document: close the window."""
        if self._guard is not None:
            self._guard.confirm()

    def _handed_over(self) -> None:
        """A second start handed files over: show the window and tell the page at once
        (it would also find them at its next look, within half a second)."""
        window = self._window
        if window is None:
            return

        def show() -> None:
            if self._presence is not None:
                self._presence.bring_to_front()
            try:
                window.run_js(instance.files_script())
            except Exception:
                log.debug("Could not tell the page about handed-over files", exc_info=True)

        # Off the server's thread: the page may take a moment to answer.
        threading.Thread(target=show, name="vettore-hand-over", daemon=True).start()

    def open_dialog(self) -> list[str]:
        import webview

        paths = self._window.create_file_dialog(
            webview.FileDialog.OPEN, directory=settings.last_folder(), allow_multiple=True,
            file_types=file_types(),
        )
        return list(paths or [])

    def save_dialog(self, suggested_name: str, directory: str = "",
                    kind: str = "genbank") -> str | None:
        """Ask where to save. ``kind`` (genbank, fasta, dna, svg, png) orders the filters; the
        dialog starts in ``directory`` or, failing that, the folder used last."""
        import webview

        start = directory if directory and Path(directory).is_dir() else settings.last_folder()
        path = self._window.create_file_dialog(
            webview.FileDialog.SAVE, directory=start, save_filename=suggested_name,
            file_types=save_file_types(kind),
        )
        if not path:
            return None
        return path if isinstance(path, str) else path[0]


def _ensure_streams() -> tuple[Path, bool]:
    """Windowed (no-console) builds have no stdout/stderr; send them to the log file. Returns
    the log file and whether stderr now writes to it."""
    log_path = settings.log_path()
    if sys.stdout is None or sys.stderr is None:
        stream = open(log_path, "a", encoding="utf-8", buffering=1)  # noqa: SIM115
        sys.stdout = sys.stdout or stream
        sys.stderr = sys.stderr or stream
        return log_path, True
    return log_path, False


def _parse(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="vettore")
    parser.add_argument("--browser", action="store_true",
                        help="open the UI in the system browser instead of a native window")
    parser.add_argument("--no-window", action="store_true",
                        help="only run the server (for tests); prints the URL")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--token", default=None, help=argparse.SUPPRESS)
    # A server started with --no-window always runs on its own; this makes it the running
    # copy that later starts hand their files to (used by the tests of that hand-over).
    parser.add_argument("--single-instance", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("files", nargs="*", help="sequence files to open")
    return parser.parse_args(argv)


def _hand_over(files: list[Path]) -> bool:
    """Give the files to a running copy of Vettore; True when it took them."""
    running = instance.find_live()
    if running is None:
        return False
    if instance.hand_over(running, files):
        log.info("Handed %d file(s) to the running Vettore (pid %d)", len(files), running.pid)
        return True
    log.warning("The running Vettore (pid %d) did not open the files; starting a new window",
                running.pid)
    return False


def main(argv: list[str] | None = None) -> None:
    log_path, redirected = _ensure_streams()
    args = _parse(argv)

    # Messages for unexpected problems tell the user the details are in the log file, so it
    # is written whether or not there is a console.
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if not redirected:
        try:
            handlers.append(logging.FileHandler(log_path, encoding="utf-8"))
        except OSError:
            pass
    logging.basicConfig(level=logging.INFO, handlers=handlers,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    files = [Path(f).resolve() for f in args.files]
    single = args.single_instance or not args.no_window
    api = DesktopApi()
    record = None
    # Starts that overlap (several files opened at once) take turns: the first becomes the
    # running copy, the next ones find it and hand their files over.
    with instance.startup_lock() if single else contextlib.nullcontext():
        if single and _hand_over(files):
            return
        token = args.token or secrets.token_urlsafe(24)
        port = args.port or free_port()
        start_server(port, token, files, api.set_unsaved, api._handed_over)
        if single:
            record = instance.register(port, token)
    try:
        _run_ui(api, f"http://127.0.0.1:{port}/?token={token}", args)
    finally:
        if record is not None:
            instance.unregister(record)


def _run_ui(api: DesktopApi, url: str, args: argparse.Namespace) -> None:
    if args.no_window:
        print(url, flush=True)
        threading.Event().wait()
        return

    if not args.browser:
        try:
            import webview

            window = webview.create_window(
                TITLE, url, js_api=api, width=1400, height=900, maximized=True,
                # Shown only when the page does not answer (desktop.close_guard): normally
                # Vettore's own dialog offers to save each document.
                localization={"global.quitConfirmation": close_guard.FALLBACK_QUESTION},
            )
            api._attach(window)
            window.events.loaded += lambda: _bind_drop(window)
            storage = Path.home() / ".vettore" / "webview"
            storage.mkdir(parents=True, exist_ok=True)
            # On Linux the bundled Qt WebEngine is used (GTK bindings are not shipped).
            gui = "qt" if sys.platform.startswith("linux") else None
            # The window icon option is honoured by the Qt backend only; on Windows the window
            # takes the executable's icon (packaging/vettore.spec).
            static = static_dir()
            icon = str(static / "vettore.png") if gui == "qt" and static else None
            webview.start(gui=gui, private_mode=False, storage_path=str(storage), icon=icon)
            return
        except Exception:  # no GUI backend available: fall back to the browser
            log.exception("Native window unavailable, opening the system browser")

    webbrowser.open(url)
    print(f"Vettore is running at {url} — close this window to quit.", flush=True)
    try:
        threading.Event().wait()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
