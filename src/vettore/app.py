"""Desktop entry point: starts the local API server and opens the UI in a native window."""

from __future__ import annotations

import argparse
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

import uvicorn

from . import app_info, fileio, settings
from .server import TOKEN_HEADER, create_app

log = logging.getLogger("vettore")

# Open dialog: every extension the engine reads (it also recognises files by their content, so
# "All files" opens files with a wrong or missing extension). "*" matches names without a dot.
FILE_TYPES = (
    "Sequence files (" + ";".join(f"*{ext}" for ext in fileio.FORMATS) + ")",
    "All files (*.*;*)",
)
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
                 on_unsaved=None) -> uvicorn.Server:
    app = create_app(token, static_dir(), startup_files or [], on_unsaved)
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

    def set_unsaved(self, count: int) -> None:
        # pywebview reads confirm_close when the window is about to close.
        if self._window is not None:
            self._window.confirm_close = count > 0

    def open_dialog(self) -> list[str]:
        import webview

        paths = self._window.create_file_dialog(
            webview.FileDialog.OPEN, directory=settings.last_folder(), allow_multiple=True,
            file_types=FILE_TYPES,
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
    log_path = app_info.log_path()
    if sys.stdout is None or sys.stderr is None:
        stream = open(log_path, "a", encoding="utf-8", buffering=1)  # noqa: SIM115
        sys.stdout = sys.stdout or stream
        sys.stderr = sys.stderr or stream
        return log_path, True
    return log_path, False


def main(argv: list[str] | None = None) -> None:
    log_path, redirected = _ensure_streams()
    parser = argparse.ArgumentParser(prog="vettore")
    parser.add_argument("--browser", action="store_true",
                        help="open the UI in the system browser instead of a native window")
    parser.add_argument("--no-window", action="store_true",
                        help="only run the server (for tests); prints the URL")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--token", default=None, help=argparse.SUPPRESS)
    parser.add_argument("files", nargs="*", help="sequence files to open")
    args = parser.parse_args(argv)

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
    token = args.token or secrets.token_urlsafe(24)
    port = args.port or free_port()
    api = DesktopApi()
    start_server(port, token, [Path(f).resolve() for f in args.files], api.set_unsaved)
    url = f"http://127.0.0.1:{port}/?token={token}"

    if args.no_window:
        print(url, flush=True)
        threading.Event().wait()
        return

    if not args.browser:
        try:
            import webview

            api._window = webview.create_window(
                "Vettore", url, js_api=api, width=1400, height=900, maximized=True,
                localization={"global.quitConfirmation":
                              "Some documents have unsaved changes. Quit Vettore and discard them?"},
            )
            window = api._window
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
