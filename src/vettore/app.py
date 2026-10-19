"""Desktop entry point: starts the local API server and opens the UI in a native window."""

from __future__ import annotations

import argparse
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

import uvicorn

from .server import TOKEN_HEADER, create_app

log = logging.getLogger("vettore")

FILE_TYPES = (
    "Sequence files (*.dna;*.gb;*.gbk;*.gbff;*.ape;*.fa;*.fasta;*.fna;*.ab1)",
    "All files (*.*)",
)


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


def start_server(port: int, token: str, startup_files: list[Path] | None = None) -> uvicorn.Server:
    app = create_app(token, static_dir(), startup_files or [])
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
        self.window = None

    def open_dialog(self) -> list[str]:
        import webview

        paths = self.window.create_file_dialog(
            webview.FileDialog.OPEN, allow_multiple=True, file_types=FILE_TYPES
        )
        return list(paths or [])

    def save_dialog(self, suggested_name: str) -> str | None:
        import webview

        path = self.window.create_file_dialog(
            webview.FileDialog.SAVE, save_filename=suggested_name,
            file_types=("GenBank (*.gb)", "FASTA (*.fasta)"),
        )
        if not path:
            return None
        return path if isinstance(path, str) else path[0]


def _ensure_streams() -> Path:
    """Windowed (no-console) builds have no stdout/stderr; send them to a log file."""
    home = Path.home() / ".vettore"
    home.mkdir(parents=True, exist_ok=True)
    log_path = home / "vettore.log"
    if sys.stdout is None or sys.stderr is None:
        stream = open(log_path, "a", encoding="utf-8", buffering=1)  # noqa: SIM115
        sys.stdout = sys.stdout or stream
        sys.stderr = sys.stderr or stream
    return log_path


def main(argv: list[str] | None = None) -> None:
    _ensure_streams()
    parser = argparse.ArgumentParser(prog="vettore")
    parser.add_argument("--browser", action="store_true",
                        help="open the UI in the system browser instead of a native window")
    parser.add_argument("--no-window", action="store_true",
                        help="only run the server (for tests); prints the URL")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--token", default=None, help=argparse.SUPPRESS)
    parser.add_argument("files", nargs="*", help="sequence files to open")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    token = args.token or secrets.token_urlsafe(24)
    port = args.port or free_port()
    start_server(port, token, [Path(f).resolve() for f in args.files])
    url = f"http://127.0.0.1:{port}/?token={token}"

    if args.no_window:
        print(url, flush=True)
        threading.Event().wait()
        return

    if not args.browser:
        try:
            import webview

            api = DesktopApi()
            api.window = webview.create_window("Vettore", url, js_api=api, width=1400,
                                               height=900)
            storage = Path.home() / ".vettore" / "webview"
            storage.mkdir(parents=True, exist_ok=True)
            # On Linux the bundled Qt WebEngine is used (GTK bindings are not shipped).
            gui = "qt" if sys.platform.startswith("linux") else None
            webview.start(gui=gui, private_mode=False, storage_path=str(storage))
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
