"""Printing in the Linux window (Qt WebEngine).

File → Print… prints through the page (``window.print()`` on a frame holding the map). Edge
WebView2, the Windows window, answers that with its own print preview, whose printers include
Save as PDF. Qt WebEngine instead emits ``QWebEnginePage.printRequested`` (main frame) or
``printRequestedByFrame`` (Qt 6.8+, any frame) and prints nothing unless the application
handles the signal; pywebview 6.2.1 does not. Until Vettore prints there itself, the page is
told, so Print… explains what to do instead of silently doing nothing.
"""

from __future__ import annotations

import json
import logging
import threading
from typing import Any

log = logging.getLogger(__name__)

#: Dispatched on the page when a print request cannot be served (frontend/src/App.tsx).
PRINT_EVENT = "vettore:print-unavailable"

SIGNALS = ("printRequested", "printRequestedByFrame")


def print_script() -> str:
    return f"window.dispatchEvent(new CustomEvent({json.dumps(PRINT_EVENT)}))"


def watch_print_requests(window: Any) -> bool:
    """Tell the page about print requests Qt WebEngine would drop. Must run on the Qt thread
    (a ``before_show`` handler does), so the signals are delivered there."""
    try:
        page = window.native.webview.page()
    except AttributeError:
        return False

    def tell(*_args: Any) -> None:
        # On the Qt thread: run_js waits for that thread, so it is called from another one.
        threading.Thread(target=lambda: window.run_js(print_script()), name="vettore-print",
                         daemon=True).start()

    connected = False
    for name in SIGNALS:
        signal = getattr(page, name, None)
        if signal is not None:
            signal.connect(tell)
            connected = True
    if not connected:
        log.warning("This Qt WebEngine has no print request signal; Print… does nothing")
    return connected
