"""Closing the window with unsaved work asks in Vettore's own dialog.

pywebview runs the window's ``closing`` handlers on the UI thread before the window closes;
a handler that returns False keeps the window open (Qt: ``closeEvent`` ignores the event,
WinForms: ``FormClosing`` sets ``Cancel``). Calling into the page from that handler would
deadlock, because the page answers on the same UI thread. So the guard:

1. keeps the window open and, from a background thread, asks the page to show its "Unsaved
   changes" dialog (a cancelable ``QUIT_EVENT`` on ``window``; the page cancels it to say it
   took over);
2. closes the window when the page calls :meth:`CloseGuard.confirm` (the documents were saved,
   or the user chose not to save them);
3. falls back to pywebview's native confirm (``confirm_close``) when the page does not answer
   within ``answer_timeout`` seconds, so a frozen page can still be closed.

``window.destroy()`` closes through the same path as the window's close button on both
backends, so the guard sees every close.
"""

from __future__ import annotations

import json
import logging
import threading
from typing import Any

log = logging.getLogger(__name__)

#: Dispatched on the page's ``window``; cancelable, the page cancels it when it will answer.
QUIT_EVENT = "vettore:quit-requested"

#: How long the page has to take over before the native confirm asks instead (seconds).
ANSWER_TIMEOUT = 3.0

#: Text of pywebview's native confirm, shown only when the page did not answer.
FALLBACK_QUESTION = ("Vettore's window is not responding, so it cannot offer to save your work. "
                     "Some documents have unsaved changes: quit Vettore and discard them?")


def quit_script() -> str:
    """JavaScript that asks the page to take over closing; evaluates to true when it does."""
    return (f"!window.dispatchEvent(new CustomEvent({json.dumps(QUIT_EVENT)}, "
            "{cancelable: true}))")


class CloseGuard:
    """Decides, for one window, whether a close goes ahead. Thread-safe."""

    def __init__(self, window: Any, answer_timeout: float = ANSWER_TIMEOUT) -> None:
        self._window = window
        self._timeout = answer_timeout
        self._lock = threading.Lock()
        self._unsaved = 0
        self._confirmed = False   # the page said: close now
        self._fallback = False    # the page did not answer: the native confirm decides
        self._asking = False      # a request to the page is waiting for its answer

    @property
    def unsaved(self) -> int:
        return self._unsaved

    def set_unsaved(self, count: int) -> None:
        """The page reports how many documents have unsaved changes."""
        with self._lock:
            self._unsaved = max(0, int(count))
        # Read by pywebview only after on_closing lets a close through: that is, in the
        # fallback, where it makes the native confirm ask before discarding anything.
        self._window.confirm_close = self._unsaved > 0

    def on_closing(self) -> bool:
        """``window.events.closing`` handler. Runs on the UI thread: never blocks it.
        Returns False to keep the window open."""
        with self._lock:
            if self._confirmed or self._unsaved == 0:
                return True
            if self._fallback:
                self._fallback = False   # this close only; the next one asks the page again
                return True
            if self._asking:             # a second click while the page is being asked
                return False
            self._asking = True
        threading.Thread(target=self._ask_page, name="vettore-close", daemon=True).start()
        return False

    def confirm(self) -> None:
        """The page saved or discarded every document: close without asking again."""
        with self._lock:
            self._confirmed = True
        self._window.confirm_close = False
        self._window.destroy()

    def _ask_page(self) -> None:
        answered = threading.Event()
        reply: dict[str, Any] = {}

        def ask() -> None:
            try:
                reply["value"] = self._window.evaluate_js(quit_script())
            except Exception:  # a page that failed to load, or a window being torn down
                log.exception("Could not ask the page about unsaved documents")
            finally:
                answered.set()

        # evaluate_js waits for the page; a frozen page never answers, so wait separately.
        threading.Thread(target=ask, name="vettore-close-ask", daemon=True).start()
        took_over = answered.wait(self._timeout) and reply.get("value") is True
        with self._lock:
            self._asking = False
            if took_over:
                return
            self._fallback = True
        log.warning("The page did not answer the close request within %.1f s; asking with the "
                    "system dialog instead", self._timeout)
        self._window.destroy()   # closes again; on_closing lets the native confirm ask
