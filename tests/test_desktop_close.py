"""Closing the window with unsaved work: the page's own dialog asks, the native one is the
fallback for a page that does not answer.

The window is a fake with pywebview's interface. As on both real backends (Qt ``closeEvent``,
WinForms ``FormClosing``), ``destroy()`` closes through the closing handlers, and the native
confirm appears after them only when ``confirm_close`` is set.
"""

import threading
import time

from vettore.app import DesktopApi
from vettore.desktop.close_guard import (
    FALLBACK_QUESTION,
    QUIT_EVENT,
    CloseGuard,
    quit_script,
)

TIMEOUT = 0.3


class Event:
    def __init__(self) -> None:
        self.handlers = []

    def __iadd__(self, handler):
        self.handlers.append(handler)
        return self


class Events:
    def __init__(self) -> None:
        for name in ("closing", "loaded", "minimized", "maximized", "restored"):
            setattr(self, name, Event())


class FakeWindow:
    """``answer`` is what the page's script evaluates to; ``hang`` makes it never answer."""

    def __init__(self, answer=True, hang=False, error=None) -> None:
        self.events = Events()
        self.confirm_close = False
        self.answer, self.hang, self.error = answer, hang, error
        self.scripts: list[str] = []
        self.asked = threading.Event()
        self.release = threading.Event()   # lets a hanging page go at the end of a test
        self.closed = threading.Event()
        self.native_confirms = 0

    def evaluate_js(self, script):
        self.scripts.append(script)
        self.asked.set()
        if self.hang:
            self.release.wait(5)
        if self.error:
            raise self.error
        return self.answer

    def run_js(self, script):
        self.scripts.append(script)

    def close_button(self) -> bool:
        """What the backend does on a close: the handlers first, then its own confirm."""
        if not all(handler() for handler in self.events.closing.handlers):
            return False
        if self.confirm_close:
            self.native_confirms += 1
        self.closed.set()
        return True

    def destroy(self):
        self.close_button()


def guarded(window: FakeWindow, unsaved: int) -> CloseGuard:
    guard = CloseGuard(window, answer_timeout=TIMEOUT)
    window.events.closing += guard.on_closing
    guard.set_unsaved(unsaved)
    return guard


def test_without_unsaved_work_the_window_closes_at_once():
    window = FakeWindow()
    guarded(window, 0)
    assert window.close_button() is True
    assert window.scripts == [] and window.native_confirms == 0


def test_unsaved_work_keeps_the_window_open_and_asks_the_page():
    window = FakeWindow(answer=True)
    guard = guarded(window, 2)
    assert window.close_button() is False
    assert window.asked.wait(2)
    assert window.scripts == [quit_script()]
    time.sleep(TIMEOUT * 2)          # the page took over: no fallback follows
    assert not window.closed.is_set() and window.native_confirms == 0

    guard.confirm()                  # the page saved everything, or the user chose Don't Save
    assert window.closed.is_set()
    assert window.native_confirms == 0
    assert window.confirm_close is False


def test_a_second_click_while_the_page_is_asked_does_not_ask_twice():
    window = FakeWindow(hang=True)
    guarded(window, 1)
    assert window.close_button() is False
    assert window.asked.wait(2)
    assert window.close_button() is False
    assert len(window.scripts) == 1
    window.release.set()


def test_a_frozen_page_falls_back_to_the_native_confirm():
    window = FakeWindow(hang=True)
    guarded(window, 1)
    started = time.monotonic()
    assert window.close_button() is False
    assert window.closed.wait(5)
    assert time.monotonic() - started >= TIMEOUT
    assert window.native_confirms == 1     # confirm_close was set: the system dialog asked
    window.release.set()


def test_the_fallback_covers_one_close_only():
    """The user answered No in the native confirm: the next close asks the page again."""
    window = FakeWindow(answer=None)       # a page that does not handle the request
    guard = guarded(window, 1)
    assert window.close_button() is False
    assert window.closed.wait(5) and window.native_confirms == 1
    window.closed.clear()
    window.asked.clear()
    window.answer = True                   # the page is back
    assert window.close_button() is False
    assert window.asked.wait(2) and len(window.scripts) == 2
    time.sleep(TIMEOUT * 2)
    assert not window.closed.is_set()
    assert guard.unsaved == 1


def test_a_page_that_cannot_be_reached_falls_back_too():
    window = FakeWindow(error=RuntimeError("page not loaded"))
    guarded(window, 3)
    assert window.close_button() is False
    assert window.closed.wait(5) and window.native_confirms == 1


def test_the_native_confirm_asks_only_while_there_is_unsaved_work():
    window = FakeWindow()
    guard = guarded(window, 2)
    assert window.confirm_close is True
    guard.set_unsaved(0)
    assert window.confirm_close is False
    assert "unsaved changes" in FALLBACK_QUESTION and "not responding" in FALLBACK_QUESTION


def test_the_request_to_the_page_is_a_cancelable_event():
    # frontend/e2e/desktop-quit.spec.ts dispatches exactly this expression.
    assert QUIT_EVENT == "vettore:quit-requested"
    assert quit_script() == ('!window.dispatchEvent(new CustomEvent("vettore:quit-requested", '
                             '{cancelable: true}))')


def test_desktop_api_guards_the_window_it_is_attached_to():
    api = DesktopApi()
    api.set_unsaved(1)                      # reported before the window exists
    window = FakeWindow(answer=True)
    api._attach(window)
    assert window.confirm_close is True
    assert window.close_button() is False
    assert window.asked.wait(2)
    api.confirm_quit()
    assert window.closed.is_set() and window.native_confirms == 0
    # Nothing but methods reaches JavaScript (pywebview walks public attributes).
    assert [n for n in vars(api) if not n.startswith("_")] == []


def test_confirm_quit_without_a_window_does_nothing():
    DesktopApi().confirm_quit()   # browser mode: no window to close
