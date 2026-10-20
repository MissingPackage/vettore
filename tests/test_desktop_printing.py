"""Print requests the Linux window (Qt WebEngine) would drop reach the page as a message.

The real window is checked in tests/test_desktop_qt.py (opt-in); here Qt's page is a fake with
the two signals Qt WebEngine emits for ``window.print()``.
"""

import threading

from vettore.desktop.printing import PRINT_EVENT, print_script, watch_print_requests


class Signal:
    def __init__(self) -> None:
        self.slots = []

    def connect(self, slot) -> None:
        self.slots.append(slot)

    def emit(self, *args) -> None:
        for slot in self.slots:
            slot(*args)


class Page:
    def __init__(self, by_frame: bool = True) -> None:
        self.printRequested = Signal()
        if by_frame:
            self.printRequestedByFrame = Signal()


class Window:
    def __init__(self, page) -> None:
        self.native = type("BrowserView", (), {})()
        self.native.webview = type("WebView", (), {"page": lambda _self: page})()
        self.scripts: list[str] = []
        self.ran = threading.Event()

    def run_js(self, script) -> None:
        self.scripts.append(script)
        self.ran.set()


def test_a_print_request_from_any_frame_tells_the_page():
    page = Page()
    window = Window(page)
    assert watch_print_requests(window) is True
    page.printRequestedByFrame.emit(object())   # the map is printed from an iframe
    assert window.ran.wait(2)
    assert window.scripts == [print_script()]
    window.ran.clear()
    page.printRequested.emit()
    assert window.ran.wait(2)
    assert len(window.scripts) == 2


def test_an_older_qt_without_frame_requests_still_reports_the_main_frame():
    page = Page(by_frame=False)
    window = Window(page)
    assert watch_print_requests(window) is True
    page.printRequested.emit()
    assert window.ran.wait(2)


def test_a_window_without_qt_is_left_alone():
    window = type("Other", (), {"native": None})()
    assert watch_print_requests(window) is False


def test_the_page_event():
    # frontend/src/files/desktop.ts listens for this event.
    assert PRINT_EVENT == "vettore:print-unavailable"
    assert print_script() == 'window.dispatchEvent(new CustomEvent("vettore:print-unavailable"))'
