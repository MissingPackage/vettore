"""Bring the window to the front when files are handed over from a second start.

Best effort, never an error: the files open either way. On Windows the window is found by
process id and restored and activated through Win32 (``ShowWindow`` restores a minimized
window to maximized if it was maximized; the second start has allowed this process to take
the foreground, see ``instance.allow_foreground``). Elsewhere pywebview un-minimizes it; a
window that is open but behind others stays where it is, as window managers decide.
"""

from __future__ import annotations

import logging
import os
import sys
from typing import Any

log = logging.getLogger(__name__)


class WindowPresence:
    """Follows whether the window is minimized, and brings it back to the front."""

    def __init__(self, window: Any, title: str) -> None:
        self._window = window
        self._title = title
        self._minimized = False
        self._maximized = True     # the window opens maximized (vettore.app)
        events = window.events
        events.minimized += self._on_minimized
        events.maximized += self._on_maximized
        events.restored += self._on_restored

    def _on_minimized(self) -> None:
        self._minimized = True

    def _on_maximized(self) -> None:
        self._minimized, self._maximized = False, True

    def _on_restored(self) -> None:
        # pywebview reports "restored" for any return to the normal state.
        self._minimized, self._maximized = False, False

    def bring_to_front(self) -> None:
        try:
            if sys.platform == "win32":
                _win32_bring_to_front(os.getpid(), self._title)
            elif self._minimized:
                if self._maximized:
                    self._window.maximize()
                else:
                    self._window.restore()
        except Exception:
            log.debug("Could not bring the window to the front", exc_info=True)


def _win32_bring_to_front(pid: int, title: str) -> None:
    import ctypes
    from ctypes import wintypes

    user32 = ctypes.WinDLL("user32", use_last_error=True)
    enum_proc = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)
    # Declared, so window handles above 2**31 are passed as handles, not as C ints.
    user32.EnumWindows.argtypes = (enum_proc, wintypes.LPARAM)
    user32.GetWindowThreadProcessId.argtypes = (wintypes.HWND, ctypes.POINTER(wintypes.DWORD))
    user32.GetWindowThreadProcessId.restype = wintypes.DWORD
    user32.GetWindowTextW.argtypes = (wintypes.HWND, wintypes.LPWSTR, ctypes.c_int)
    for name in ("IsWindowVisible", "IsIconic", "SetForegroundWindow"):
        getattr(user32, name).argtypes = (wintypes.HWND,)
    user32.ShowWindow.argtypes = (wintypes.HWND, ctypes.c_int)
    found: list[int] = []

    @enum_proc
    def visit(hwnd, _lparam):
        owner = wintypes.DWORD()
        user32.GetWindowThreadProcessId(hwnd, ctypes.byref(owner))
        if owner.value == pid and user32.IsWindowVisible(hwnd):
            buffer = ctypes.create_unicode_buffer(256)
            user32.GetWindowTextW(hwnd, buffer, 256)
            if buffer.value == title:
                found.append(hwnd)
                return False
        return True

    user32.EnumWindows(visit, 0)
    for hwnd in found:
        if user32.IsIconic(hwnd):
            user32.ShowWindow(hwnd, 9)   # SW_RESTORE: back to its size before minimizing
        user32.SetForegroundWindow(hwnd)
