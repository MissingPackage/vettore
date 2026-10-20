"""The desktop behaviours in the real Qt WebEngine window (Linux), under a virtual display.

Opt-in: these start the full application with its window, which needs Qt WebEngine's system
libraries and takes tens of seconds. Run them with

    VETTORE_QT_TESTS=1 uv run pytest tests/test_desktop_qt.py

with a display, or without one when ``xvfb-run`` is installed. A driver thread inside the app
process plays the user: it makes a document unsaved through the page, closes the window the
way the close button does (``window.destroy()`` reaches Qt's ``closeEvent``, as the button
does), answers dialogs, and starts a second ``vettore <file>`` for the hand-over.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"

pytestmark = pytest.mark.skipif(
    os.environ.get("VETTORE_QT_TESTS") != "1" or not sys.platform.startswith("linux"),
    reason="opt-in (VETTORE_QT_TESTS=1): starts the real Qt window",
)

DRIVER = r'''
import json, os, subprocess, sys, threading, time

from qtpy.QtCore import QObject, Signal
from qtpy.QtWidgets import QApplication, QMessageBox

import webview

scenario = sys.argv[1]


def out(*parts):
    print("DRIVER", *parts, flush=True)


class OnMain(QObject):
    """Runs a function on the Qt thread (created here, on the main thread)."""
    call = Signal(object)


on_main = OnMain()
on_main.call.connect(lambda f: f())


def modal_text(answer=False):
    """Text of the system dialog on screen, if any; answers Yes when asked to."""
    found, done = [], threading.Event()

    def look():
        box = QApplication.activeModalWidget()
        if isinstance(box, QMessageBox):
            found.append(box.text())
            if answer:
                box.button(QMessageBox.Yes).click()
        done.set()

    on_main.call.emit(look)
    done.wait(10)
    return found[0] if found else None


def wait_for(window, script, timeout=90):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            if window.evaluate_js(script):
                return
        except Exception:
            pass
        time.sleep(0.2)
    out("TIMEOUT", script)
    os._exit(3)


def q(testid):
    return f"document.querySelector('[data-testid={testid}]')"


def make_unsaved(window):
    wait_for(window, f"!!{q('new')}")
    window.evaluate_js(f"{q('new')}.click()")
    wait_for(window, f"!!{q('new-seq-text')}")
    window.evaluate_js(
        "(() => { const t = " + q("new-seq-text") + ";"
        " Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value')"
        ".set.call(t, 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACC');"
        " t.dispatchEvent(new Event('input', {bubbles: true})); return true; })()")
    window.evaluate_js(f"{q('new-seq-create')}.click()")
    wait_for(window, "!!document.querySelector('[data-testid=doctab][data-unsaved=true]')")
    deadline = time.monotonic() + 30
    while window._js_api._guard.unsaved != 1:
        if time.monotonic() > deadline:
            out("TIMEOUT unsaved count")
            os._exit(3)
        time.sleep(0.1)
    out("UNSAVED", window._js_api._guard.unsaved)


def drive():
    while not webview.windows:
        time.sleep(0.1)
    window = webview.windows[0]
    window.events.loaded.wait(90)
    make_unsaved(window)
    if scenario == "dialog":
        window.destroy()                      # Qt closeEvent, as the close button
        wait_for(window, f"!!{q('quit-dialog')}")
        out("DIALOG", window.evaluate_js(f"{q('quit-dialog')}.innerText").replace("\n", " | "))
        window.evaluate_js(f"{q('quit-cancel')}.click()")
        wait_for(window, f"!{q('quit-dialog')}")
        time.sleep(5)                         # longer than the fallback's wait
        out("AFTER_CANCEL", "open" if window in webview.windows else "closed",
            "system-dialog" if modal_text() else "no-system-dialog")
        window.destroy()
        wait_for(window, f"!!{q('quit-dialog')}")
        out("DIALOG_AGAIN")
        window.evaluate_js(f"{q('quit-discard')}.click()")
    elif scenario == "frozen":
        window.evaluate_js("setTimeout(() => { const end = Date.now() + 120000; "
                           "while (Date.now() < end) {} }, 100); true")
        time.sleep(1)
        started = time.monotonic()
        window.destroy()
        while time.monotonic() - started < 60:
            text = modal_text(answer=True)
            if text:
                out("SYSTEM_DIALOG", f"{time.monotonic() - started:.1f}", text)
                break
            time.sleep(0.3)
    elif scenario == "handover":
        path = sys.argv[2]
        second = subprocess.run([sys.executable, "-m", "vettore.app", path], timeout=120)
        out("SECOND_EXIT", second.returncode)
        wait_for(window, "[...document.querySelectorAll('[data-testid=doctab]')]"
                         f".some(t => t.title === {json.dumps(path)})")
        out("TAB_OPENED")
        window.destroy()
        wait_for(window, f"!!{q('quit-dialog')}")
        window.evaluate_js(f"{q('quit-discard')}.click()")
    elif scenario == "print":
        # What File → Print… does in the end: print() on a frame holding the map.
        window.evaluate_js("(() => { const f = document.createElement('iframe');"
                           " document.body.appendChild(f);"
                           " setTimeout(() => f.contentWindow.print(), 10); return true; })()")
        wait_for(window, "[...document.querySelectorAll('[data-testid=notice]')]"
                         ".some(n => n.innerText.includes('cannot print on Linux'))")
        out("PRINT_EXPLAINED")
        window.destroy()
        wait_for(window, f"!!{q('quit-dialog')}")
        window.evaluate_js(f"{q('quit-discard')}.click()")


threading.Thread(target=drive, daemon=True).start()
from vettore.app import main

main([])
out("EXITED")
'''


def run(scenario: str, tmp_path: Path, *args: str) -> list[str]:
    env = {**os.environ, "VETTORE_HOME": str(tmp_path / "home"),
           "QTWEBENGINE_DISABLE_SANDBOX": "1", "QT_QPA_PLATFORM": "xcb"}
    cmd = [sys.executable, "-c", DRIVER, scenario, *args]
    if not os.environ.get("DISPLAY"):
        if not shutil.which("xvfb-run"):
            pytest.skip("no display and no xvfb-run")
        cmd = ["xvfb-run", "-a", "-s", "-screen 0 1400x900x24", *cmd]
    done = subprocess.run(cmd, check=False, env=env, capture_output=True, text=True,
                          timeout=300)
    lines = [line.removeprefix("DRIVER ") for line in done.stdout.splitlines()
             if line.startswith("DRIVER ")]
    print(scenario, lines)
    assert done.returncode == 0, (lines, done.stderr[-3000:])
    # A copy that closed normally no longer claims to be the running one.
    assert not (tmp_path / "home" / "instance.json").exists()
    return lines


def test_closing_the_window_asks_in_vettores_dialog(tmp_path):
    lines = run("dialog", tmp_path)
    assert lines[0] == "UNSAVED 1"
    assert lines[1].startswith("DIALOG ")
    assert "“New DNA” has changes that are not saved." in lines[1]
    assert lines[2] == "AFTER_CANCEL open no-system-dialog"
    assert lines[3:] == ["DIALOG_AGAIN", "EXITED"]


def test_a_frozen_page_gets_the_system_dialog(tmp_path):
    lines = run("frozen", tmp_path)
    assert lines[0] == "UNSAVED 1"
    assert lines[1].startswith("SYSTEM_DIALOG ")
    seconds = float(lines[1].split()[1])
    assert 3.0 <= seconds < 15, seconds
    assert "not responding" in lines[1]
    assert lines[-1] == "EXITED"


def test_a_second_start_opens_its_file_in_the_qt_window(tmp_path):
    sample = tmp_path / "from-file-manager.gb"
    sample.write_bytes((DATA / "MT612431.1.gb").read_bytes())
    lines = run("handover", tmp_path, str(sample))
    assert lines == ["UNSAVED 1", "SECOND_EXIT 0", "TAB_OPENED", "EXITED"]


def test_printing_in_the_qt_window_explains_what_to_do_instead(tmp_path):
    """Qt WebEngine drops print requests that the application does not handle; without the
    message, File → Print… did nothing at all."""
    lines = run("print", tmp_path)
    assert lines == ["UNSAVED 1", "PRINT_EXPLAINED", "EXITED"]
