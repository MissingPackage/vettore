"""Open sequence files with Vettore on double-click: what is written, and that undo removes it.

Windows is tested against a fake ``winreg`` with the registry's behaviour (case-insensitive
key names, a key with subkeys cannot be deleted, missing keys and values raise OSError);
Linux against temporary data and config folders. The packaged smoke test checks that the
real functions are reachable without writing anything.
"""

import configparser
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from vettore import fileio
from vettore.desktop import associations
from vettore.desktop.associations import (
    EXTENSIONS,
    LINUX_TYPES,
    PROGID,
    LinuxAssociations,
    WindowsAssociations,
)
from vettore.server import TOKEN_HEADER, create_app

EXE = r"C:\Users\Lab User\Downloads\Vettore\Vettore.exe"
CLASSES = r"Software\Classes"
H = {TOKEN_HEADER: "t0ken"}


class FakeWinreg:
    """The part of ``winreg`` the associations use, over a dictionary."""

    HKEY_CURRENT_USER = "HKCU"
    KEY_READ, KEY_WRITE = 0x20019, 0x20006
    REG_NONE, REG_SZ = 0, 1

    def __init__(self, fail_on: str | None = None) -> None:
        self.keys: dict[str, dict[str, tuple[object, int]]] = {}
        self.fail_on = fail_on.lower() if fail_on else None

    class Handle:
        def __init__(self, path: str) -> None:
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    @staticmethod
    def _norm(path: str) -> str:
        return path.strip("\\").lower()

    def CreateKeyEx(self, root, path, reserved=0, access=0):
        assert root == self.HKEY_CURRENT_USER
        parts = self._norm(path).split("\\")
        for i in range(1, len(parts) + 1):
            self.keys.setdefault("\\".join(parts[:i]), {})
        return self.Handle(self._norm(path))

    def OpenKey(self, root, path, reserved=0, access=0):
        assert root == self.HKEY_CURRENT_USER
        if self._norm(path) not in self.keys:
            raise FileNotFoundError(2, "The system cannot find the file specified")
        return self.Handle(self._norm(path))

    def SetValueEx(self, key, name, reserved, kind, value):
        if self.fail_on and self.fail_on in key.path:
            raise PermissionError(5, "Access is denied")
        self.keys[key.path][name.lower()] = (value, kind)

    def QueryValueEx(self, key, name):
        try:
            return self.keys[key.path][name.lower()]
        except KeyError:
            raise FileNotFoundError(2, "The system cannot find the file specified") from None

    def DeleteValue(self, key, name):
        try:
            del self.keys[key.path][name.lower()]
        except KeyError:
            raise FileNotFoundError(2, "The system cannot find the file specified") from None

    def _children(self, path: str) -> list[str]:
        prefix = path + "\\"
        return sorted({k[len(prefix):].split("\\")[0] for k in self.keys if k.startswith(prefix)})

    def DeleteKey(self, root, path):
        path = self._norm(path)
        if self._children(path):
            raise PermissionError(5, "Access is denied")
        del self.keys[path]

    def EnumKey(self, key, index):
        children = self._children(key.path)
        if index >= len(children):
            raise OSError(259, "No more data is available")
        return children[index]

    def EnumValue(self, key, index):
        values = list(self.keys[key.path].items())
        if index >= len(values):
            raise OSError(259, "No more data is available")
        name, (value, kind) = values[index]
        return name, value, kind

    # test helpers
    def value(self, path: str, name: str = ""):
        return self.keys.get(self._norm(path), {}).get(name.lower(), (None, None))[0]

    def has(self, path: str) -> bool:
        return self._norm(path) in self.keys


def windows(reg: FakeWinreg, exe: str = EXE) -> tuple[WindowsAssociations, list[int]]:
    notified: list[int] = []
    return WindowsAssociations(exe, winreg=reg, notify=lambda: notified.append(1)), notified


def test_the_extensions_are_formats_vettore_opens():
    assert set(EXTENSIONS) <= set(fileio.FORMATS)
    assert {".dna", ".gb", ".gbk", ".gbff", ".ape", ".embl", ".fa", ".fasta", ".xdna", ".gck",
            ".ab1"} == set(EXTENSIONS)
    linux = [ext for _, extensions in LINUX_TYPES.values() for ext in extensions]
    assert sorted(linux) == sorted(EXTENSIONS)


# ---- Windows ----------------------------------------------------------------------------

def test_windows_writes_a_progid_and_points_every_extension_at_it():
    reg = FakeWinreg()
    assoc, notified = windows(reg)
    assoc.enable()
    progid = rf"{CLASSES}\{PROGID}"
    assert reg.value(progid) == "DNA sequence file"
    assert reg.value(rf"{progid}\DefaultIcon") == f'"{EXE}",0'
    assert reg.value(rf"{progid}\shell\open\command") == f'"{EXE}" "%1"'
    for ext in EXTENSIONS:
        assert reg.value(rf"{CLASSES}\{ext}") == PROGID, ext
        assert reg.value(rf"{CLASSES}\{ext}\OpenWithProgids", PROGID) == "", ext
    assert all(not k.startswith("software\\microsoft") for k in reg.keys)   # only Classes
    assert notified == [1]                     # Explorer is told
    assert assoc.status() == {
        "available": True, "enabled": True, "current": True, "registered_executable": EXE,
        "extensions": list(EXTENSIONS), "overridden": [], "platform": "win32"}


def test_windows_undo_restores_what_was_there_and_removes_the_rest():
    reg = FakeWinreg()
    with reg.CreateKeyEx(reg.HKEY_CURRENT_USER, rf"{CLASSES}\.fasta") as key:
        reg.SetValueEx(key, "", 0, reg.REG_SZ, "OtherEditor.fasta")
        reg.SetValueEx(key, "Content Type", 0, reg.REG_SZ, "text/plain")
    assoc, notified = windows(reg)
    assoc.enable()
    assert reg.value(rf"{CLASSES}\.fasta") == PROGID
    assoc.enable()                             # twice: the previous default is not lost
    assoc.disable()
    assert reg.value(rf"{CLASSES}\.fasta") == "OtherEditor.fasta"
    assert reg.value(rf"{CLASSES}\.fasta", "Content Type") == "text/plain"
    assert not reg.has(rf"{CLASSES}\.fasta\OpenWithProgids")
    assert not reg.has(rf"{CLASSES}\{PROGID}")
    assert not reg.has(rf"{CLASSES}\.dna")     # created by Vettore, gone again
    assert notified == [1, 1, 1]
    assert assoc.status()["enabled"] is False
    assert assoc.status()["registered_executable"] is None


def test_windows_undo_leaves_an_extension_another_program_took_since():
    reg = FakeWinreg()
    assoc, _ = windows(reg)
    assoc.enable()
    with reg.CreateKeyEx(reg.HKEY_CURRENT_USER, rf"{CLASSES}\.gb") as key:
        reg.SetValueEx(key, "", 0, reg.REG_SZ, "Newer.Program")
    assoc.disable()
    assert reg.value(rf"{CLASSES}\.gb") == "Newer.Program"


def test_windows_status_after_the_folder_was_moved():
    reg = FakeWinreg()
    windows(reg, r"D:\old\Vettore\Vettore.exe")[0].enable()
    status = windows(reg)[0].status()
    assert (status["enabled"], status["current"]) == (False, False)
    assert status["registered_executable"] == r"D:\old\Vettore\Vettore.exe"
    windows(reg)[0].enable()                   # this copy takes over
    assert windows(reg)[0].status()["enabled"] is True


def test_windows_reports_extensions_the_user_gave_to_another_app():
    """"Always use this app" is stored by Windows in UserChoice, which programs cannot set."""
    reg = FakeWinreg()
    explorer = r"Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts"
    for ext, progid in ((".dna", "OtherEditor.dna"), (".gb", r"Applications\Vettore.exe"),
                        (".ab1", PROGID)):
        with reg.CreateKeyEx(reg.HKEY_CURRENT_USER, rf"{explorer}\{ext}\UserChoice") as key:
            reg.SetValueEx(key, "ProgId", 0, reg.REG_SZ, progid)
    with reg.CreateKeyEx(reg.HKEY_CURRENT_USER, rf"{explorer}\.fa\UserChoiceLatest") as key:
        reg.SetValueEx(key, "ProgId", 0, reg.REG_SZ, "Notepad.fa")
    assoc, _ = windows(reg)
    assoc.enable()
    assert assoc.status()["overridden"] == [".dna", ".fa"]


def test_windows_a_refused_write_is_undone_and_reported():
    reg = FakeWinreg(fail_on=r"classes\.embl")
    assoc, notified = windows(reg)
    with pytest.raises(PermissionError):
        assoc.enable()
    assert not reg.has(rf"{CLASSES}\{PROGID}")
    assert all(reg.value(rf"{CLASSES}\{ext}") is None for ext in EXTENSIONS)
    assert notified


def test_windows_plan_lists_the_writes_without_writing():
    reg = FakeWinreg()
    assoc, notified = windows(reg)
    plan = assoc.plan()
    assert rf'{CLASSES}\{PROGID}\shell\open\command = "{EXE}" "%1"' in plan
    for ext in EXTENSIONS:
        assert rf"{CLASSES}\{ext} = {PROGID}" in plan
    assert reg.keys == {} and notified == []


# ---- Linux ------------------------------------------------------------------------------

LINUX_EXE = "/home/lab user/Vettore $new/Vettore"


def linux(tmp_path: Path, tools: set[str] = frozenset()) -> tuple[LinuxAssociations, list]:
    calls: list[list[str]] = []
    assoc = LinuxAssociations(LINUX_EXE, icon=tmp_path / "vettore.png",
                              data_home=tmp_path / "share", config_home=tmp_path / "config",
                              which=lambda tool: f"/usr/bin/{tool}" if tool in tools else None,
                              run=lambda cmd, **kw: calls.append(cmd))
    return assoc, calls


def test_linux_writes_a_desktop_entry_a_mime_package_and_the_defaults(tmp_path):
    mimeapps = tmp_path / "config" / "mimeapps.list"
    mimeapps.parent.mkdir(parents=True)
    mimeapps.write_text("[Default Applications]\ntext/plain=org.kde.kate.desktop;\n",
                        encoding="utf-8")
    assoc, calls = linux(tmp_path)
    assoc.enable()

    entry = (tmp_path / "share" / "applications" / "vettore.desktop").read_text(encoding="utf-8")
    assert 'Exec="/home/lab user/Vettore \\\\$new/Vettore" %F' in entry
    assert f"Icon={tmp_path / 'vettore.png'}" in entry
    assert "MimeType=" + "".join(f"{t};" for t in LINUX_TYPES) in entry

    package = ET.parse(tmp_path / "share" / "mime" / "packages" / "vettore-sequences.xml")
    ns = {"m": "http://www.freedesktop.org/standards/shared-mime-info"}
    globs = {g.get("pattern"): g.get("weight") for g in package.findall(".//m:glob", ns)}
    assert globs == {f"*{ext}": "80" for ext in EXTENSIONS}

    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read(mimeapps, encoding="utf-8")
    defaults = dict(parser.items("Default Applications"))
    assert defaults["text/plain"] == "org.kde.kate.desktop;"      # the user's own, kept
    assert all(defaults[t] == "vettore.desktop;" for t in LINUX_TYPES)
    assert calls == []                         # no tools installed: nothing run
    assert assoc.status() == {
        "available": True, "enabled": True, "current": True,
        "registered_executable": LINUX_EXE, "extensions": list(EXTENSIONS), "overridden": [],
        "platform": "linux"}


def test_linux_uses_the_desktop_tools_when_installed(tmp_path):
    assoc, calls = linux(tmp_path, {"xdg-mime", "update-mime-database",
                                    "update-desktop-database"})
    assoc.enable()
    assert calls == [
        ["update-mime-database", str(tmp_path / "share" / "mime")],
        ["update-desktop-database", str(tmp_path / "share" / "applications")],
        ["xdg-mime", "default", "vettore.desktop", *LINUX_TYPES],
    ]
    assert not (tmp_path / "config" / "mimeapps.list").exists()   # xdg-mime's job


def test_linux_undo_removes_only_what_vettore_added(tmp_path):
    assoc, _ = linux(tmp_path)
    assoc.enable()
    mimeapps = tmp_path / "config" / "mimeapps.list"
    mimeapps.write_text(mimeapps.read_text(encoding="utf-8")
                        + "[Added Associations]\nchemical/x-genbank=vettore.desktop;other.desktop;\n"
                        + "text/plain=org.kde.kate.desktop;\n", encoding="utf-8")
    assoc.disable()
    assert not (tmp_path / "share" / "applications" / "vettore.desktop").exists()
    assert not (tmp_path / "share" / "mime" / "packages" / "vettore-sequences.xml").exists()
    text = mimeapps.read_text(encoding="utf-8")
    assert "vettore.desktop" not in text
    assert "chemical/x-genbank=other.desktop;" in text
    assert "text/plain=org.kde.kate.desktop;" in text
    assert "[Added Associations]" in text
    assert assoc.status()["enabled"] is False


def test_linux_undo_removes_a_defaults_file_that_held_only_vettore(tmp_path):
    assoc, _ = linux(tmp_path)
    assoc.enable()
    assoc.disable()
    assert not (tmp_path / "config" / "mimeapps.list").exists()


def test_linux_status_after_the_folder_was_moved(tmp_path):
    old, _ = linux(tmp_path)
    old.enable()
    moved = LinuxAssociations("/opt/Vettore/Vettore", data_home=tmp_path / "share",
                              config_home=tmp_path / "config", which=lambda _: None)
    status = moved.status()
    assert (status["enabled"], status["current"]) == (False, False)
    assert status["registered_executable"] == LINUX_EXE


# ---- through the API --------------------------------------------------------------------

def client() -> TestClient:
    return TestClient(create_app("t0ken", None), base_url="http://127.0.0.1")


def test_a_copy_run_from_source_does_not_offer_it():
    c = client()
    status = c.get("/api/desktop/file-associations", headers=H).json()
    assert status["available"] is False
    assert "downloaded" in status["reason"]
    res = c.post("/api/desktop/file-associations", headers=H, json={"enabled": True})
    assert res.status_code == 409 and res.json()["detail"] == status["reason"]


def test_the_api_switches_the_associations_on_and_off(monkeypatch):
    reg = FakeWinreg()
    assoc, _ = windows(reg)
    monkeypatch.setattr(associations, "for_this_system", lambda: assoc)
    c = client()
    dry = c.post("/api/desktop/file-associations", headers=H,
                 json={"enabled": True, "dry_run": True}).json()
    assert dry["enabled"] is False and len(dry["plan"]) == 3 + 2 * len(EXTENSIONS)
    assert reg.keys == {}
    on = c.post("/api/desktop/file-associations", headers=H, json={"enabled": True}).json()
    assert on["enabled"] is True
    assert c.get("/api/desktop/file-associations", headers=H).json() == on
    off = c.post("/api/desktop/file-associations", headers=H, json={"enabled": False}).json()
    assert off["enabled"] is False and not reg.has(rf"{CLASSES}\{PROGID}")


def test_a_refusal_is_explained(monkeypatch):
    assoc, _ = windows(FakeWinreg(fail_on=PROGID))
    monkeypatch.setattr(associations, "for_this_system", lambda: assoc)
    res = client().post("/api/desktop/file-associations", headers=H, json={"enabled": True})
    assert res.status_code == 422
    detail = res.json()["detail"]
    assert detail.startswith("Vettore could not change which program opens sequence files: "
                             "Access is denied. Nothing was changed.")
    assert "log file" in detail
