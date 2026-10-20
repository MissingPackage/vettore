"""Open sequence files with Vettore when they are double-clicked.

Offered only by the packaged application (``sys.frozen``): a copy run from source code has no
executable worth pointing the system at. Everything is written for the current user, so no
administrator rights are needed, and undoing it removes exactly what was written.

Windows (under ``HKEY_CURRENT_USER\\Software\\Classes``, which overrides the machine-wide
entries for this user only):

- a ProgID, ``Vettore.SequenceFile``, with the executable's icon and its open command;
- for each extension, the ProgID as the default (the previous default is kept under the
  ProgID and put back on undo) and in ``OpenWithProgids``, so Vettore is listed in Open with;
- then ``SHChangeNotify(SHCNE_ASSOCCHANGED)`` so Explorer picks the change up.

Windows keeps a choice the user made with "Always use this app" (the ``UserChoice`` key,
which programs cannot write). Such extensions are reported as ``overridden``, so the UI can
say how to change them by hand.

Linux (freedesktop.org): a ``vettore.desktop`` entry and a MIME package with the extensions
under ``~/.local/share``, the databases updated when their tools are installed, and Vettore
made the default with ``xdg-mime`` or, without it, in ``~/.config/mimeapps.list``. The MIME
globs carry a higher weight than the shared database's ``*.gb`` (a game cartridge image) and
``*.ape`` (an audio format), which the user has chosen to give to Vettore.
"""

from __future__ import annotations

import configparser
import importlib
import logging
import os
import shutil
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

#: The extensions handed to Vettore (all of them are formats Vettore opens, see fileio.FORMATS).
EXTENSIONS = (".dna", ".gb", ".gbk", ".gbff", ".ape", ".embl", ".fa", ".fasta", ".xdna", ".gck",
              ".ab1")

PROGID = "Vettore.SequenceFile"
PROGID_LABEL = "DNA sequence file"
CLASSES = r"Software\Classes"
#: Where Windows keeps "Always use this app" (the second key on recent Windows 11 builds).
FILE_EXTS = r"Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts"
USER_CHOICES = (FILE_EXTS + r"\{ext}\UserChoice", FILE_EXTS + r"\{ext}\UserChoiceLatest")
PREVIOUS = "PreviousDefaults"     # subkey of the ProgID: the defaults Vettore replaced

DESKTOP_FILE = "vettore.desktop"
MIME_PACKAGE = "vettore-sequences.xml"
#: MIME types per extension group, with the description file managers show.
LINUX_TYPES: dict[str, tuple[str, tuple[str, ...]]] = {
    "chemical/x-genbank": ("GenBank sequence", (".gb", ".gbk", ".gbff", ".ape")),
    "chemical/x-embl-dl-nucleotide": ("EMBL sequence", (".embl",)),
    "chemical/seq-na-fasta": ("FASTA sequence", (".fa", ".fasta")),
    "application/x-dna-sequence": ("DNA sequence file", (".dna",)),
    "application/x-xdna-sequence": ("DNA sequence file (.xdna)", (".xdna",)),
    "application/x-gck-sequence": ("DNA sequence file (.gck)", (".gck",)),
    "application/x-abi-trace": ("Sanger sequencing trace", (".ab1",)),
}


class Unavailable(Exception):
    """This copy of Vettore cannot be set to open files (not the packaged application)."""


def is_packaged() -> bool:
    return bool(getattr(sys, "frozen", False))


def for_this_system() -> WindowsAssociations | LinuxAssociations | None:
    """The associations of this copy, or None when it cannot offer them."""
    if not is_packaged():
        return None
    if sys.platform == "win32":
        return WindowsAssociations(sys.executable)
    if sys.platform.startswith("linux"):
        return LinuxAssociations(sys.executable, _bundled_icon())
    return None


def _bundled_icon() -> Path | None:
    bundle = getattr(sys, "_MEIPASS", None)
    icon = Path(bundle) / "vettore" / "static" / "vettore.png" if bundle else None
    return icon if icon and icon.is_file() else None


def unavailable_status() -> dict[str, Any]:
    return {
        "available": False,
        "enabled": False,
        "current": False,
        "registered_executable": None,
        "extensions": list(EXTENSIONS),
        "overridden": [],
        "platform": sys.platform,
        "reason": ("Only the Vettore application downloaded from the Releases page can be set "
                   "to open sequence files; this copy runs from source code."),
    }


def _same_file(a: str | None, b: str) -> bool:
    return bool(a) and os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


# ---- Windows ----------------------------------------------------------------------------

def _shell_changed() -> None:
    """Tell Explorer that file associations changed (SHCNE_ASSOCCHANGED, SHCNF_IDLIST)."""
    import ctypes

    ctypes.windll.shell32.SHChangeNotify(0x08000000, 0x0000, None, None)


class WindowsAssociations:
    """Per-user file associations in the registry. ``winreg`` and ``notify`` are replaceable
    for tests."""

    platform = "win32"

    def __init__(self, executable: str, winreg: Any = None,
                 notify: Callable[[], None] | None = None) -> None:
        self.executable = executable
        self.reg = winreg if winreg is not None else importlib.import_module("winreg")
        self.notify = notify or _shell_changed

    @property
    def command(self) -> str:
        return f'"{self.executable}" "%1"'

    def plan(self) -> list[str]:
        """What :meth:`enable` writes, one line each, under HKEY_CURRENT_USER."""
        progid = rf"{CLASSES}\{PROGID}"
        lines = [rf"{progid} = {PROGID_LABEL}",
                 rf"{progid}\DefaultIcon = " + f'"{self.executable}",0',
                 rf"{progid}\shell\open\command = {self.command}"]
        for ext in EXTENSIONS:
            lines.append(rf"{CLASSES}\{ext} = {PROGID}")
            lines.append(rf"{CLASSES}\{ext}\OpenWithProgids\{PROGID}")
        return lines

    # registry helpers
    def _key(self, path: str):
        return self.reg.CreateKeyEx(self.reg.HKEY_CURRENT_USER, path, 0, self.reg.KEY_WRITE)

    def _read(self, path: str, name: str = "") -> str | None:
        try:
            with self.reg.OpenKey(self.reg.HKEY_CURRENT_USER, path, 0, self.reg.KEY_READ) as key:
                value, _kind = self.reg.QueryValueEx(key, name)
        except OSError:
            return None
        return value if isinstance(value, str) else None

    def _set(self, path: str, value: str, name: str = "") -> None:
        with self._key(path) as key:
            self.reg.SetValueEx(key, name, 0, self.reg.REG_SZ, value)

    def _delete_value(self, path: str, name: str = "") -> None:
        try:
            with self.reg.OpenKey(self.reg.HKEY_CURRENT_USER, path, 0,
                                  self.reg.KEY_WRITE) as key:
                self.reg.DeleteValue(key, name)
        except OSError:
            pass

    def _delete_tree(self, path: str) -> None:
        try:
            with self.reg.OpenKey(self.reg.HKEY_CURRENT_USER, path, 0,
                                  self.reg.KEY_READ) as key:
                children = []
                while True:
                    try:
                        children.append(self.reg.EnumKey(key, len(children)))
                    except OSError:
                        break
        except OSError:
            return
        for child in children:
            self._delete_tree(rf"{path}\{child}")
        try:
            self.reg.DeleteKey(self.reg.HKEY_CURRENT_USER, path)
        except OSError:
            log.warning("Could not remove HKEY_CURRENT_USER\\%s", path)

    def _delete_key_if_empty(self, path: str) -> None:
        try:
            with self.reg.OpenKey(self.reg.HKEY_CURRENT_USER, path, 0,
                                  self.reg.KEY_READ) as key:
                for enum in (self.reg.EnumKey, self.reg.EnumValue):
                    try:
                        enum(key, 0)
                        return            # still holds something
                    except OSError:
                        pass
            self.reg.DeleteKey(self.reg.HKEY_CURRENT_USER, path)
        except OSError:
            pass

    # the feature
    def enable(self) -> None:
        progid = rf"{CLASSES}\{PROGID}"
        try:
            self._set(progid, PROGID_LABEL)
            self._set(rf"{progid}\DefaultIcon", f'"{self.executable}",0')
            self._set(rf"{progid}\shell\open\command", self.command)
            for ext in EXTENSIONS:
                previous = self._read(rf"{CLASSES}\{ext}")
                if previous and previous != PROGID:
                    self._set(rf"{progid}\{PREVIOUS}", previous, ext)
                self._set(rf"{CLASSES}\{ext}", PROGID)
                self._set(rf"{CLASSES}\{ext}\OpenWithProgids", "", PROGID)
        except OSError:
            log.exception("Could not write the file associations; undoing what was written")
            self.disable()
            raise
        finally:
            self.notify()

    def disable(self) -> None:
        progid = rf"{CLASSES}\{PROGID}"
        for ext in EXTENSIONS:
            key = rf"{CLASSES}\{ext}"
            if self._read(key) == PROGID:
                previous = self._read(rf"{progid}\{PREVIOUS}", ext)
                if previous:
                    self._set(key, previous)
                else:
                    self._delete_value(key)
            self._delete_value(rf"{key}\OpenWithProgids", PROGID)
            self._delete_key_if_empty(rf"{key}\OpenWithProgids")
            self._delete_key_if_empty(key)
        self._delete_tree(progid)
        self.notify()

    def status(self) -> dict[str, Any]:
        command = self._read(rf"{CLASSES}\{PROGID}\shell\open\command")
        registered = _command_executable(command)
        ours = [ext for ext in EXTENSIONS if self._read(rf"{CLASSES}\{ext}") == PROGID]
        overridden = []
        for ext in EXTENSIONS:
            choices = [self._read(path.format(ext=ext), "ProgId") for path in USER_CHOICES]
            if any(c and c != PROGID and "vettore" not in c.lower() for c in choices):
                overridden.append(ext)
        current = _same_file(registered, self.executable)
        return {
            "available": True,
            "enabled": current and len(ours) == len(EXTENSIONS),
            "current": current,
            "registered_executable": registered,
            "extensions": list(EXTENSIONS),
            "overridden": overridden,
            "platform": self.platform,
        }


def _command_executable(command: str | None) -> str | None:
    """The program of an open command such as ``"C:\\Vettore\\Vettore.exe" "%1"``."""
    if not command:
        return None
    command = command.strip()
    if command.startswith('"'):
        end = command.find('"', 1)
        return command[1:end] if end > 0 else None
    return command.split(" ", 1)[0]


# ---- Linux ------------------------------------------------------------------------------

def _desktop_quote(arg: str) -> str:
    """An Exec argument as the Desktop Entry Specification writes it: quoted, with ``"`$\\``
    escaped inside the quotes and ``%`` doubled (field codes), then backslashes doubled again
    (the escape every string value of the file gets)."""
    quoted = '"' + "".join("\\" + c if c in '"`$\\' else c for c in arg).replace("%", "%%") + '"'
    return quoted.replace("\\", "\\\\")


def _desktop_unquote_first(value: str) -> str | None:
    """The first argument of an Exec value written by :func:`_desktop_quote`."""
    value = value.replace("\\\\", "\\").strip()
    if not value.startswith('"'):
        return value.split(" ", 1)[0].replace("%%", "%") or None
    out, i = [], 1
    while i < len(value):
        c = value[i]
        if c == "\\" and i + 1 < len(value):
            out.append(value[i + 1])
            i += 2
            continue
        if c == '"':
            return "".join(out).replace("%%", "%")
        out.append(c)
        i += 1
    return None


class LinuxAssociations:
    """A desktop entry, a MIME package and the default application, per user. The folders,
    ``which`` and ``run`` are replaceable for tests."""

    platform = "linux"

    def __init__(self, executable: str, icon: Path | None = None,
                 data_home: Path | None = None, config_home: Path | None = None,
                 which: Callable[[str], str | None] = shutil.which,
                 run: Callable[..., Any] = subprocess.run) -> None:
        self.executable = executable
        self.icon = icon
        home = Path.home()
        self.data_home = data_home or Path(os.environ.get("XDG_DATA_HOME") or home / ".local/share")
        self.config_home = config_home or Path(os.environ.get("XDG_CONFIG_HOME") or home / ".config")
        self.which = which
        self.run = run

    @property
    def desktop_file(self) -> Path:
        return self.data_home / "applications" / DESKTOP_FILE

    @property
    def mime_package(self) -> Path:
        return self.data_home / "mime" / "packages" / MIME_PACKAGE

    @property
    def mimeapps(self) -> Path:
        return self.config_home / "mimeapps.list"

    def plan(self) -> list[str]:
        return [f"{self.desktop_file}: Exec={_desktop_quote(self.executable)} %F",
                f"{self.mime_package}: " + " ".join(f"*{e}" for e in EXTENSIONS),
                f"default application for {', '.join(LINUX_TYPES)}: {DESKTOP_FILE}"]

    def _desktop_entry(self) -> str:
        lines = ["[Desktop Entry]", "Type=Application", "Name=Vettore",
                 "GenericName=Plasmid editor",
                 "Comment=Open plasmid maps and sequence files, and simulate cloning",
                 f"Exec={_desktop_quote(self.executable)} %F",
                 "Terminal=false", "Categories=Science;Biology;",
                 "MimeType=" + "".join(f"{t};" for t in LINUX_TYPES)]
        if self.icon:
            lines.append(f"Icon={self.icon}")
        return "\n".join(lines) + "\n"

    @staticmethod
    def _mime_xml() -> str:
        out = ['<?xml version="1.0" encoding="UTF-8"?>',
               '<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">']
        for mime, (comment, extensions) in LINUX_TYPES.items():
            out.append(f'  <mime-type type="{mime}">')
            out.append(f"    <comment>{comment}</comment>")
            out += [f'    <glob pattern="*{ext}" weight="80"/>' for ext in extensions]
            out.append("  </mime-type>")
        out.append("</mime-info>")
        return "\n".join(out) + "\n"

    def _tool(self, *cmd: str) -> None:
        if not self.which(cmd[0]):
            return
        try:
            self.run(list(cmd), check=False, timeout=60, stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL)
        except (OSError, subprocess.SubprocessError):
            log.warning("%s failed", cmd[0], exc_info=True)

    def _update_databases(self) -> None:
        self._tool("update-mime-database", str(self.data_home / "mime"))
        self._tool("update-desktop-database", str(self.data_home / "applications"))

    def _read_mimeapps(self) -> configparser.ConfigParser:
        parser = configparser.ConfigParser(interpolation=None, delimiters=("=",))
        parser.optionxform = str   # MIME types and file names are case-sensitive
        try:
            parser.read(self.mimeapps, encoding="utf-8")
        except (OSError, configparser.Error):
            log.warning("Could not read %s", self.mimeapps, exc_info=True)
        return parser

    def _write_mimeapps(self, parser: configparser.ConfigParser) -> None:
        self.mimeapps.parent.mkdir(parents=True, exist_ok=True)
        with open(self.mimeapps, "w", encoding="utf-8") as stream:
            parser.write(stream, space_around_delimiters=False)

    def enable(self) -> None:
        self.desktop_file.parent.mkdir(parents=True, exist_ok=True)
        self.desktop_file.write_text(self._desktop_entry(), encoding="utf-8")
        self.mime_package.parent.mkdir(parents=True, exist_ok=True)
        self.mime_package.write_text(self._mime_xml(), encoding="utf-8")
        self._update_databases()
        if self.which("xdg-mime"):
            self._tool("xdg-mime", "default", DESKTOP_FILE, *LINUX_TYPES)
        else:
            parser = self._read_mimeapps()
            if not parser.has_section("Default Applications"):
                parser.add_section("Default Applications")
            for mime in LINUX_TYPES:
                parser.set("Default Applications", mime, f"{DESKTOP_FILE};")
            self._write_mimeapps(parser)

    def disable(self) -> None:
        for path in (self.desktop_file, self.mime_package):
            path.unlink(missing_ok=True)
        self._update_databases()
        if not self.mimeapps.is_file():
            return
        parser = self._read_mimeapps()
        changed = False
        for section in ("Default Applications", "Added Associations"):
            if not parser.has_section(section):
                continue
            touched = False
            for mime, value in list(parser.items(section)):
                apps = [a for a in value.split(";") if a and a != DESKTOP_FILE]
                if len(apps) != len([a for a in value.split(";") if a]):
                    touched = True
                    if apps:
                        parser.set(section, mime, "".join(f"{a};" for a in apps))
                    else:
                        parser.remove_option(section, mime)
            if touched and not parser.items(section):
                parser.remove_section(section)
            changed = changed or touched
        if not changed:
            return
        if parser.sections():
            self._write_mimeapps(parser)
        else:                      # it held only what Vettore wrote
            self.mimeapps.unlink(missing_ok=True)

    def status(self) -> dict[str, Any]:
        registered = None
        try:
            for line in self.desktop_file.read_text(encoding="utf-8").splitlines():
                if line.startswith("Exec="):
                    registered = _desktop_unquote_first(line[5:])
        except OSError:
            pass
        current = _same_file(registered, self.executable) and self.mime_package.is_file()
        return {
            "available": True,
            "enabled": current,
            "current": current,
            "registered_executable": registered,
            "extensions": list(EXTENSIONS),
            "overridden": [],
            "platform": self.platform,
        }
