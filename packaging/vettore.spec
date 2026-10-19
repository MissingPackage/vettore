# PyInstaller spec for Vettore (one-folder build). Run from the repo root:
#   pyinstaller packaging/vettore.spec --noconfirm
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules

root = Path(SPECPATH).parent
static = root / "src" / "vettore" / "static"
if not (static / "index.html").is_file():
    raise SystemExit("Build the frontend first and copy frontend/dist to src/vettore/static")

datas = [(str(static), "vettore/static")]
binaries = []
hiddenimports = collect_submodules("uvicorn") + collect_submodules("Bio")
for package in ("pydna", "primer3"):
    d, b, h = collect_all(package)
    datas += d
    binaries += b
    hiddenimports += h

a = Analysis(
    [str(root / "packaging" / "launcher.py")],
    pathex=[str(root / "src")],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=["tkinter", "matplotlib", "IPython", "notebook"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Vettore",
    console=False,
    icon=None,
)
coll = COLLECT(exe, a.binaries, a.datas, name="Vettore")
