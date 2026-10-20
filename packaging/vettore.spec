# PyInstaller spec for Vettore (one-folder build). Run from the repo root:
#   pyinstaller packaging/vettore.spec --noconfirm
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_submodules, copy_metadata

root = Path(SPECPATH).parent
static = root / "src" / "vettore" / "static"
if not (static / "index.html").is_file():
    raise SystemExit("Build the frontend first and copy frontend/dist to src/vettore/static")

datas = [(str(static), "vettore/static"),
         # Engine data read through importlib.resources: common features, sample constructs.
         (str(root / "src" / "vettore" / "data"), "vettore/data")]
binaries = []
hiddenimports = collect_submodules("uvicorn") + collect_submodules("Bio")

# pydna reads these distributions' metadata at import time (importlib.metadata.version).
for distribution in ("seguid", "pydna", "biopython", "primer3-py", "sgffp", "pydivsufsort",
                     "networkx", "fastapi", "uvicorn", "pywebview"):
    try:
        datas += copy_metadata(distribution)
    except Exception as exc:  # a distribution that is not installed is not needed either
        print(f"vettore.spec: no metadata for {distribution}: {exc}")
for package in ("pydna", "primer3", "sgffp"):
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
    icon=str(root / "frontend" / "public" / "vettore.ico"),
)
coll = COLLECT(exe, a.binaries, a.datas, name="Vettore")
