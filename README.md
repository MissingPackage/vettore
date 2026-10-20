# Vettore

Offline desktop plasmid editor and cloning simulator with a familiar workflow.
Runs on Windows 10/11 and Linux, without cloud services, accounts or Docker.

### ⬇ Download

**[Vettore for Windows](https://github.com/MissingPackage/vettore/releases/latest/download/Vettore-windows-x64.zip)**
 · [Vettore for Linux](https://github.com/MissingPackage/vettore/releases/latest/download/Vettore-linux-x64.zip)

Unzip anywhere, then run `Vettore.exe` (Windows) or `Vettore` (Linux). Nothing is installed and
no account is needed. On Windows, SmartScreen warns about an unsigned app the first time:
choose **More info → Run anyway**. Start here: **[user guide](docs/USER-GUIDE.md)**.

Vettore is built on open,
tested libraries — [Open Vector Editor](https://github.com/TeselaGen/tg-oss) for the editor UI and
[pydna](https://github.com/pydna-group/pydna), [Biopython](https://biopython.org) and
[primer3-py](https://github.com/libnano/primer3-py) for the biology.

**[User guide](docs/USER-GUIDE.md)** · status and roadmap: [`specs/vettore.md`](specs/vettore.md).

## Install (users)

Download the zip for your system from the latest [Release](../../releases), unzip it anywhere
and run `Vettore.exe` (Windows) or `Vettore` (Linux). No installation or admin rights are
needed. On Windows, SmartScreen may warn about an unsigned app: choose *More info → Run anyway*.

What it does: opens `.dna`, GenBank, FASTA and `.ab1` files; map, sequence, enzymes,
features, primers and history views; restriction cloning, PCR, primer design, mutagenesis and
Gibson/HiFi/In-Fusion assembly; saves GenBank and exports maps as SVG/PNG.

Logs are written to `~/.vettore/vettore.log` (`%USERPROFILE%\.vettore\vettore.log` on Windows).

## Develop

Requirements: Python 3.12 with [uv](https://docs.astral.sh/uv/), Node.js 24.

```sh
uv sync --extra dev
(cd frontend && npm ci && npm run build)
uv run vettore path/to/file.dna     # native window
uv run vettore --browser            # system browser instead of a native window
```

Tests:

```sh
uv run pytest                                   # Python engine and API
(cd frontend && npx playwright test)            # end-to-end, browser mode
```

Packaging (one-folder build in `dist/Vettore`):

```sh
rm -rf src/vettore/static && cp -r frontend/dist src/vettore/static
uv run pyinstaller packaging/vettore.spec --noconfirm
uv run python packaging/smoke_test.py dist/Vettore/Vettore
```

On Linux the native window uses Qt WebEngine and needs the usual X11/xcb libraries
(on Debian/Ubuntu: `libxcb-cursor0 libxcb-icccm4 libxkbcommon-x11-0 libegl1 libnss3`).

## Architecture

| Layer | Technology |
|---|---|
| UI | React + Open Vector Editor (`frontend/`) |
| Bridge | FastAPI on `127.0.0.1`, random port, per-session token (`src/vettore/server.py`) |
| Biology engine | pydna, Biopython, primer3-py (`src/vettore/`) |
| Window | pywebview: Edge WebView2 on Windows, Qt WebEngine on Linux (`src/vettore/app.py`) |
| Packaging | PyInstaller (`packaging/`), built and smoke-tested in GitHub Actions |

## Licence

GPL-3.0-or-later — see [`LICENSE`](LICENSE). The packaged build links primer3-py
(GPL-2.0-or-later) in-process, so the whole distributed work is under the GNU GPL.
Third-party components and data sources are listed in
[`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md).

