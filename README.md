# Vettore

A desktop plasmid editor and cloning simulator for the bench: open a construct, read its map,
plan the cloning, order the primers, check the clone against its Sanger reads. It runs
entirely on your own computer — Windows 10/11 or Linux — with no account, no cloud service and
nothing to install.

### ⬇ Download

**[Vettore for Windows](https://github.com/MissingPackage/vettore/releases/latest/download/Vettore-windows-x64.zip)**
 · [Vettore for Linux](https://github.com/MissingPackage/vettore/releases/latest/download/Vettore-linux-x64.zip)

Unzip anywhere, then run `Vettore.exe` (Windows) or `Vettore` (Linux). On Windows, SmartScreen
warns about an unsigned app the first time: choose **More info → Run anyway**.
Start here: **[user guide](docs/USER-GUIDE.md)**.

## What it does

**Sequences and maps**
- Opens `.dna`, GenBank (`.gb`, `.gbk`, `.ape`), FASTA and `.ab1` files, keeping features,
  colours, primers (with their 5′ tails) and topology. Saves GenBank; exports FASTA and maps as
  SVG or PNG for figures.
- Map, Sequence, Enzymes, Features, Primers and History views, with editing, undo and find
  (DNA or protein, with mismatches).

**Enzymes**
- REBASE commercial enzymes, one per specificity; sets for unique, unique-and-dual, 6+, type IIS
  and Golden Gate cutters.
- Cut positions and ends per site, and a warning when Dam or Dcm methylation falls inside a
  site that REBASE records as blocked.

**Cloning simulation** — each with a live preview, warnings and the primers to order
- Restriction cloning (insert or delete a fragment, either orientation, directionality and
  self-ligation checks).
- Gibson, NEBuilder HiFi and In-Fusion assembly of one or more inserts, with the homology arms
  designed for you.
- Golden Gate with junction-overhang checks; annealed-oligo cloning for gRNAs, shRNAs and
  linkers, with the oligo pair written from the cut vector's own overhangs.
- PCR (tailed primers, amplicons across the origin), primer design to a target Tm, and
  site-directed mutagenesis with back-to-back primers.

**Verification**
- Sanger reads aligned to the expected construct on either strand, with chromatograms and a
  per-read list of every difference in plasmid coordinates.
- Simulated agarose gels against NEB, Thermo and Invitrogen ladders.
- Automatic annotation of AAV ITRs, common promoters, WPRE, polyA signals, markers, origins,
  fluorescent proteins, lox/FRT sites, Cas9 and epitope tags — each traceable to a public
  GenBank record.

Primer Tm is always shown with the model and conditions behind it, with presets matching NEB's
calculator for Q5 and Taq.

## How it is built

The biology runs on open, tested libraries — [pydna](https://github.com/pydna-group/pydna) and
[Biopython](https://biopython.org) for cloning and sequence handling,
[primer3-py](https://github.com/libnano/primer3-py) for thermodynamics — and every simulation is
tested against sequences worked out by hand. The editor is built on
[Open Vector Editor](https://github.com/TeselaGen/tg-oss). Status and roadmap:
[`specs/vettore.md`](specs/vettore.md).

| Layer | Technology |
|---|---|
| UI | React + Open Vector Editor (`frontend/`) |
| Bridge | FastAPI on `127.0.0.1`, random port, per-session token (`src/vettore/server.py`) |
| Biology engine | pydna, Biopython, primer3-py (`src/vettore/`) |
| Window | pywebview: Edge WebView2 on Windows, Qt WebEngine on Linux (`src/vettore/app.py`) |
| Packaging | PyInstaller (`packaging/`), built and smoke-tested in GitHub Actions |

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
(cd frontend && npm test)                       # unit tests
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

## Licence

GPL-3.0-or-later — see [`LICENSE`](LICENSE). The packaged build links primer3-py
(GPL-2.0-or-later) in-process, so the whole distributed work is under the GNU GPL.
Third-party components and data sources are listed in
[`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md).
