# Third-party components

Vettore is distributed as a bundle that contains the components below. Their licences apply to
their own code; Vettore's own source is under GPL-3.0-or-later (see `LICENSE`).

## Why GPL-3.0-or-later

The packaged application loads **primer3-py** (GPL-2.0-or-later) in the same process, so the
distributed work as a whole is released under the GNU GPL. Version 3 is used because the bundle
also contains Apache-2.0 components, which are compatible with GPLv3 but not with GPLv2.

## Python

| Component | Licence | Used for |
|---|---|---|
| [pydna](https://github.com/pydna-group/pydna) | BSD-3-Clause | Cloning simulation (digestion, ligation, assembly, PCR) |
| [Biopython](https://biopython.org) | Biopython License (BSD-style) | Sequence parsing, restriction enzymes, Tm |
| [primer3-py](https://github.com/libnano/primer3-py) | **GPL-2.0-or-later** | Primer thermodynamics |
| [sgffp](https://github.com/merv1n34k/sgffp) | MIT | `.dna` details (colours, primer tails) |
| [pywebview](https://pywebview.flowrl.com) | BSD-3-Clause | Native window |
| [FastAPI](https://fastapi.tiangolo.com) / [Starlette](https://www.starlette.io) | MIT / BSD-3-Clause | Local API |
| [uvicorn](https://www.uvicorn.org) | BSD-3-Clause | Local server |
| [pydantic](https://docs.pydantic.dev) | MIT | Request validation |
| [PyQt6 / Qt WebEngine](https://www.qt.io) (Linux build only) | GPL-3.0 / LGPL-3.0 | Browser engine for the native window |

## JavaScript

| Component | Licence | Used for |
|---|---|---|
| [Open Vector Editor and @teselagen/*](https://github.com/TeselaGen/tg-oss) | MIT (`licenses/teselagen-tg-oss-MIT.txt`) | Map, sequence and annotation editor |
| [Blueprint](https://blueprintjs.com) | Apache-2.0 | UI components |
| [React](https://react.dev), [Redux](https://redux.js.org) and related packages | MIT | UI framework |

The published npm tarballs of `@teselagen/*` carry no licence file; the text included here is
the one published in their source repository, `TeselaGen/tg-oss`.

## Data

- **Restriction enzymes** — derived from REBASE through Biopython's `Bio.Restriction`.
  REBASE: Roberts RJ, Vincze T, Posfai J, Macelis D. *REBASE — a database for DNA restriction
  and modification: enzymes, genes and genomes.* Nucleic Acids Res. <http://rebase.neb.com>.
  REBASE data may be redistributed free of charge with attribution; Vettore is free of charge.
  The per-enzyme methylation sensitivities in `frontend/src/data/methylation.json` come from
  REBASE's `msbuffmore.txt` (see `scripts/build_methylation_table.py`).
- **DNA ladders** — band sizes taken from the vendors' own public product documentation (New
  England Biolabs, Thermo Fisher, Invitrogen). Product names and trademarks belong to their
  owners and are used only to identify the ladders.
- **Common features** — `src/vettore/data/common_features.json` holds subsequences of public
  GenBank/RefSeq records fetched from NCBI, each entry citing its accession and coordinate
  range (see `scripts/build_feature_database.py`). NCBI places no restriction on reuse of
  these records.
- **Test fixtures** in `tests/data/` — copied from the Biopython test suite (Biopython License /
  BSD-3-Clause); see `tests/data/README.md`.
