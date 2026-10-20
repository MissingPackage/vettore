# HANDOFF — Vettore

**Objective:** an offline plasmid editor for daily bench work, with a familiar workflow and
look; biology computed by tested engines (pydna/Biopython).


**Repo:** github.com/MissingPackage/vettore (private). Spec: `specs/vettore.md`.
Research: `docs/research/`.

## Current state (2026-10-20 05:30 UTC)

| Slice | State | Evidence |
|---|---|---|
| S0 spike | done | CI green; Windows screenshot shows the native window with a .dna map |
| S1 daily driver | done, not yet user-tested | 36 pytest, 13 vitest, 7 e2e green; Windows window verified after the js_api freeze fix (2c102e3) |
| S2 cloning | in progress | engine contract `src/vettore/cloning/base.py`; dialog framework `frontend/src/cloning/`; engine modules still placeholders |
| S3 verification tools | not started | gel, Sanger alignment, common features, primer3 design |
| S4 polish | not started | |

## Next steps (in order)

1. S2 engine: restriction.py → overlap.py (PCR, Gibson/HiFi/In-Fusion, mutagenesis) →
   golden_gate.py → gateway_ta.py, each with hand-computed expected products in tests.
2. Wire the Actions menu to `frontend/src/cloning/registry.ts`; products open as new tabs with
   History entries; e2e for restriction cloning and Gibson.
3. Release v0.2.0 (tag → GitHub Release with Windows/Linux zips) + short user guide.
4. S3: Simulate Agarose Gel, Align to Reference (.ab1), Detect Common Features (curated DB),
   primer design (primer3).

## Decisions taken (see spec "Architecture")

- UI in English. GenBank is the save format; .dna is read (Biopython + sgffp).
- Enzymes: nonredundant commercial REBASE set, cross-checked OVE vs Biopython.
- Tm presets: Standard (default), NEB Q5, NEB Taq.

## Open questions (spec §Open questions)

Cloning methods she uses; whether writing .dna matters; Sanger priority; network policy for
NCBI/Addgene import; licence (GPL-3.0-or-later implied by bundling primer3).

## Gotchas

- pywebview exposes public js_api attributes recursively: keep native objects private (`_window`).
- Blueprint 3 typings lack `children` with React 18: cast container components.
- Vite 8's CSS minifier rejects Blueprint 3 CSS: `cssMinify: false`.
- Background subagents hit the session limit on 2026-10-19 night; prefer sequential work.
