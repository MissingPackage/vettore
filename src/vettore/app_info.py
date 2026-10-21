"""Application-level routes under /api/app: version, licence, sample constructs.

The sample constructs are public GenBank records shipped with the application and read from
the installed package, so the welcome screen works offline and in a packaged build. Their
provenance (accession, version, retrieval date, licence note) is in
``data/samples/README.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

from fastapi import APIRouter, HTTPException

from . import fileio, settings

router = APIRouter(prefix="/api/app")

LICENCE = "GPL-3.0-or-later"

#: Components a user may want to know about, with the licence that covers them. The full list
#: with links is in THIRD-PARTY-NOTICES.md; this is the short version shown in the About box.
THIRD_PARTY = [
    {"name": "pydna", "licence": "BSD-3-Clause", "used_for": "cloning simulation"},
    {"name": "Biopython", "licence": "Biopython License (BSD-style)",
     "used_for": "sequence parsing, restriction enzymes, melting temperature"},
    {"name": "primer3-py", "licence": "GPL-2.0-or-later", "used_for": "primer thermodynamics"},
    {"name": "Open Vector Editor", "licence": "MIT", "used_for": "map and sequence editor"},
    {"name": "Blueprint", "licence": "Apache-2.0", "used_for": "user-interface components"},
    {"name": "React and Redux", "licence": "MIT", "used_for": "user-interface framework"},
    {"name": "FastAPI and uvicorn", "licence": "MIT / BSD-3-Clause", "used_for": "local engine"},
    {"name": "pywebview", "licence": "BSD-3-Clause", "used_for": "native window"},
    {"name": "REBASE (through Biopython)", "licence": "free with attribution",
     "used_for": "restriction enzyme data"},
    {"name": "NCBI GenBank records", "licence": "public domain",
     "used_for": "common-feature catalogue and sample constructs"},
]


@dataclass(frozen=True)
class Sample:
    """A construct offered on the welcome screen."""

    id: str
    label: str            # the name the document gets when opened
    description: str      # one line, in the words a bench scientist would use
    filename: str
    accession: str


#: Provenance for every entry is in data/samples/README.md.
SAMPLES: tuple[Sample, ...] = (
    Sample("puc19", "pUC19",
           "The classic cloning vector, 2686 bp. The GenBank record is lightly annotated: "
           "a good first try for Features → Detect Common Features",
           "M77789.2.gb", "M77789.2"),
    Sample("paav-ctnt-sacas9", "pAAV-cTNT-SaCas9-U6",
           "AAV transfer plasmid, 7079 bp: SaCas9 and a guide RNA between two annotated ITRs",
           "PZ629150.1.gb", "PZ629150.1"),
    Sample("paav-ef1a-dio", "pAAV-EF1A-DIO-Gephyrin.FingR-GFP",
           "Cre-dependent AAV transfer plasmid, 7204 bp: a GFP reporter in a DIO cassette",
           "MT612431.1.gb", "MT612431.1"),
)


def app_version() -> str:
    try:
        return version("vettore")
    except PackageNotFoundError:
        return "dev"


def log_path() -> Path:
    """Where the app writes its log (see vettore.app.main): ~/.vettore/vettore.log, or under
    $VETTORE_HOME when that is set."""
    return settings.log_path()


def sample_bytes(sample: Sample) -> bytes:
    # vettore.data is the package; samples/ is a plain directory inside it.
    return resources.files("vettore.data").joinpath("samples", sample.filename).read_bytes()


def sample_documents(sample: Sample) -> list[dict[str, Any]]:
    """Read a sample into editor documents, named after the construct rather than after the
    record's LOCUS line (an old record's locus name is an internal identifier)."""
    documents = fileio.read_documents(sample.filename, sample_bytes(sample))
    for document in documents:
        document["name"] = sample.label
    return documents


@router.get("/about")
def about() -> dict[str, Any]:
    return {
        "version": app_version(),
        "licence": LICENCE,
        "log_path": str(log_path()),
        "third_party": THIRD_PARTY,
    }


# The welcome screen draws a small map of each sample from its longest features.
MINI_MAP_FEATURES = 12


def sample_outline(sample: Sample) -> dict[str, Any]:
    """Length, topology and the longest features of a sample, for its card on the welcome screen."""
    document = sample_documents(sample)[0]
    length = len(document["sequence"])
    # A feature across the origin has end < start: the modulo gives its real length.
    span = {id(f): (f["end"] - f["start"]) % length + 1 for f in document.get("features", [])}
    # Records that annotate their whole backbone (pUC19's "pBR322") would draw a closed ring.
    features = [f for f in document.get("features", [])
                if f.get("type") != "source" and span[id(f)] <= length * 0.6]
    features.sort(key=lambda f: span[id(f)], reverse=True)
    return {
        "length": length,
        "circular": bool(document.get("circular")),
        "features": [{"start": f["start"], "end": f["end"], "type": f.get("type", "misc_feature")}
                     for f in features[:MINI_MAP_FEATURES]],
    }


@router.get("/samples")
def samples() -> dict[str, Any]:
    return {"samples": [{"id": s.id, "label": s.label, "description": s.description,
                         "accession": s.accession, **sample_outline(s)} for s in SAMPLES]}


@router.get("/samples/{sample_id}")
def open_sample(sample_id: str) -> dict[str, Any]:
    sample = next((s for s in SAMPLES if s.id == sample_id), None)
    if sample is None:
        raise HTTPException(404, f"Unknown sample construct: {sample_id}")
    # path stays None: a sample opens as an unsaved document, so "Save" asks where to put it
    # instead of writing inside the installed application.
    return {"path": None, "documents": sample_documents(sample)}
