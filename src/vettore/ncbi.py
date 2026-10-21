"""File → Open from NCBI…: a public GenBank record fetched by its accession number.

This is the only place where Vettore uses the internet, and only when the user asks. The
request carries the accession and nothing else (no sequence, no file, no account): it goes to
NCBI's E-utilities ``efetch`` service, which answers with the GenBank flat file. The record is
then read like an opened file and becomes an unsaved document.

Offline, or when NCBI does not answer, the user gets one sentence saying so; nothing else in
Vettore depends on this module.
"""

from __future__ import annotations

import logging
import re
import ssl
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from . import fileio

router = APIRouter(prefix="/api/ncbi")

log = logging.getLogger(__name__)

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TIMEOUT = 30            # seconds; a large record on a slow line still arrives
MAX_BYTES = 16_000_000  # a GenBank file this size holds several megabases: a genome, not a vector

#: GenBank, EMBL, DDBJ and RefSeq nucleotide accessions, with or without the version:
#: MT612431.1, M77789, NC_001422.1, NZ_CP009072.1, and WGS contigs such as AAAA02000001.
ACCESSION = re.compile(r"^(?:[A-Z]{2}_[A-Z]{0,6}|[A-Z]{1,6})\d{5,12}(?:\.\d{1,3})?$")


class FetchRequest(BaseModel):
    accession: str


class Offline(Exception):
    """NCBI could not be reached: no connection, a proxy that refused, or a timeout."""


def clean_accession(raw: str) -> str:
    """The accession as NCBI writes it, or an error that says what an accession looks like."""
    accession = raw.strip().upper()
    if not accession:
        raise HTTPException(422, "Type an accession number, for example MT612431.1.")
    if not ACCESSION.match(accession):
        raise HTTPException(422, f"“{raw.strip()}” is not a GenBank accession number. It looks "
                                 "like letters followed by digits, for example MT612431.1 or "
                                 "NC_001422.1.")
    return accession


def _ssl_context() -> ssl.SSLContext:
    """The system's certificates (on Windows, those of the certificate store, so a university
    proxy that inspects traffic keeps working), plus certifi's when it is installed."""
    context = ssl.create_default_context()
    try:
        import certifi
    except ImportError:
        return context
    context.load_verify_locations(certifi.where())
    return context


def download(accession: str) -> bytes:
    """The GenBank flat file of ``accession`` (``gbwithparts``: a record assembled from parts
    comes with its sequence). Raises :class:`Offline` when NCBI cannot be reached and
    ``HTTPException`` when it answers that there is no such record."""
    query = urllib.parse.urlencode({"db": "nuccore", "id": accession,
                                    "rettype": "gbwithparts", "retmode": "text"})
    try:
        with urllib.request.urlopen(f"{EFETCH_URL}?{query}", timeout=TIMEOUT,
                                    context=_ssl_context()) as response:
            content = response.read(MAX_BYTES + 1)
    except urllib.error.HTTPError as exc:
        log.warning("NCBI answered %s for %s", exc.code, accession)
        if exc.code in (400, 404):
            raise HTTPException(404, not_found(accession)) from exc
        raise HTTPException(502, f"NCBI could not send {accession} right now (it answered "
                                 f"with error {exc.code}). Try again in a few minutes.") from exc
    except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
        log.warning("NCBI could not be reached for %s: %r", accession, exc)
        raise Offline(str(exc)) from exc
    if len(content) > MAX_BYTES:
        raise HTTPException(422, f"{accession} is too large to open in Vettore (more than "
                                 f"{MAX_BYTES // 1_000_000} MB: a chromosome or a genome). "
                                 "Download the region you need from the NCBI website and open "
                                 "that file.")
    return content


def not_found(accession: str) -> str:
    return (f"NCBI has no nucleotide record {accession}. Check the accession number (for "
            "example MT612431.1); a protein accession cannot be opened as DNA.")


def fetch(accession: str) -> dict[str, Any]:
    """Open result for the UI, like an opened file with no path."""
    content = download(accession)
    # efetch answers an unknown id with an error text, sometimes with status 200.
    if not content.lstrip().startswith(b"LOCUS"):
        log.warning("NCBI sent no GenBank record for %s: %r", accession, content[:200])
        raise HTTPException(404, not_found(accession))
    result = fileio.read_file(f"{accession}.gb", content).as_dict()
    return {"path": None, "accession": accession, **result}


@router.post("/fetch")
def fetch_record(body: FetchRequest) -> dict[str, Any]:
    accession = clean_accession(body.accession)
    try:
        return fetch(accession)
    except Offline as exc:
        raise HTTPException(503, "Vettore could not reach NCBI: this computer seems to be "
                                 "offline, or the network does not let Vettore through. Check "
                                 "the internet connection and try again; everything else in "
                                 "Vettore works without it.") from exc
