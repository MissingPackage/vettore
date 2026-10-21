"""File → Open from NCBI…: a public GenBank record fetched by its accession number.

This is the only place where Vettore uses the internet, and only when the user asks. The
request carries the accession and nothing else (no sequence, no file, no account): it goes to
NCBI's E-utilities ``efetch`` service, which answers with the GenBank flat file. The record is
then read like an opened file and becomes an unsaved document.

Offline, or when NCBI does not answer, the user gets one sentence saying so; nothing else in
Vettore depends on this module.
"""

from __future__ import annotations

import http.client
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
TIMEOUT = 30              # seconds; a large record on a slow line still arrives
MAX_BP = 2_000_000        # a BAC or a viral genome opens; a bacterial genome would freeze the map
MAX_BYTES = 16_000_000    # the download stops here whatever the LOCUS line says

#: GenBank, EMBL, DDBJ and RefSeq nucleotide accessions, with or without the version:
#: MT612431.1, M77789, NC_001422.1, NZ_CP009072.1, and WGS contigs such as AAAA02000001.
ACCESSION = re.compile(r"^(?:[A-Z]{2}_[A-Z]{0,6}|[A-Z]{1,6})\d{5,12}(?:\.\d{1,3})?$")
#: The accession inside pasted text such as "GenBank: MT612431.1" or "MT612431.1 (circular)".
ACCESSION_IN_TEXT = re.compile(
    r"\b(?:[A-Z]{2}_[A-Z]{0,6}|[A-Z]{1,6})\d{5,12}(?:\.\d{1,3})?(?![.\d])")
#: A region written after the accession, as NCBI and BLAST print it: "NC_045512.2:21563-25384",
#: "NC_000011.10:c27722058-27654893".
REGION = re.compile(r":\s*C?\d[\d,]*\s*(?:-|\.\.)\s*\d")
#: Identifiers of other NCBI and EBI databases that look like accessions: SRA runs,
#: experiments, studies and samples, BioProjects, BioSamples, PubMed Central.
OTHER_DATABASE = re.compile(r"^(?:[SED]R[RXPSA]|PRJ[A-Z]{1,2}|SAM[NED]|PMC)\d")
#: Protein accessions: RefSeq NP_/XP_/WP_/YP_/AP_, and INSDC three letters and five or seven
#: digits (nucleotide accessions have one, two, four, five or six letters).
PROTEIN = re.compile(r"^(?:[NXWYA]P_\d+|[A-Z]{3}(?:\d{5}|\d{7}))(?:\.\d+)?$")
#: What efetch writes, with status 200 or 400, when the id names no nucleotide record.
#: Compared without spaces: efetch sometimes spaces the letters out ("F a i l e d  t o …").
NOT_FOUND_TEXTS = (b"failedtounderstandid", b"cannotretrieve", b"idlistisempty",
                   b"invaliduid", b"cannotgetdocumentsummary")


class FetchRequest(BaseModel):
    accession: str


class Offline(Exception):
    """NCBI could not be reached: no connection, a proxy that refused, or a timeout."""


class Unverified(Offline):
    """A connection was made but its certificate could not be checked (a network that
    inspects traffic with a certificate this computer does not trust)."""


def clean_accession(raw: str) -> str:
    """The accession as NCBI writes it, or an error that says what an accession looks like.
    An accession copied with words around it ("GenBank: MT612431.1") is picked out."""
    typed = raw.strip().upper()
    if not typed:
        raise HTTPException(422, "Type an accession number, for example MT612431.1.")
    if REGION.search(typed):
        # Opening the whole record instead would hand over a different sequence than the one
        # asked for, without a word.
        raise HTTPException(422, f"“{raw.strip()}” asks for a region of a record; Vettore "
                                 "opens whole records only. Type the accession alone and "
                                 "select the region in Vettore, or, for a region of a "
                                 "chromosome, download that region from the NCBI website as a "
                                 "GenBank file and open it.")
    if not ACCESSION.match(typed):
        found = sorted(set(ACCESSION_IN_TEXT.findall(typed)))
        if len(found) != 1:
            raise HTTPException(422, f"“{raw.strip()}” is not a GenBank accession number. It "
                                     "looks like letters followed by digits, for example "
                                     "MT612431.1 or NC_001422.1.")
        typed = found[0]
    if OTHER_DATABASE.match(typed):
        raise HTTPException(422, f"{typed} is not a GenBank record: it looks like an identifier "
                                 "of another database (sequencing reads, a project, a sample or "
                                 "an article). Use the accession of a nucleotide record, for "
                                 "example MT612431.1.")
    if PROTEIN.match(typed):
        raise HTTPException(422, f"{typed} is a protein accession: Vettore opens DNA. Use the "
                                 "accession of the nucleotide record (on the protein's NCBI "
                                 "page, under “Coded by” or “DBSOURCE”).")
    return typed


def _ssl_context() -> ssl.SSLContext:
    """The system's certificates (on Windows, those of the certificate store, so a university
    proxy that inspects traffic keeps working), plus certifi's."""
    context = ssl.create_default_context()
    try:
        import certifi
    except ImportError:   # a dependency; kept optional so a stripped install still runs
        return context
    context.load_verify_locations(certifi.where())
    return context


def _too_large(accession: str, what: str) -> HTTPException:
    return HTTPException(422, f"{accession} is too large to open in Vettore ({what}: a "
                              "chromosome or a genome). Download the region you need from the "
                              "NCBI website and open that file.")


def _length_on_locus_line(line: bytes) -> int | None:
    match = re.search(rb"\s(\d+)\s+bp\b", line)
    return int(match.group(1)) if match else None


def download(accession: str) -> bytes:
    """The GenBank flat file of ``accession`` (``gbwithparts``: a record assembled from parts
    comes with its sequence). Raises :class:`Offline` when NCBI cannot be reached and
    ``HTTPException`` when it answers that there is no such record."""
    # "tool" names the program, as NCBI's usage policy asks; nothing about the user is sent.
    query = urllib.parse.urlencode({"db": "nuccore", "id": accession, "rettype": "gbwithparts",
                                    "retmode": "text", "tool": "vettore"})
    try:
        with urllib.request.urlopen(f"{EFETCH_URL}?{query}", timeout=TIMEOUT,
                                    context=_ssl_context()) as response:
            # The LOCUS line gives the length: a genome is refused before it is downloaded.
            first = response.readline(4096)
            length = _length_on_locus_line(first) if first.startswith(b"LOCUS") else None
            if length is not None and length > MAX_BP:
                raise _too_large(accession, f"{length:,} bp")
            content = first + response.read(max(0, MAX_BYTES + 1 - len(first)))
    except urllib.error.HTTPError as exc:
        log.warning("NCBI answered %s for %s", exc.code, accession)
        if exc.code in (400, 404):
            raise HTTPException(404, not_found(accession)) from exc
        if exc.code in (401, 403, 407):   # a proxy or filter of the local network said no
            raise Offline(f"HTTP {exc.code}") from exc
        raise _busy(accession, f"it answered with error {exc.code}") from exc
    except urllib.error.URLError as exc:
        log.warning("NCBI could not be reached for %s: %r", accession, exc)
        if isinstance(exc.reason, ssl.SSLCertVerificationError):
            raise Unverified(str(exc)) from exc
        raise Offline(str(exc)) from exc
    except ssl.SSLCertVerificationError as exc:
        log.warning("NCBI's certificate could not be verified for %s: %r", accession, exc)
        raise Unverified(str(exc)) from exc
    except (http.client.HTTPException, OSError) as exc:   # timeouts, resets, a cut transfer
        log.warning("The download of %s was interrupted: %r", accession, exc)
        raise Offline(str(exc)) from exc
    if len(content) > MAX_BYTES:
        raise _too_large(accession, f"more than {MAX_BYTES // 1_000_000} MB")
    return content


def _busy(accession: str, why: str) -> HTTPException:
    return HTTPException(502, f"NCBI could not send {accession} right now ({why}). Try again "
                              "in a few minutes.")


def not_found(accession: str) -> str:
    return (f"NCBI has no nucleotide record {accession}. Check the accession number (for "
            "example MT612431.1); a protein accession cannot be opened as DNA.")


def _version(content: bytes, accession: str) -> str:
    """The record's accession.version (``VERSION`` line), so a tab opened as "M77789" is named
    after the version it holds, "M77789.2"."""
    match = re.search(rb"^VERSION\s+(\S+)", content, re.MULTILINE)
    version = match.group(1).decode("ascii", "replace").upper() if match else ""
    return version if ACCESSION.match(version) else accession


def fetch(accession: str) -> dict[str, Any]:
    """Open result for the UI, like an opened file with no path."""
    content = download(accession)
    if not content.lstrip().startswith(b"LOCUS"):
        log.warning("NCBI sent no GenBank record for %s: %r", accession, content[:200])
        # efetch answers an unknown id with an error text, sometimes with status 200; anything
        # else (an HTML page during an outage, an empty answer) is NCBI not working right now.
        squeezed = re.sub(rb"\s+", b"", content.lower())
        if any(text in squeezed for text in NOT_FOUND_TEXTS):
            raise HTTPException(404, not_found(accession))
        raise _busy(accession, "it did not send a GenBank record")
    version = _version(content, accession)
    if not re.search(rb"^ORIGIN", content, re.MULTILINE) and re.search(
            rb"^(?:WGS|TSA|TLS)\s", content, re.MULTILINE):
        raise HTTPException(422, f"{version} is the master record of a sequencing project: it "
                                 "lists contigs but holds no sequence. Open one of its contig "
                                 "accessions instead (on its NCBI page, under “WGS”).")
    if not re.search(rb"^ORIGIN", content, re.MULTILINE) and re.search(
            rb"^CONTIG\s", content, re.MULTILINE):
        raise HTTPException(422, f"NCBI sent {version} as a list of its parts, without the "
                                 "sequence. Download it from the NCBI website as a GenBank file "
                                 "(“Show sequence”) and open that file.")
    try:
        result = fileio.read_file(f"{version}.gb", content).as_dict()
    except fileio.UnsupportedFile as exc:
        log.warning("Could not read the record NCBI sent for %s: %s", version, exc)
        raise HTTPException(422, f"Vettore could not read the record NCBI sent for {version}. "
                                 "Try again; if it happens again, download it from the NCBI "
                                 "website as a GenBank file and open that file.") from exc
    # The LOCUS line was checked before the download; the sequence itself is the last word.
    longest = max((len(d.get("sequence", "")) for d in result["documents"]), default=0)
    if longest > MAX_BP:
        raise _too_large(version, f"{longest:,} bp")
    return {"path": None, "accession": version, **result}


@router.post("/fetch")
def fetch_record(body: FetchRequest) -> dict[str, Any]:
    accession = clean_accession(body.accession)
    try:
        return fetch(accession)
    except Unverified as exc:
        raise HTTPException(503, "Vettore reached the network but could not verify NCBI's "
                                 "certificate: the network may inspect secure connections. "
                                 "Ask your IT service whether eutils.ncbi.nlm.nih.gov is "
                                 "allowed; everything else in Vettore works without it.") from exc
    except Offline as exc:
        raise HTTPException(503, "Vettore could not reach NCBI: this computer seems to be "
                                 "offline, or the network does not let Vettore through. Check "
                                 "the internet connection and try again; everything else in "
                                 "Vettore works without it.") from exc
