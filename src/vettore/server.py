"""Local HTTP API used by the UI. Bound to 127.0.0.1 and protected by a per-session token."""

from __future__ import annotations

import base64
import binascii
import errno
import logging
import os
import secrets
from collections.abc import Callable
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import alignment, app_info, fileio, gel, primer_design, sanger, settings, tm
from . import features as common_features
from .cloning.api import router as cloning_router
from .cloning.base import CloningError

TOKEN_HEADER = "X-Vettore-Token"
ALLOWED_HOSTS = {"127.0.0.1", "localhost"}

log = logging.getLogger(__name__)


class OpenPath(BaseModel):
    path: str


class OpenContent(BaseModel):
    filename: str
    content_b64: str


class OpenText(BaseModel):
    """Text pasted from the clipboard: a sequence, a FASTA record or a GenBank/EMBL record."""

    text: str
    name: str = ""


class SavePath(BaseModel):
    path: str
    sequenceData: dict[str, Any]
    # Used only when the path has no extension Vettore writes: genbank, fasta or dna.
    format: str | None = None


class SaveBinary(BaseModel):
    path: str
    content_b64: str


class DetectRequest(BaseModel):
    sequenceData: dict[str, Any]
    min_identity: float = common_features.DEFAULT_MIN_IDENTITY


class SangerRead(BaseModel):
    filename: str
    content_b64: str


class SangerRequest(BaseModel):
    reference: dict[str, Any]
    reads: list[SangerRead]
    cutoff: float = sanger.DEFAULT_CUTOFF


class GelLane(BaseModel):
    name: str
    sequenceData: dict[str, Any]
    enzymes: list[str] = []


class GelRequest(BaseModel):
    lanes: list[GelLane]
    ladder: str = gel.DEFAULT_LADDER
    percent: float = gel.DEFAULT_PERCENT


class Unsaved(BaseModel):
    count: int


class TmRequest(BaseModel):
    sequence: str


class Export(BaseModel):
    sequenceData: dict[str, Any]
    format: str = "genbank"


def _format(name: str | None) -> str | None:
    """The UI says "dna" for the .dna format; the engine uses Biopython's identifier."""
    return fileio.DNA_FILE_FORMAT if name == "dna" else name


def _os_reason(exc: OSError) -> str:
    """Why a file operation failed, in words; the exception itself goes to the log."""
    if isinstance(exc, PermissionError):
        return ("Vettore is not allowed to use it there, or another program has the file open. "
                "Close it in the other program, or choose another folder such as Documents.")
    if isinstance(exc, IsADirectoryError):
        return "that name belongs to a folder, not a file."
    if isinstance(exc, FileNotFoundError):
        return "the folder does not exist any more."
    if exc.errno == errno.ENOSPC:
        return "the disk is full. Free some space or choose another disk."
    if exc.errno == errno.EROFS:
        return "the disk or folder is read-only. Choose another folder."
    return "the system reported a problem. The details were written to Vettore's log file."


def _write_atomically(path: Path, data: bytes) -> None:
    """Write next to the target, then replace it, so a failed save never truncates a file."""
    tmp = path.with_name(f".{path.name}.vettore-tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    except OSError:
        tmp.unlink(missing_ok=True)
        raise


def _save_failed(path: Path, exc: OSError) -> HTTPException:
    log.error("Could not save %s: %r", path, exc)
    return HTTPException(422, f"“{path.name}” could not be saved in {path.parent}: "
                              f"{_os_reason(exc)}")


def create_app(token: str, static_dir: Path | None,
               startup_files: list[Path] | None = None,
               on_unsaved: Callable[[int], None] | None = None) -> FastAPI:
    app = FastAPI(title="Vettore", docs_url=None, redoc_url=None, openapi_url=None)

    @app.middleware("http")
    async def guard(request: Request, call_next):
        host = (request.headers.get("host") or "").rsplit(":", 1)[0]
        if host not in ALLOWED_HOSTS:
            return JSONResponse({"detail": "Forbidden host"}, status_code=403)
        if request.url.path.startswith("/api/"):
            supplied = request.headers.get(TOKEN_HEADER, "")
            if not secrets.compare_digest(supplied, token):
                return JSONResponse({"detail": "Invalid token"}, status_code=403)
        return await call_next(request)

    @app.exception_handler(CloningError)
    async def cloning_error(_request: Request, exc: CloningError):
        return JSONResponse({"detail": str(exc)}, status_code=422)

    app.include_router(cloning_router)
    # Feature modules own their routes; each exposes ``router``.
    for module in (alignment, primer_design, app_info):
        app.include_router(module.router)

    @app.exception_handler(fileio.UnsupportedFile)
    async def unsupported(_request: Request, exc: fileio.UnsupportedFile):
        return JSONResponse({"detail": str(exc)}, status_code=422)

    @app.exception_handler(fileio.SaveError)
    async def unsaveable(_request: Request, exc: fileio.SaveError):
        return JSONResponse({"detail": str(exc)}, status_code=422)

    @app.get("/api/health")
    def health() -> dict[str, str]:
        return {"status": "ok"}

    pending_startup = list(startup_files or [])

    @app.get("/api/startup-files")
    def startup() -> dict[str, Any]:
        """Files passed on the command line (e.g. double-click in the file manager), once."""
        results, errors = [], []
        while pending_startup:
            path = pending_startup.pop(0)
            try:
                results.append(_open_path(path))
            except HTTPException as exc:
                errors.append(str(exc.detail))
            except fileio.UnsupportedFile as exc:
                errors.append(str(exc))
        return {"results": results, "errors": errors}

    def _open_path(path: Path) -> dict[str, Any]:
        if not path.is_file():
            raise HTTPException(404, f"“{path.name}” was not found in {path.parent}. It may have "
                                     "been moved, renamed or deleted.")
        try:
            content = path.read_bytes()
        except OSError as exc:
            log.error("Could not read %s: %r", path, exc)
            raise HTTPException(422, f"“{path.name}” could not be opened: "
                                     f"{_os_reason(exc)}") from exc
        result = _read(path.name, content)
        settings.add_recent(path)
        settings.remember_folder(path.parent)
        return {"path": str(path), **result}

    @app.post("/api/open-path")
    def open_path(body: OpenPath) -> dict[str, Any]:
        return _open_path(Path(body.path))

    @app.post("/api/unsaved")
    def unsaved(body: Unsaved) -> dict[str, int]:
        """The UI reports how many documents have unsaved changes (used to confirm quitting)."""
        if on_unsaved:
            on_unsaved(body.count)
        return {"count": body.count}

    @app.post("/api/features/detect")
    def detect_features(body: DetectRequest) -> dict[str, Any]:
        data = body.sequenceData
        sequence = (data.get("sequence") or "").upper()
        if len(sequence) < 20:
            raise HTTPException(422, "The sequence is too short to annotate")
        existing = data.get("features") or []
        if isinstance(existing, dict):
            existing = list(existing.values())
        return common_features.detect_payload(sequence, bool(data.get("circular")),
                                              body.min_identity, existing)

    @app.post("/api/sanger/align")
    def sanger_align(body: SangerRequest) -> dict[str, Any]:
        reference = body.reference.get("sequence") or ""
        if len(reference) < 50:
            raise HTTPException(422, "The reference sequence is too short to align reads to")
        if not body.reads:
            raise HTTPException(422, "Add at least one .ab1 trace")
        reads = [(r.filename, base64.b64decode(r.content_b64)) for r in body.reads]
        return sanger.alignment_payload(
            body.reference.get("name") or "Reference", reference.upper(),
            bool(body.reference.get("circular")), reads, body.cutoff)

    @app.get("/api/gel/ladders")
    def gel_ladders() -> dict[str, Any]:
        return {"ladders": gel.ladders(), "percentages": sorted(gel.RESOLVING_RANGE),
                "default_ladder": gel.DEFAULT_LADDER, "default_percent": gel.DEFAULT_PERCENT}

    @app.post("/api/gel")
    def simulate_gel(body: GelRequest) -> dict[str, Any]:
        from .cloning.base import SequenceDoc
        from .cloning.restriction import digest

        lanes = []
        for lane in body.lanes:
            doc = SequenceDoc(**lane.sequenceData)
            if lane.enzymes:
                sizes = [f["length"] for f in digest(doc, lane.enzymes)]
                supercoiled = False
                title = f"{lane.name} / {' + '.join(lane.enzymes)}"
            else:
                sizes = [len(doc.sequence)]
                supercoiled = bool(doc.circular)
                title = f"{lane.name} / uncut"
            bands = gel.lane_bands(sizes, body.percent, supercoiled)
            lanes.append({"name": title, "bands": [b.as_dict() for b in bands],
                          "sizes": sorted(sizes, reverse=True)})
        return {"lanes": lanes, "ladder": gel.ladder_bands(body.ladder, body.percent),
                "ladder_label": gel.LADDERS[body.ladder]["label"], "percent": body.percent}

    @app.get("/api/recent")
    def recent() -> dict[str, list[str]]:
        return {"paths": settings.recent_files()}

    @app.post("/api/tm")
    def melting_temperature(body: TmRequest) -> dict[str, Any]:
        try:
            presets = tm.all_presets(body.sequence)
        except ValueError as exc:
            raise HTTPException(422, str(exc)) from exc
        default = next(p for p in presets if p["key"] == tm.DEFAULT_PRESET)
        return {"tm": default["tm"], "method": default["description"], "presets": presets}

    @app.post("/api/open-content")
    def open_content(body: OpenContent) -> dict[str, Any]:
        try:
            content = base64.b64decode(body.content_b64, validate=True)
        except (binascii.Error, ValueError) as exc:
            raise HTTPException(422, f"“{body.filename}” did not reach the engine intact; "
                                     "try opening it again.") from exc
        return {"path": None, **_read(body.filename, content)}

    @app.post("/api/open-text")
    def open_text(body: OpenText) -> dict[str, Any]:
        """Pasted text, recognised the same way as file content."""
        name = body.name.strip() or "Pasted sequence"
        return {"path": None, **_read(name, body.text.encode("utf-8"), pasted=True)}

    @app.post("/api/save-path")
    def save_path(body: SavePath) -> dict[str, str]:
        """Save in the format of the path's extension, appending one when there is none."""
        path, fmt = fileio.save_target(body.path, _format(body.format))
        content = fileio.encode(body.sequenceData, fmt)
        try:
            _write_atomically(path, content)
        except OSError as exc:
            raise _save_failed(path, exc) from exc
        if fmt != "fasta":  # a FASTA copy is an export, not a document to reopen
            settings.add_recent(path)
        settings.remember_folder(path.parent)
        return {"path": str(path), "format": fmt, "format_label": fileio.LABELS[fmt]}

    @app.post("/api/save-binary")
    def save_binary(body: SaveBinary) -> dict[str, str]:
        path = Path(body.path)
        try:
            _write_atomically(path, base64.b64decode(body.content_b64))
        except OSError as exc:
            raise _save_failed(path, exc) from exc
        settings.remember_folder(path.parent)
        return {"path": str(path)}

    @app.post("/api/export")
    def export(body: Export) -> dict[str, str]:
        """File content for a browser download (no desktop save dialog)."""
        fmt = _format(body.format)
        if fmt not in fileio.DEFAULT_EXTENSION:
            raise HTTPException(400, f"Unknown export format: {body.format}")
        content = fileio.encode(body.sequenceData, fmt)
        if fmt in fileio.BINARY_FORMATS:
            return {"content_b64": base64.b64encode(content).decode("ascii")}
        return {"content": content.decode("utf-8")}

    if static_dir and static_dir.is_dir():
        app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")

    return app


def _read(filename: str, content: bytes, pasted: bool = False) -> dict[str, Any]:
    """Open result for the UI: documents, format, warnings and whether Save may overwrite."""
    try:
        return fileio.read_file(filename, content, pasted=pasted).as_dict()
    except fileio.UnsupportedFile:
        raise
    except Exception as exc:  # anything unexpected: plain words for the user, details in the log
        log.exception("Unexpected error opening %s", filename)
        what = "The pasted text" if pasted else f"“{Path(filename).name}”"
        raise fileio.UnsupportedFile(
            f"{what} could not be opened because of an unexpected problem. The details were "
            "written to Vettore's log file; please send it along with the file.") from exc
