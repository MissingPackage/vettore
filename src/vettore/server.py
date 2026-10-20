"""Local HTTP API used by the UI. Bound to 127.0.0.1 and protected by a per-session token."""

from __future__ import annotations

import base64
import secrets
from collections.abc import Callable
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import fileio, gel, sanger, settings, tm
from .cloning.api import router as cloning_router
from .cloning.base import CloningError

TOKEN_HEADER = "X-Vettore-Token"
ALLOWED_HOSTS = {"127.0.0.1", "localhost"}


class OpenPath(BaseModel):
    path: str


class OpenContent(BaseModel):
    filename: str
    content_b64: str


class SavePath(BaseModel):
    path: str
    sequenceData: dict[str, Any]


class SaveBinary(BaseModel):
    path: str
    content_b64: str


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


def _export(data: dict[str, Any], fmt: str) -> str:
    if fmt == "genbank":
        return fileio.write_genbank(data)
    if fmt == "fasta":
        return fileio.write_fasta(data)
    raise HTTPException(400, f"Unknown export format: {fmt}")


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

    @app.exception_handler(fileio.UnsupportedFile)
    async def unsupported(_request: Request, exc: fileio.UnsupportedFile):
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
                results.append({"path": str(path), "documents": _read(path.name, path.read_bytes())})
            except (OSError, fileio.UnsupportedFile) as exc:
                errors.append(str(exc))
        return {"results": results, "errors": errors}

    @app.post("/api/open-path")
    def open_path(body: OpenPath) -> dict[str, Any]:
        path = Path(body.path)
        if not path.is_file():
            raise HTTPException(404, f"File not found: {path}")
        documents = _read(path.name, path.read_bytes())
        settings.add_recent(path)
        return {"path": str(path), "documents": documents}

    @app.post("/api/unsaved")
    def unsaved(body: Unsaved) -> dict[str, int]:
        """The UI reports how many documents have unsaved changes (used to confirm quitting)."""
        if on_unsaved:
            on_unsaved(body.count)
        return {"count": body.count}

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
        content = base64.b64decode(body.content_b64)
        return {"path": None, "documents": _read(body.filename, content)}

    @app.post("/api/save-path")
    def save_path(body: SavePath) -> dict[str, str]:
        path = Path(body.path)
        fmt = "fasta" if path.suffix.lower() in {".fa", ".fasta", ".fna", ".fas"} else "genbank"
        path.write_text(_export(body.sequenceData, fmt), encoding="utf-8")
        if fmt == "genbank":
            settings.add_recent(path)
        return {"path": str(path)}

    @app.post("/api/save-binary")
    def save_binary(body: SaveBinary) -> dict[str, str]:
        path = Path(body.path)
        path.write_bytes(base64.b64decode(body.content_b64))
        return {"path": str(path)}

    @app.post("/api/export")
    def export(body: Export) -> dict[str, str]:
        return {"content": _export(body.sequenceData, body.format)}

    if static_dir and static_dir.is_dir():
        app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")

    return app


def _read(filename: str, content: bytes) -> list[dict[str, Any]]:
    try:
        return fileio.read_documents(filename, content)
    except fileio.UnsupportedFile:
        raise
    except Exception as exc:  # parser errors from Biopython are heterogeneous
        raise fileio.UnsupportedFile(f"Could not read {filename}: {exc}") from exc
