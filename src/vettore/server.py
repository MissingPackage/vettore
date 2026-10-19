"""Local HTTP API used by the UI. Bound to 127.0.0.1 and protected by a per-session token."""

from __future__ import annotations

import base64
import secrets
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import fileio

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
               startup_files: list[Path] | None = None) -> FastAPI:
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
        return {"path": str(path), "documents": _read(path.name, path.read_bytes())}

    @app.post("/api/open-content")
    def open_content(body: OpenContent) -> dict[str, Any]:
        content = base64.b64decode(body.content_b64)
        return {"path": None, "documents": _read(body.filename, content)}

    @app.post("/api/save-path")
    def save_path(body: SavePath) -> dict[str, str]:
        path = Path(body.path)
        fmt = "fasta" if path.suffix.lower() in {".fa", ".fasta", ".fna", ".fas"} else "genbank"
        path.write_text(_export(body.sequenceData, fmt), encoding="utf-8")
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
