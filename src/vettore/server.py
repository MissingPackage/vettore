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
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import alignment, app_info, fileio, gel, primer_design, sanger, settings, tm
from . import features as common_features
from .cloning.api import router as cloning_router
from .cloning.base import CloningError, dna_sequence
from .desktop import associations
from .desktop.instance import HAND_OVER_WAIT, TOKEN_HEADER, Inbox

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


class HandOver(BaseModel):
    """Files from a second start of Vettore (a double-click while this one runs)."""

    paths: list[str] = []


class FileAssociationsChange(BaseModel):
    """Open sequence files with Vettore on double-click (``enabled``), or undo it.
    ``dry_run`` only lists what would be written."""

    enabled: bool
    dry_run: bool = False


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


# What the user was doing, by route, for the messages about problems inside Vettore.
ACTIONS: tuple[tuple[str, str], ...] = (
    ("/api/open-", "Opening the file"),
    ("/api/startup-files", "Opening the file"),
    ("/api/save-", "Saving the file"),
    ("/api/export", "Saving the file"),
    ("/api/features/", "Detecting common features"),
    ("/api/sanger/", "Aligning the reads"),
    ("/api/gel", "Simulating the gel"),
    ("/api/tm", "Calculating the Tm"),
    ("/api/alignment/", "Aligning the sequences"),
    ("/api/primer-design/", "Designing primers"),
    ("/api/cloning/", "The simulation"),
    ("/api/app/samples", "Opening the sample"),
    ("/api/unsaved", "Counting the documents with unsaved changes"),
    ("/api/hand-over", "Opening the file"),
    ("/api/desktop/file-associations", "Changing which program opens sequence files"),
)


def _action(path: str) -> str:
    return next((label for prefix, label in ACTIONS if path.startswith(prefix)), "That action")


def unexpected_message(path: str) -> str:
    """One sentence for a problem inside Vettore; the traceback goes to the log file."""
    return (f"{_action(path)} failed because of an unexpected problem inside Vettore. The "
            f"details were written to the log file {app_info.log_path()}; please send it to "
            "the developers with a note of what you were doing.")


def malformed_message(path: str) -> str:
    """One sentence for a request the engine cannot read (a fault of the window, not of the
    user's data); the field-by-field report goes to the log file."""
    return (f"{_action(path)} failed because Vettore's window sent a request its engine could "
            f"not read. Try again; if it happens again, please send the log file "
            f"{app_info.log_path()} to the developers.")


def _save_target(raw: str) -> Path:
    """The path the user chose, or an error a person can act on: no name, a folder, or a name
    with a character no file system accepts."""
    path = Path(raw)
    if not raw.strip() or not path.name:
        raise HTTPException(422, "No file name was given. Choose where to save the file and "
                                 "type a name for it.")
    if "\x00" in raw:
        raise HTTPException(422, f"“{path.name}” cannot be used as a file name. Choose another "
                                 "name.")
    if path.is_dir():
        raise HTTPException(422, f"“{path.name}” is a folder. Choose a file name inside it and "
                                 "save again.")
    return path


def _write_or_explain(path: Path, content: bytes) -> None:
    try:
        _write_atomically(path, content)
    except OSError as exc:
        raise _save_failed(path, exc) from exc


def _save_failed(path: Path, exc: OSError) -> HTTPException:
    log.error("Could not save %s: %r", path, exc)
    return HTTPException(422, f"“{path.name}” could not be saved in {path.parent}: "
                              f"{_os_reason(exc)}")


def create_app(token: str, static_dir: Path | None,
               startup_files: list[Path] | None = None,
               on_unsaved: Callable[[int], None] | None = None,
               on_hand_over: Callable[[], None] | None = None,
               hand_over_wait: float = HAND_OVER_WAIT) -> FastAPI:
    """The API and the UI's files. ``on_unsaved`` hears the count of unsaved documents;
    ``on_hand_over`` is called when a second start hands files over (to bring the window
    forward and tell the page); ``hand_over_wait`` is how long the page has to take them."""
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

    @app.exception_handler(RequestValidationError)
    async def malformed_request(request: Request, exc: RequestValidationError):
        # A malformed request is a fault of the UI, never of the user's data: one sentence for
        # the dialog, the field-by-field report for the log.
        log.error("Malformed request to %s: %s", request.url.path, exc.errors())
        return JSONResponse({"detail": malformed_message(request.url.path)},
                            status_code=422)

    @app.exception_handler(Exception)
    async def unexpected(request: Request, exc: Exception):
        log.error("Unexpected error in %s %s", request.method, request.url.path, exc_info=exc)
        return JSONResponse({"detail": unexpected_message(request.url.path)}, status_code=500)

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

    inbox = Inbox(startup_files or [])

    @app.get("/api/startup-files")
    def startup() -> dict[str, Any]:
        """Files waiting to be opened, each served once: those on the command line (a
        double-click in the file manager), then those handed over by later starts. The page
        asks every half second."""
        results, errors = [], []
        for path in inbox.take():
            try:
                results.append(_open_path(path))
            except HTTPException as exc:
                errors.append(str(exc.detail))
            except fileio.UnsupportedFile as exc:
                errors.append(str(exc))
        return {"results": results, "errors": errors}

    def _open_path(path: Path) -> dict[str, Any]:
        if not path.name:
            raise HTTPException(422, "No file was chosen. Choose a sequence file to open.")
        try:
            is_file, is_folder = path.is_file(), path.is_dir()
        except (OSError, ValueError) as exc:   # e.g. a folder Vettore may not look into
            log.error("Could not look at %s: %r", path, exc)
            reason = (_os_reason(exc) if isinstance(exc, OSError)
                      else "the name contains a character no file name can have.")
            raise HTTPException(422, f"“{path.name}” could not be opened: {reason}") from exc
        if is_folder:
            raise HTTPException(422, f"“{path.name}” is a folder, not a file. Open it in the "
                                     "file dialog and choose a sequence file inside it.")
        if not is_file:
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

    @app.post("/api/hand-over")
    def hand_over(body: HandOver) -> dict[str, bool]:
        """Queue files from a second start and wait until the page has taken them.
        ``delivered`` is false when it did not within ``hand_over_wait`` seconds (a frozen or
        closed page): the files are withdrawn and the caller opens its own window."""
        ticket = inbox.put(Path(p) for p in body.paths)
        if on_hand_over:
            on_hand_over()
        return {"delivered": inbox.wait_or_withdraw(ticket, hand_over_wait)}

    @app.get("/api/desktop/file-associations")
    def file_associations() -> dict[str, Any]:
        """Whether sequence files open with this copy of Vettore on double-click, and whether
        it can offer that at all (only the packaged application can)."""
        target = associations.for_this_system()
        return target.status() if target else associations.unavailable_status()

    @app.post("/api/desktop/file-associations")
    def set_file_associations(body: FileAssociationsChange) -> dict[str, Any]:
        target = associations.for_this_system()
        if target is None:
            raise HTTPException(409, associations.unavailable_status()["reason"])
        if body.dry_run:
            return {**target.status(), "plan": target.plan()}
        try:
            if body.enabled:
                target.enable()
            else:
                target.disable()
        except OSError as exc:
            log.error("Could not change the file associations: %r", exc)
            # enable() undoes a partial change itself; disable() may stop half way.
            outcome = ("Nothing was changed." if body.enabled else
                       "Some sequence files may still open with Vettore: try again.")
            raise HTTPException(
                422, "Vettore could not change which program opens sequence files"
                     f"{': ' + exc.strerror if exc.strerror else ''}. {outcome} If it happens "
                     f"again, please send the log file {app_info.log_path()} to the "
                     "developers.") from exc
        return target.status()

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
        sequence = data.get("sequence")
        sequence = sequence.upper() if isinstance(sequence, str) else ""
        if len(sequence) < 20:
            raise HTTPException(422, "The sequence is too short to annotate: Detect Common "
                                     "Features needs at least 20 bp.")
        return common_features.detect_payload(sequence, bool(data.get("circular")),
                                              body.min_identity, data.get("features"))

    @app.post("/api/sanger/align")
    def sanger_align(body: SangerRequest) -> dict[str, Any]:
        name = str(body.reference.get("name") or "Reference")
        reference = body.reference.get("sequence")
        reference = reference if isinstance(reference, str) else ""
        if len(reference) < 50:
            raise HTTPException(422, "The reference sequence is too short to align reads to: "
                                     "it needs at least 50 bp.")
        try:
            reference = dna_sequence(name, reference).upper()
        except CloningError as exc:
            raise HTTPException(422, str(exc)) from exc
        if not body.reads:
            raise HTTPException(422, "Add at least one .ab1 trace")
        reads = []
        for read in body.reads:
            try:
                reads.append((read.filename, base64.b64decode(read.content_b64, validate=True)))
            except (binascii.Error, ValueError) as exc:
                raise HTTPException(422, f"“{read.filename}” did not reach the engine intact; "
                                         "add it again.") from exc
        return sanger.alignment_payload(name, reference, bool(body.reference.get("circular")),
                                        reads, body.cutoff)

    @app.get("/api/gel/ladders")
    def gel_ladders() -> dict[str, Any]:
        return {"ladders": gel.ladders(), "percentages": sorted(gel.RESOLVING_RANGE),
                "default_ladder": gel.DEFAULT_LADDER, "default_percent": gel.DEFAULT_PERCENT}

    @app.post("/api/gel")
    def simulate_gel(body: GelRequest) -> dict[str, Any]:
        """Each lane digested, or drawn uncut with a note (``lanes[i].note``) when its enzyme
        does not cut: a non-cutter is a normal control, not a reason to fail the gel."""
        if body.ladder not in gel.LADDERS:
            raise HTTPException(422, gel.unknown_ladder(body.ladder))
        return gel.simulate([(lane.name, lane.sequenceData, lane.enzymes) for lane in body.lanes],
                            body.ladder, body.percent)

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
    def save_path(body: SavePath) -> dict[str, Any]:
        """Save in the format of the path's extension, appending one when there is none.
        ``warnings`` (only when there are any) says what the file leaves out."""
        path, fmt = fileio.save_target(_save_target(body.path), _format(body.format))
        content = fileio.encode(body.sequenceData, fmt)
        _write_or_explain(path, content)
        if fmt != "fasta":  # a FASTA copy is an export, not a document to reopen
            settings.add_recent(path)
        settings.remember_folder(path.parent)
        result: dict[str, Any] = {"path": str(path), "format": fmt,
                                  "format_label": fileio.LABELS[fmt]}
        if warnings := fileio.save_warnings(body.sequenceData, fmt):
            result["warnings"] = warnings
        return result

    @app.post("/api/save-binary")
    def save_binary(body: SaveBinary) -> dict[str, str]:
        path = _save_target(body.path)
        try:
            content = base64.b64decode(body.content_b64, validate=True)
        except (binascii.Error, ValueError) as exc:
            raise HTTPException(422, f"“{path.name}” did not reach the engine intact; save it "
                                     "again.") from exc
        _write_or_explain(path, content)
        settings.remember_folder(path.parent)
        return {"path": str(path)}

    @app.post("/api/export")
    def export(body: Export) -> dict[str, Any]:
        """File content for a browser download (no desktop save dialog)."""
        fmt = _format(body.format)
        if fmt not in fileio.DEFAULT_EXTENSION:
            raise HTTPException(400, f"Vettore cannot export the “{body.format}” format; "
                                     "choose GenBank, FASTA or .dna.")
        content = fileio.encode(body.sequenceData, fmt)
        result: dict[str, Any] = (
            {"content_b64": base64.b64encode(content).decode("ascii")}
            if fmt in fileio.BINARY_FORMATS else {"content": content.decode("utf-8")})
        if warnings := fileio.save_warnings(body.sequenceData, fmt):
            result["warnings"] = warnings
        return result

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
