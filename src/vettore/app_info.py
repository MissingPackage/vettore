"""Application-level routes under /api/app: version, sample constructs, preferences."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from typing import Any

from fastapi import APIRouter, HTTPException

router = APIRouter(prefix="/api/app")


def app_version() -> str:
    try:
        return version("vettore")
    except PackageNotFoundError:
        return "dev"


@router.get("/about")
def about() -> dict[str, Any]:
    return {"version": app_version()}


@router.get("/samples")
def samples() -> dict[str, Any]:
    return {"samples": []}


@router.get("/samples/{sample_id}")
def open_sample(sample_id: str) -> dict[str, Any]:
    raise HTTPException(404, f"Unknown sample: {sample_id}")
