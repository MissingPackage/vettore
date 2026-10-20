"""Multiple sequence alignment of open documents. Routes under /api/alignment."""

from __future__ import annotations

from fastapi import APIRouter

router = APIRouter(prefix="/api/alignment")
