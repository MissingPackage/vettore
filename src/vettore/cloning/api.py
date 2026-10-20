"""HTTP routes for cloning simulations, mounted under /api/cloning by vettore.server."""

from __future__ import annotations

from fastapi import APIRouter

from . import gateway_ta, golden_gate, overlap, restriction

router = APIRouter(prefix="/api/cloning")
for module in (restriction, overlap, golden_gate, gateway_ta):
    router.include_router(module.router)
