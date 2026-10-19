"""User-level persistent state (recent files) under ~/.vettore or $VETTORE_HOME."""

from __future__ import annotations

import json
import os
from pathlib import Path

MAX_RECENT = 10


def home() -> Path:
    path = Path(os.environ.get("VETTORE_HOME") or Path.home() / ".vettore")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _recent_file() -> Path:
    return home() / "recent.json"


def recent_files() -> list[str]:
    try:
        paths = json.loads(_recent_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return [p for p in paths if isinstance(p, str) and Path(p).is_file()][:MAX_RECENT]


def add_recent(path: Path) -> None:
    entry = str(path)
    paths = [entry] + [p for p in recent_files() if p != entry]
    _recent_file().write_text(json.dumps(paths[:MAX_RECENT], indent=1), encoding="utf-8")
