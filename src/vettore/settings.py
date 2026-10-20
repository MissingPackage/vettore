"""User-level persistent state (recent files, last folder) under ~/.vettore or $VETTORE_HOME."""

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


def _state_file() -> Path:
    return home() / "state.json"


def recent_files() -> list[str]:
    try:
        paths = json.loads(_recent_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    if not isinstance(paths, list):
        return []
    return [p for p in paths if isinstance(p, str) and Path(p).is_file()][:MAX_RECENT]


def add_recent(path: Path) -> None:
    entry = str(path)
    paths = [entry] + [p for p in recent_files() if p != entry]
    try:
        _recent_file().write_text(json.dumps(paths[:MAX_RECENT], indent=1), encoding="utf-8")
    except OSError:  # the list is a convenience; never fail an open or a save because of it
        pass


def _state() -> dict:
    try:
        state = json.loads(_state_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return state if isinstance(state, dict) else {}


def last_folder() -> str:
    """The folder of the last file opened or saved, so file dialogs start there."""
    folder = _state().get("last_folder")
    return folder if isinstance(folder, str) and Path(folder).is_dir() else ""


def remember_folder(folder: Path) -> None:
    state = _state()
    state["last_folder"] = str(folder)
    try:
        _state_file().write_text(json.dumps(state, indent=1), encoding="utf-8")
    except OSError:
        pass
