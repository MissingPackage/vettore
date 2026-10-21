"""The document's History, kept inside the saved file so that lineage survives the file.

A saved GenBank record carries the History in its COMMENT block, after the comment the record
already had. Every reader that understands GenBank keeps COMMENT text and shows it as a
comment, so the block is safe: it says in its first line that other programs can ignore it.

The History itself is JSON (the steps, each with snapshots of the documents the step ran on),
compressed and written as base64 in short lines. Compression is what makes the snapshots
affordable; base64 in short lines is what makes the block survive a round trip through another
program, because a reader that re-wraps or re-indents the comment changes only white space, and
white space is thrown away when the block is read back.

A History that is still too large for a file after compression is written without the snapshots
(the steps are kept), and if even that does not fit, it is left out: a file that cannot be
written is worse than a file without its lineage.
"""

from __future__ import annotations

import base64
import binascii
import gzip
import json
import re
from typing import Any

#: First line of the block: what the block is, and that it can be ignored.
NOTICE = "Vettore History (the lineage of this construct), safe to ignore."
BEGIN = "##vettore-lineage-v1 (gzip, base64)"
END = "##vettore-lineage-end"
#: Dropped snapshots are recorded, so the Lineage view can say why a snapshot is missing.
TRUNCATED = "snapshots_dropped"

#: Characters of base64 per line: 60 plus GenBank's 12-column indent stays inside the 79
#: columns a GenBank line has, so a reader that re-wraps the comment has nothing to re-wrap.
#: Re-wrapping would be harmless anyway, because white space is thrown away when the block is
#: read back.
LINE_LENGTH = 60
#: The largest History block written into a file, in characters of base64. A plasmid file is
#: tens of kilobytes; a megabyte of lineage would make the file something else.
MAX_BLOB = 1_000_000

_BLOCK = re.compile(rf"^[ \t]*{re.escape(BEGIN)}\s*(?P<body>.*?)^[ \t]*{re.escape(END)}[ \t]*$",
                    re.DOTALL | re.MULTILINE)
# The notice as another program may have re-wrapped it: its words, any white space between.
_NOTICE_LINE = re.compile(r"[ \t]*" + r"\s+".join(re.escape(word) for word in NOTICE.split())
                          + r"[ \t]*\n?")
_BASE64_ONLY = re.compile(r"[^A-Za-z0-9+/=]")


def _blob(history: list[Any]) -> str:
    """The History as base64 of its compressed JSON, in lines of ``LINE_LENGTH``."""
    raw = json.dumps(history, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    # mtime=0: the same History gives the same bytes, so saving twice gives the same file.
    text = base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")
    return "\n".join(text[i:i + LINE_LENGTH] for i in range(0, len(text), LINE_LENGTH))


def without_snapshots(history: list[Any]) -> list[Any]:
    """The same steps without the copies of their input documents (names and paths stay)."""
    out = []
    for step in history:
        if not isinstance(step, dict):
            continue
        step = dict(step)
        if isinstance(step.get("parents"), list):
            parents = []
            for parent in step["parents"]:
                if not isinstance(parent, dict):
                    continue
                parent = {k: v for k, v in parent.items() if k != "data"}
                if isinstance(parent.get("history"), list):
                    parent["history"] = without_snapshots(parent["history"])
                parent[TRUNCATED] = True
                parents.append(parent)
            step["parents"] = parents
        out.append(step)
    return out


def encode(history: list[Any] | None) -> str | None:
    """The block to write into a file, or None when there is no History worth writing."""
    if not history or not isinstance(history, list):
        return None
    blob = _blob(history)
    if len(blob) > MAX_BLOB:
        blob = _blob(without_snapshots(history))
    if len(blob) > MAX_BLOB:
        return None
    return f"{NOTICE}\n{BEGIN}\n{blob}\n{END}"


def decode(block: str) -> list[Any] | None:
    """The History of a block as it was written, or None when the block is not one Vettore
    wrote or cannot be read (a damaged block never keeps a file from opening)."""
    text = _BASE64_ONLY.sub("", block)
    if not text:
        return None
    try:
        history = json.loads(gzip.decompress(base64.b64decode(text, validate=True)))
    except (binascii.Error, ValueError, OSError, gzip.BadGzipFile, UnicodeDecodeError):
        return None
    if not isinstance(history, list) or not all(isinstance(step, dict) for step in history):
        return None
    return history


def add_to_comment(comment: str, history: list[Any] | None) -> str:
    """``comment`` with the History block after it (and any earlier block replaced)."""
    kept = strip(comment)
    block = encode(history)
    if not block:
        return kept
    return f"{kept}\n\n{block}" if kept else block


def strip(comment: str) -> str:
    """``comment`` without a History block, as the user wrote it."""
    without = _BLOCK.sub("", comment or "")
    return _NOTICE_LINE.sub("", without).strip()


def from_comment(comment: str) -> list[Any] | None:
    """The History a comment carries, or None."""
    if not comment:
        return None
    match = _BLOCK.search(comment)
    return decode(match.group("body")) if match else None


def payload(history: list[Any] | None) -> str | None:
    """The History as one line of base64, for a format that keeps notes rather than a comment
    (the ``.dna`` format). None when there is nothing to write or it does not fit."""
    if not history or not isinstance(history, list):
        return None
    blob = _blob(history).replace("\n", "")
    if len(blob) > MAX_BLOB:
        blob = _blob(without_snapshots(history)).replace("\n", "")
    return blob if len(blob) <= MAX_BLOB else None
