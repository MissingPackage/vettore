"""The document's History, kept inside the saved file so that lineage survives the file.

A saved GenBank record carries the History in its COMMENT block, after the comment the record
already had. Every reader that understands GenBank keeps COMMENT text and shows it as a
comment, so the block is safe: it says in its first line that other programs can ignore it.

The History itself is JSON (the steps, each with snapshots of the documents the step ran on),
compressed and written as base64 in short lines. Compression is what makes the snapshots
affordable; base64 in short lines is what makes the block survive a round trip through another
program, because a reader that re-wraps or re-indents the comment changes only white space, and
white space is thrown away when the block is read back.

A History that is still too large for a file after compression is written with fewer
generations of snapshots, nearest kept first (the immediate inputs of a step are the ones worth
reopening); if even the steps alone do not fit, the History is left out, because a file that
cannot be written is worse than a file without its lineage.

**Reading a block is reading outside content.** A file may come from anywhere, and a block in it
is whatever its writer put there: compressed data that expands without end, nesting deep enough
to exhaust the stack, or more steps than any laboratory ever performed. So the block is read
under limits at every stage — the base64 itself, what it expands to, and the shape of what it
expands to — and anything that does not fit, or cannot be read, is ignored: a bad block never
keeps a sequence file from opening.
"""

from __future__ import annotations

import base64
import gzip
import json
import logging
import zlib
from typing import Any

log = logging.getLogger(__name__)

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
#: The largest block read back. A little above MAX_BLOB, so a block another version of Vettore
#: wrote is still read, and far below what would cost noticeable memory.
MAX_READ = 1_200_000
#: The most JSON a block may expand to. Compression pays for the snapshots, so the JSON is
#: legitimately many times the block; past this the block is not a laboratory's history.
MAX_JSON = 16_000_000
#: Limits on the shape of a History read from a file: a cloning has a handful of inputs and a
#: construct a few dozen steps, so these are wide enough never to be met by real work.
MAX_STEPS = 500
MAX_PARENTS = 100
#: How deep a History read from a file is followed; the Lineage view shows 20 steps.
MAX_DEPTH = 20


def _collapse(line: str) -> str:
    """A line with its white space squeezed, as another program may have re-indented it."""
    return " ".join(line.split())


def _blob(history: list[Any]) -> str:
    """The History as base64 of its compressed JSON, in lines of ``LINE_LENGTH``."""
    raw = json.dumps(history, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    # mtime=0: the same History gives the same bytes, so saving twice gives the same file.
    text = base64.b64encode(gzip.compress(raw, mtime=0)).decode("ascii")
    return "\n".join(text[i:i + LINE_LENGTH] for i in range(0, len(text), LINE_LENGTH))


def without_snapshots(history: list[Any], keep: int = 0) -> list[Any]:
    """The steps with the copies of their inputs kept for ``keep`` generations and dropped
    below that: the immediate inputs of a step are the ones worth reopening, so they are the
    last to go. Names, paths and the steps themselves always stay."""
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
                if keep > 0:
                    parent = dict(parent)
                else:
                    parent = {k: v for k, v in parent.items()
                              if k not in {"data", "wrongClones"}}
                    parent[TRUNCATED] = True
                if isinstance(parent.get("history"), list):
                    parent["history"] = without_snapshots(parent["history"], keep - 1)
                parents.append(parent)
            step["parents"] = parents
        if keep <= 0:
            # A step's other outcomes carry sequences too, so they go with the snapshots.
            step.pop("wrongClones", None)
        out.append(step)
    return out


def encode(history: list[Any] | None) -> str | None:
    """The block to write into a file, or None when there is no History worth writing."""
    if not history or not isinstance(history, list):
        return None
    blob = _blob(history)
    # Too large: drop the furthest generations of snapshots first, then all of them.
    for keep in (3, 1, 0):
        if len(blob) <= MAX_BLOB:
            break
        blob = _blob(without_snapshots(history, keep))
    if len(blob) > MAX_BLOB:
        return None
    return f"{NOTICE}\n{BEGIN}\n{blob}\n{END}"


def _inflate(raw: bytes) -> bytes | None:
    """The compressed bytes expanded, or None when they expand past ``MAX_JSON`` (a small
    block can describe an endless stream) or are not one gzip member."""
    box = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)      # 16: gzip, not raw deflate
    try:
        out = box.decompress(raw, MAX_JSON)
    except zlib.error:
        return None
    return None if box.unconsumed_tail or not box.eof else out


def _limit(history: Any, depth: int = MAX_DEPTH) -> list[Any]:
    """The History cut to the limits above: too many steps, too many inputs of one step, or a
    tree deeper than the Lineage view follows is a file to read in part, not to choke on."""
    if not isinstance(history, list) or depth <= 0:
        return []
    out = []
    for step in history[:MAX_STEPS]:
        if not isinstance(step, dict):
            continue
        step = dict(step)
        if isinstance(step.get("parents"), list):
            parents = []
            for parent in step["parents"][:MAX_PARENTS]:
                if not isinstance(parent, dict):
                    continue
                parent = dict(parent)
                if "history" in parent:
                    parent["history"] = _limit(parent["history"], depth - 1)
                parents.append(parent)
            step["parents"] = parents
        out.append(step)
    return out


def decode(block: str) -> list[Any] | None:
    """The History of a block as it was written, or None when the block is not one Vettore
    wrote or cannot be read (a bad block never keeps a file from opening)."""
    text = "".join(block.split()) if block else ""
    if not text or len(text) > MAX_READ:
        return None
    try:
        raw = base64.b64decode(text, validate=True)
        expanded = _inflate(raw)
        if expanded is None:
            return None
        history = json.loads(expanded)
    except Exception:
        # Damaged base64, a stream that is not gzip, JSON that is not JSON, nesting deep
        # enough to exhaust the stack (RecursionError, which is not an error about data at
        # all): whatever a block turns out to be, it is only a block, and the sequence file
        # around it must still open. Listing the exceptions instead would mean the next kind
        # nobody thought of keeps a scientist from opening their own plasmid.
        log.debug("the History block of this file could not be read; ignoring it",
                  exc_info=True)
        return None
    if not isinstance(history, list) or not all(isinstance(step, dict) for step in history):
        return None
    return _limit(history)


def add_to_comment(comment: str, history: list[Any] | None) -> str:
    """``comment`` with the History block after it (and any earlier block replaced)."""
    kept = strip(comment)
    block = encode(history)
    if not block:
        return kept
    return f"{kept}\n\n{block}" if kept else block


def _block_lines(comment: str) -> tuple[int, int] | None:
    """Which lines of ``comment`` the History block occupies, notice included, or None.

    Found by reading the lines rather than by one regular expression over the whole comment: a
    comment is outside content, and a pattern that searches across lines can be made to cost
    minutes of processor time by a file that repeats the opening line.
    """
    lines = comment.split("\n")
    begin = end = None
    for i, line in enumerate(lines):
        collapsed = _collapse(line)
        if begin is None:
            if collapsed.startswith(BEGIN):
                begin = i
        elif collapsed == END:
            end = i
            break
    if begin is None or end is None:
        return None
    # The notice belongs to the block, but only where it sits on the line just above it.
    first = begin
    words = []
    while first > 0 and len(words) < len(NOTICE.split()):
        above = _collapse(lines[first - 1])
        if not above:
            break
        words = above.split() + words
        first -= 1
        if " ".join(words) == NOTICE:
            break
    return (first if " ".join(words) == NOTICE else begin), end


def strip(comment: str) -> str:
    """``comment`` without a History block, as the user wrote it."""
    if not comment:
        return ""
    found = _block_lines(comment)
    if not found:
        return comment.strip()
    first, last = found
    lines = comment.split("\n")
    return "\n".join(lines[:first] + lines[last + 1:]).strip()


def from_comment(comment: str) -> list[Any] | None:
    """The History a comment carries, or None."""
    if not comment:
        return None
    found = _block_lines(comment)
    if not found:
        return None
    _first, last = found
    lines = comment.split("\n")
    body_starts = next(i for i in range(len(lines)) if _collapse(lines[i]).startswith(BEGIN)) + 1
    return decode("\n".join(lines[body_starts:last]))


def payload(history: list[Any] | None) -> str | None:
    """The History as one line of base64, for a format that keeps notes rather than a comment
    (the ``.dna`` format). None when there is nothing to write or it does not fit."""
    block = encode(history)
    if not block:
        return None
    lines = block.split("\n")
    return "".join(lines[2:-1])
