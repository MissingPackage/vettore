"""The document's History inside the file it is saved in (S6.1).

What is checked: the History comes back exactly as it was written, it survives another program
reading and writing the file (Biopython here), a History too large for a file is written
without its snapshots rather than making the file unreasonable, a damaged block never keeps a
file from opening, and the comment the file already had is kept.
"""

from __future__ import annotations

import io
import random

import pytest
from Bio import SeqIO

from vettore import fileio, lineage

HISTORY = [
    {"time": "21/10/2026, 09:12:00", "action": "Opened", "detail": "backbone.gb (GenBank)"},
    {
        "time": "21/10/2026, 10:30:00",
        "action": "Restriction cloning",
        "detail": "EcoRI + BamHI. From: backbone, insert.",
        "parents": [
            {
                "name": "backbone",
                "path": "/lab/plasmids/backbone.gb",
                "data": {"name": "backbone", "sequence": "ACGT" * 200, "circular": True,
                         "features": [], "primers": []},
                "history": [{"time": "20/10/2026, 08:00:00", "action": "Opened",
                             "detail": "backbone.gb"}],
            },
            {
                "name": "insert",
                "path": None,
                "data": {"name": "insert", "sequence": "TTGCA" * 40, "circular": False,
                         "features": [], "primers": []},
                "history": [],
            },
        ],
    },
]


def document(**extra) -> dict:
    return {"name": "pLineage", "sequence": "ACGT" * 150, "circular": True,
            "features": [], "primers": [], **extra}


def reopen(data: dict, name: str = "pLineage.gb") -> dict:
    [doc] = fileio.read_documents(name, fileio.encode(data, fileio.SAVE_FORMATS[".gb"]))
    return doc


def test_genbank_keeps_the_history() -> None:
    doc = reopen(document(history=HISTORY))
    assert doc["history"] == HISTORY


def test_the_block_says_it_can_be_ignored() -> None:
    text = fileio.write_genbank(document(history=HISTORY))
    assert "COMMENT" in text
    assert lineage.NOTICE.split(",")[0] in text.replace("\n", " ").replace("            ", "")
    # The sequence and the features are written as they were before, after the comment.
    assert "ORIGIN" in text and text.index("COMMENT") < text.index("ORIGIN")


def test_a_file_without_a_history_gets_no_block() -> None:
    text = fileio.write_genbank(document())
    assert "vettore-lineage" not in text
    assert "history" not in reopen(document())


def test_the_comment_of_the_file_is_kept_and_stays_the_users() -> None:
    doc = reopen(document(comment="Clone 3 grew badly; use clone 7.", history=HISTORY))
    assert doc["comment"] == "Clone 3 grew badly; use clone 7."
    assert doc["history"] == HISTORY


def test_another_program_can_read_and_write_the_file() -> None:
    """Biopython reads the file and writes it again: the lineage survives as the comment does."""
    written = fileio.write_genbank(document(comment="A note.", history=HISTORY))
    record = SeqIO.read(io.StringIO(written), "genbank")
    out = io.StringIO()
    SeqIO.write(record, out, "genbank")
    [doc] = fileio.read_documents("pLineage.gb", out.getvalue().encode())
    assert doc["history"] == HISTORY
    assert doc["comment"] == "A note."


def test_re_wrapped_lines_are_still_read() -> None:
    """A reader that re-indents or re-wraps the comment changes only white space."""
    block = lineage.encode(HISTORY)
    assert block is not None
    squeezed = "   " + block.replace("\n", "\n      ")
    assert lineage.from_comment(squeezed) == HISTORY


def test_a_damaged_block_is_ignored_and_the_file_still_opens() -> None:
    text = fileio.write_genbank(document(history=HISTORY))
    broken = text.replace(lineage.BEGIN, lineage.BEGIN).replace("H4sI", "ZZZZ", 1)
    [doc] = fileio.read_documents("pLineage.gb", broken.encode())
    assert doc["sequence"]
    assert "history" not in doc


def test_a_block_from_another_version_is_ignored() -> None:
    assert lineage.from_comment("##vettore-lineage-v9 (something else)\nAAAA\n") is None
    assert lineage.from_comment("A plain comment.") is None


def test_a_history_that_is_too_large_is_written_without_its_snapshots() -> None:
    """A snapshot per step would make the file something other than a plasmid file."""
    huge = [{"time": "t", "action": "Gibson assembly", "detail": "d",
             "parents": [{"name": f"fragment {i}", "path": None,
                          "data": {"name": f"fragment {i}",
                                   # Random bases: there is nothing for compression to save.
                                   "sequence": "".join(random.Random(i).choices("ACGT",
                                                                                k=500_000)),
                                   "circular": False, "features": [], "primers": []},
                          "history": []}]}
            for i in range(1, 9)]
    assert lineage.encode(huge) == lineage.encode(lineage.without_snapshots(huge))
    doc = reopen(document(history=huge))
    assert [step["action"] for step in doc["history"]] == ["Gibson assembly"] * 8
    parents = [parent for step in doc["history"] for parent in step["parents"]]
    assert [p["name"] for p in parents] == [f"fragment {i}" for i in range(1, 9)]
    assert all("data" not in parent and parent[lineage.TRUNCATED] is True for parent in parents)


def test_the_same_history_gives_the_same_bytes() -> None:
    """Saving twice without a change must not make the file look changed."""
    assert lineage.encode(HISTORY) == lineage.encode(HISTORY)


@pytest.mark.parametrize("history", [None, [], "not a history", 42])
def test_nothing_is_written_for_what_is_not_a_history(history: object) -> None:
    assert lineage.encode(history) is None  # type: ignore[arg-type]
    assert lineage.payload(history) is None  # type: ignore[arg-type]


def test_a_dna_file_keeps_the_history_in_its_notes() -> None:
    content = fileio.encode(document(history=HISTORY), fileio.DNA_FILE_FORMAT)
    [doc] = fileio.read_documents("pLineage.dna", content)
    assert doc["history"] == HISTORY


def test_fasta_keeps_the_sequence_only() -> None:
    text = fileio.write_fasta(document(history=HISTORY))
    assert "vettore-lineage" not in text
