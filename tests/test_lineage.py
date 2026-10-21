"""The document's History inside the file it is saved in (S6.1).

What is checked: the History comes back exactly as it was written, it survives another program
reading and writing the file (Biopython here), a History too large for a file is written
without its snapshots rather than making the file unreasonable, a damaged block never keeps a
file from opening, and the comment the file already had is kept.
"""

from __future__ import annotations

import base64
import gzip
import io
import json
import random
import time

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


# ---- a block is outside content ----------------------------------------------------------
# A file can come from anywhere. These are the blocks a file should never be able to write:
# each one made Vettore hang, exhaust memory or refuse to open the sequence at all.


def file_with_block(body: str, name: str = "evil.gb") -> bytes:
    """A GenBank file whose COMMENT carries ``body`` as its History block."""
    plain = fileio.write_genbank(document())
    block = "".join(f"            {body[i:i + 60]}\n" for i in range(0, len(body), 60))
    comment = (f"COMMENT     {lineage.NOTICE}\n            {lineage.BEGIN}\n"
               f"{block}            {lineage.END}\n")
    return plain.replace("FEATURES", comment + "FEATURES", 1).encode()


def block_of(payload: object) -> str:
    """``payload`` as the base64 of its compressed JSON, however large it is."""
    raw = json.dumps(payload, separators=(",", ":")).encode()
    return base64.b64encode(gzip.compress(raw)).decode()


def test_a_block_that_expands_without_end_is_refused(monkeypatch) -> None:
    """Compression means a small block can describe a very large JSON ("a gzip bomb"): the
    limit is on what it expands to, not only on the block."""
    monkeypatch.setattr(lineage, "MAX_JSON", 10_000)
    blob = block_of([{"time": "t", "action": "A" * 100_000}])
    assert len(blob) <= lineage.MAX_READ          # the block itself looks ordinary
    assert lineage.decode(blob) is None
    [doc] = fileio.read_documents("evil.gb", file_with_block(blob))
    assert doc["sequence"]                        # and the sequence still opens
    assert "history" not in doc


def test_a_block_larger_than_a_file_should_hold_is_not_read() -> None:
    assert lineage.decode("A" * (lineage.MAX_READ + 4)) is None


def test_nesting_deep_enough_to_exhaust_the_stack_does_not_keep_the_file_closed() -> None:
    """Nesting raises RecursionError, which is not an error about data: the file still opens."""
    deep = '[{"action":"x","parents":[' * 20_000 + '{"name":"leaf"}' + "]}]" * 20_000
    blob = base64.b64encode(gzip.compress(deep.encode())).decode()
    assert lineage.decode(blob) is None
    [doc] = fileio.read_documents("deep.gb", file_with_block(blob))
    assert doc["sequence"] and "history" not in doc


def test_a_comment_that_repeats_the_opening_line_costs_no_time() -> None:
    """Searching across lines for the block can be made quadratic by a file that repeats the
    opening line and never closes it, so the block is found by reading lines."""
    comment = "\n".join([lineage.BEGIN] * 32_000)
    started = time.monotonic()
    assert lineage.from_comment(comment) is None
    assert lineage.strip(comment) == comment.strip()
    assert time.monotonic() - started < 1.0


def test_more_steps_or_inputs_than_any_laboratory_performed_are_cut() -> None:
    crowd = [{"time": "t", "action": "Gibson assembly",
              "parents": [{"name": f"fragment {i}"} for i in range(150_000)]}]
    [step] = lineage.decode(block_of(crowd))
    assert len(step["parents"]) == lineage.MAX_PARENTS
    many = [{"time": "t", "action": f"step {i}"} for i in range(lineage.MAX_STEPS + 50)]
    assert len(lineage.decode(block_of(many))) == lineage.MAX_STEPS


def test_a_tree_deeper_than_the_view_follows_is_cut() -> None:
    history: list = [{"time": "t", "action": "Opened"}]
    for i in range(40):
        history = [{"time": "t", "action": f"step {i}",
                    "parents": [{"name": f"input {i}", "history": history}]}]
    depth, step = 0, lineage.decode(block_of(history))
    while step:
        depth += 1
        step = step[0].get("parents", [{}])[0].get("history")
    assert depth == lineage.MAX_DEPTH


def test_a_comment_that_merely_quotes_the_notice_is_left_alone() -> None:
    """The notice is part of the block, not a sentence to delete wherever it appears: a lab's
    own note that quotes it (or a file saved by a Vettore whose block was stripped) stays."""
    note = f"Our SOP says: {lineage.NOTICE}\nKeep clone 7."
    assert lineage.strip(note) == note
    # Through a file, GenBank's own line wrapping may move a line break; the words all stay.
    assert " ".join(reopen(document(comment=note))["comment"].split()) == " ".join(note.split())
    # With a real block after it, the notice goes with the block and the note stays whole.
    doc = reopen(document(comment=note, history=HISTORY))
    assert " ".join(doc["comment"].split()) == " ".join(note.split())
    assert doc["history"] == HISTORY


def test_the_nearest_snapshots_are_the_last_to_be_dropped(monkeypatch) -> None:
    """When the History does not fit, the inputs of the newest step are the ones worth
    keeping: those are what a scientist reopens."""
    monkeypatch.setattr(lineage, "MAX_BLOB", 30_000)
    history: list = [{"time": "t", "action": "Opened"}]
    for i in range(6):
        history = [{"time": "t", "action": f"assembly {i}", "parents": [
            {"name": f"input {i}",
             "data": {"name": f"input {i}",
                      "sequence": "".join(random.Random(i).choices("ACGT", k=20_000)),
                      "circular": False, "features": [], "primers": []},
             "history": history}]}]
    doc = reopen(document(history=history))
    [step] = doc["history"]
    [nearest] = step["parents"]
    assert "data" in nearest                      # the immediate input is still there
    furthest = nearest
    while furthest.get("history") and furthest["history"][0].get("parents"):
        furthest = furthest["history"][0]["parents"][0]
    assert furthest is not nearest
    assert "data" not in furthest and furthest[lineage.TRUNCATED] is True


def test_a_records_structured_comment_survives_a_save() -> None:
    """The ##Assembly-Data## tables NCBI records carry are provenance: a save must not drop
    them, with or without a History of our own next to them."""
    gb = ("LOCUS       test                      20 bp    DNA     linear   SYN 21-OCT-2026\n"
          "DEFINITION  t.\n"
          "COMMENT     A plain note.\n"
          "            ##Assembly-Data-START##\n"
          "            Assembly Method       :: SPAdes v. 3.9\n"
          "            Sequencing Technology :: Illumina\n"
          "            ##Assembly-Data-END##\n"
          "FEATURES             Location/Qualifiers\n"
          "ORIGIN\n"
          "        1 acgtacgtac gtacgtacgt\n"
          "//\n")
    assembly = {"Assembly-Data": {"Assembly Method": "SPAdes v. 3.9",
                                  "Sequencing Technology": "Illumina"}}
    [doc] = fileio.read_documents("t.gb", gb.encode())
    assert doc["comment"] == "A plain note."
    assert doc["structured_comment"] == assembly

    doc["history"] = HISTORY
    [again] = fileio.read_documents("t.gb", fileio.write_genbank(doc).encode())
    assert again["comment"] == "A plain note."
    assert again["structured_comment"] == assembly
    assert again["history"] == HISTORY
