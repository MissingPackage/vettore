"""Opening the files colleagues send: every format, the wrong extension, and other tools' quirks.

Reference values come from the public records themselves (NCBI J01749.1 / ENA J01749 for pBR322,
NCBI L09137.2 for pUC19c) and from Biopython's own tests of its .xdna/.gck test files.
"""

from pathlib import Path

import pytest

from vettore import fileio

DATA = Path(__file__).parent / "data"


def read(name: str, as_name: str | None = None) -> fileio.ReadResult:
    return fileio.read_file(as_name or name, (DATA / name).read_bytes())


def only(result: fileio.ReadResult) -> dict:
    [doc] = result.documents
    return doc


def plain(doc: dict) -> str:
    return doc["sequence"].upper()


PBR322 = plain(only(read("exchange-pBR322.embl")))
PUC19 = plain(read("exchange-plasmids.fasta").documents[0])


def test_reference_sequences_agree_across_formats():
    # pBR322 from ENA (EMBL), NCBI (GenBank) and NCBI (FASTA) is the same 4361-bp molecule.
    assert len(PBR322) == 4361
    assert plain(only(read("exchange-quirks-header.gb"))) == PBR322
    assert plain(read("exchange-plasmids.fasta").documents[1]) == PBR322
    assert len(PUC19) == 2686 and PUC19.startswith("TCGCGCGTTTCGGTGATGAC")


def test_embl_opens_with_features_and_topology():
    result = read("exchange-pBR322.embl")
    doc = only(result)
    assert result.format == "embl" and result.label == "EMBL"
    assert doc["name"] == "exchange-pBR322"  # a single record takes the file name
    assert doc["circular"] is True
    assert len(doc["features"]) == 48  # 50 in the record, minus its two source features
    bla = [f for f in doc["features"] if f["type"] == "CDS" and f["strand"] == -1]
    assert any((f["start"], f["end"]) == (3292, 4152) for f in bla)  # complement(3293..4153)
    assert result.writable is False  # Vettore does not write EMBL


@pytest.mark.parametrize(("name", "length", "circular", "features"), [
    # Bio.SeqIO test suite expectations for its sample files.
    ("exchange-sample-a.xdna", 1000, False,
     [("promoter", "FeatureA", 49, 149, 1), ("misc_binding", "FeatureB", 499, 699, -1)]),
    ("exchange-sample-b.xdna", 1000, True,
     [("rep_origin", "FeatureA", 160, 240, 1), ("terminator", "FeatureB", 399, 749, -1)]),
    ("exchange-artificial.gck", 12, False,
     [("misc_feature", "FeatureA", 2, 5, 1), ("CDS", "FeatureB", 7, 10, -1)]),
])
def test_binary_formats_from_other_programs(name, length, circular, features):
    doc = only(read(name))
    assert len(doc["sequence"]) == length
    assert doc["circular"] is circular
    got = sorted((f["type"], f["name"], f["start"], f["end"], f["strand"]) for f in doc["features"])
    assert got == sorted(features)


def test_multi_record_fasta_opens_one_document_per_record():
    result = read("exchange-plasmids.fasta")
    assert [d["name"] for d in result.documents] == ["L09137.2", "J01749.1", "X06403.1"]
    assert [len(d["sequence"]) for d in result.documents] == [2686, 4361, 4245]
    assert result.writable is False  # saving one tab must not overwrite the other two


def test_records_with_the_same_name_get_distinct_names():
    text = (DATA / "MT612431.1.gb").read_bytes()
    result = fileio.read_file("two.gb", text + b"\n" + text)
    assert [d["name"] for d in result.documents] == ["MT612431", "MT612431 (2)"]
    result = fileio.read_file("x", b">a\nACGT\n>a\nGGCC\n>a\nTTTT\n", pasted=True)
    assert [d["name"] for d in result.documents] == ["a", "a (2)", "a (3)"]


@pytest.mark.parametrize("name", ["exchange-pUC19-numbered.seq", "exchange-pUC19-raw.txt"])
def test_plain_sequence_text(name):
    result = read(name)
    assert result.format == fileio.RAW
    doc = only(result)
    assert plain(doc) == PUC19
    assert doc["circular"] is False and doc["features"] == []


@pytest.mark.parametrize("source", ["exchange-pBR322.embl", "exchange-plasmids.fasta",
                                    "exchange-sample-b.xdna", "exchange-artificial.gck",
                                    "MT612431.1.gb", "pFA-KanMX4.dna", "3730.ab1",
                                    "exchange-pUC19-raw.txt"])
@pytest.mark.parametrize("as_name", ["no-extension", "from-a-colleague.txt", "wrong.seq",
                                     "wrong.gb", "wrong.fasta"])
def test_format_is_recognised_from_content_not_extension(source, as_name):
    expected = read(source)
    result = read(source, as_name)
    assert result.format == expected.format
    assert [plain(d) for d in result.documents] == [plain(d) for d in expected.documents]


def test_writable_only_when_saving_back_loses_nothing():
    assert read("MT612431.1.gb").writable is True
    assert read("MT612431.1.gb", "colleague.txt").writable is False  # GenBank under .txt
    assert read("pFA-KanMX4.dna").writable is False
    assert read("3730.ab1").writable is False
    assert read("exchange-pUC19-raw.txt").writable is False
    assert fileio.read_file("one.fasta", b">one\nACGT\n").writable is True


def test_genbank_quirks_of_other_tools_open():
    # Long LOCUS name, no topology or date, Windows line endings, a cp1252 byte, no final //.
    result = read("exchange-quirks-header.gb")
    doc = only(result)
    assert plain(doc) == PBR322
    assert doc["circular"] is False  # no topology on the LOCUS line: linear, as GenBank defines
    notes = [f["notes"].get("note", [""])[0] for f in doc["features"]]
    assert "non-ASCII test: 37 °C, 5 µl" in notes
    assert len(doc["features"]) == 48
    assert result.warnings == []


def test_feature_quirks_are_repaired_or_reported():
    result = read("exchange-quirks-features.gb")
    doc = only(result)
    by_name: dict[str, list[dict]] = {}
    for f in doc["features"]:
        by_name.setdefault(f["name"], []).append(f)
    # Both ways of writing a feature across the origin give OVE's start > end form.
    for name in ("test region across origin (join)", "test region across origin (range)"):
        [f] = by_name[name]
        assert (f["start"], f["end"], f["strand"]) == (4299, 49, 1)
    # A feature running past the end is cut at the end of the sequence.
    [past] = by_name["test region past the end"]
    assert (past["start"], past["end"]) == (4299, 4360)
    # Unusual qualifiers are kept as notes; duplicate names stay as separate features.
    odd = by_name["test region with odd qualifiers"]
    assert len(odd) == 2
    assert odd[0]["notes"]["note"] == ['a "quoted" word']
    assert odd[0]["notes"]["pseudo"] == [""]
    assert odd[0]["notes"]["vntifkey"] == ["21"]
    assert len(by_name["tet"]) > 1
    # Features the editor cannot place are left out, and the user is told which and why.
    for name in ("test site between bases", "test region beyond the end",
                 "test region in another record"):
        assert name not in by_name
    assert len(result.warnings) == 2
    assert "3 features could not be placed" in result.warnings[0]
    assert "“test site between bases” (between two bases)" in result.warnings[0]
    assert "“test region in another record” (it refers to another sequence)" in result.warnings[0]
    assert "“test region past the end”" in result.warnings[1]
    # Saving over the file would silently drop those features: Save makes a copy instead.
    assert result.writable is False


def test_genbank_without_final_terminator_and_with_old_mac_line_endings():
    text = (DATA / "MT612431.1.gb").read_text().rstrip().removesuffix("//").replace("\n", "\r")
    doc = only(fileio.read_file("x.gb", text.encode()))
    assert len(doc["sequence"]) == 7204


@pytest.mark.parametrize("encode", [
    lambda t: t.encode("utf-8-sig"),          # BOM
    lambda t: t.encode("utf-16"),             # Windows "Unicode" text
    lambda t: ("\n\n" + t).encode("utf-8"),   # blank lines before the first record
])
def test_text_encodings_and_leading_blank_lines(encode):
    doc = only(fileio.read_file("x.fasta", encode(">pUC19c\nTCGCGCGTTT\nCGGTGATGAC\n")))
    assert doc["sequence"] == "TCGCGCGTTTCGGTGATGAC"


def test_bad_character_is_named_with_line_and_column():
    with pytest.raises(fileio.UnsupportedFile) as err:
        fileio.read_file("notes.txt", b"ACGTACGT\nACGTXACGT\n")
    message = str(err.value)
    assert message.startswith("“notes.txt” could not be read as a DNA sequence")
    assert "line 2, column 5 has “X”" in message


def test_name_on_the_first_line_gets_a_hint():
    with pytest.raises(fileio.UnsupportedFile, match="start it with “>”"):
        fileio.read_file("x", b"pUC19 sequence\nTCGCGCGTTT\n", pasted=True)


def test_protein_and_alignment_gaps_are_explained():
    with pytest.raises(fileio.UnsupportedFile, match="looks like a protein sequence"):
        fileio.read_file("p.fasta", b">p\nMKVLAAGIVGLLLQEEFFPPWWLLEEQQ\n")
    with pytest.raises(fileio.UnsupportedFile, match="“-” at position 3.*alignment gaps"):
        fileio.read_file("a.fasta", b">aln\nAC--GT\n")


def test_good_records_open_when_one_record_is_bad():
    result = fileio.read_file("mixed.fasta", b">good\nACGT\n>bad\nAC-GT\n")
    assert [d["name"] for d in result.documents] == ["good"]
    assert "1 of 2 sequences could not be opened" in result.warnings[0]
    assert result.writable is False


@pytest.mark.parametrize(("name", "content", "expected"), [
    ("notes.docx", b"PK\x03\x04\x00\x00\x00binary", "is not a sequence file Vettore can read"),
    ("empty.gb", b"", "“empty.gb” is empty."),
    ("broken.dna", b"junk", "looks like a .dna file but could not be read"),
    ("broken.ab1", b"\x00\x01junk", "Ask for the trace file again"),
    ("broken.gb", b"LOCUS       x 10 bp DNA\nORIGIN\n        1 acgt\nnot genbank", ""),
])
def test_unreadable_files_get_plain_messages(name, content, expected):
    with pytest.raises(fileio.UnsupportedFile) as err:
        fileio.read_file(name, content)
    message = str(err.value)
    assert expected in message
    assert "Traceback" not in message and "Error" not in message


def test_pasted_text_keeps_the_record_name_or_takes_the_given_one():
    genbank = (DATA / "MT612431.1.gb").read_bytes()
    assert only(fileio.read_file("Pasted sequence", genbank, pasted=True))["name"] == "MT612431"
    raw = fileio.read_file("My insert", b"1 acgtacgtac gt\n", pasted=True)
    assert only(raw)["name"] == "My insert"
    assert only(raw)["sequence"] == "acgtacgtacgt"
    assert raw.writable is False


def test_save_target_follows_the_extension():
    assert fileio.save_target("/x/a.gb") == (Path("/x/a.gb"), "genbank")
    assert fileio.save_target("/x/a.GBK") == (Path("/x/a.GBK"), "genbank")
    assert fileio.save_target("/x/a.fa", "genbank") == (Path("/x/a.fa"), "fasta")
    assert fileio.save_target("/x/a.dna") == (Path("/x/a.dna"), fileio.DNA_FILE_FORMAT)
    # No extension Vettore writes: the requested format's extension is appended.
    assert fileio.save_target("/x/a") == (Path("/x/a.gb"), "genbank")
    assert fileio.save_target("/x/a", "fasta") == (Path("/x/a.fasta"), "fasta")
    assert fileio.save_target("/x/pAAV.CMV.GFP") == (Path("/x/pAAV.CMV.GFP.gb"), "genbank")
    assert fileio.save_target("/x/a.txt") == (Path("/x/a.txt.gb"), "genbank")


def test_fasta_keeps_the_whole_name():
    text = fileio.write_fasta({"name": "pAAV CAG GFP WPRE bGH", "sequence": "ACGT",
                               "circular": False, "features": [], "primers": []})
    assert text.splitlines()[0].startswith(">pAAV_CAG_GFP_WPRE_bGH")
