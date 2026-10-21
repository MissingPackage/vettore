"""Lab primer library (S5.3): the shared .csv, the one binding rule, products and sequencing
primers, and the HTTP routes. Binding is checked on sequences built so the annealing footprint
is known exactly."""


import pytest
from Bio.Seq import Seq
from fastapi.testclient import TestClient

from vettore import primer_library as pl
from vettore import settings
from vettore.server import TOKEN_HEADER, create_app

TOKEN = "test-token"


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


# A 24-mer: the last 15 (the anchor) are "TGACCTGATTAACGT"; the first 9 can carry tails/mismatches.
PRIMER = "GCATTCGAA" + "TGACCTGATTAACGT"
assert len(PRIMER) == 24


def linear(before: str, insert: str, after: str) -> str:
    return (before + insert + after).upper()


# ---- the one binding rule ----------------------------------------------------------------

def test_forward_exact_match_reports_the_footprint_and_strand():
    seq = linear("AAAACCCCGGGGTTTT", PRIMER, "ACACACACAC")
    [b] = pl.find_bindings(PRIMER, seq, circular=False)
    assert b.strand == 1 and b.mismatches == 0 and b.tail == 0
    assert b.footprint == len(PRIMER)
    assert seq[b.start:b.end + 1] == PRIMER


def test_reverse_match_is_found_on_the_other_strand():
    seq = linear("AAAACCCCGGGGTTTT", rc(PRIMER), "ACACACACAC")
    [b] = pl.find_bindings(PRIMER, seq, circular=False)
    assert b.strand == -1 and b.mismatches == 0 and b.tail == 0
    # Like every primer in the editor: start is the lower coordinate, on both strands.
    assert b.start < b.end and b.three_prime == b.start
    assert rc(seq[b.start:b.end + 1]) == PRIMER


def test_a_5_prime_tail_does_not_count_against_the_binding():
    # The whole primer with a non-matching 5' tail: only the body (without the tail) is on the map.
    body = PRIMER
    tailed = "TTTTT" + body             # five bases that are not on the template
    seq = linear("GAGAGAGAGA", body, "CTCTCTCTCT")
    [b] = pl.find_bindings(tailed, seq, circular=False)
    assert b.tail == 5 and b.footprint == len(body) and b.mismatches == 0
    assert seq[b.start:b.end + 1] == body


def test_up_to_two_mismatches_inside_the_footprint_still_bind():
    # Two interior mismatches (positions 3 and 6, both 5' of the anchor); the 5' end still
    # matches, so the whole primer anneals with two mismatches counted.
    templ = "GCAATCAAA" + "TGACCTGATTAACGT"           # PRIMER with positions 3 and 6 changed
    assert templ[-15:] == PRIMER[-15:]
    assert sum(a != b for a, b in zip(PRIMER, templ)) == 2
    assert templ[0] == PRIMER[0]                     # the 5' base matches
    seq = linear("GGGGGGGG", templ, "AAAAAAAA")
    [b] = pl.find_bindings(PRIMER, seq, circular=False)
    assert b.strand == 1 and b.mismatches == 2 and b.footprint == len(PRIMER) and b.tail == 0


def test_a_run_of_mismatches_at_the_5_prime_end_is_the_tail():
    # Three mismatches at the very 5' end: they anneal to nothing, so they are the tail, not
    # scored as mismatches.
    templ = "TTG" + PRIMER[3:]                        # positions 0,1,2 differ from PRIMER
    assert sum(a != b for a, b in zip(PRIMER, templ)) == 3
    seq = linear("GGGGGGGG", templ, "AAAAAAAA")
    [b] = pl.find_bindings(PRIMER, seq, circular=False)
    assert b.mismatches == 0 and b.tail == 3 and b.footprint == len(PRIMER) - 3


def test_a_changed_anchor_base_means_no_binding():
    templ = PRIMER[:-1] + ("A" if PRIMER[-1] != "A" else "C")   # last (3') base changed
    seq = linear("GGGGGGGG", templ, "AAAAAAAA")
    assert pl.find_bindings(PRIMER, seq, circular=False) == []


def test_a_primer_across_the_origin_binds_once_on_a_circular_template():
    # Split the primer so half sits at the end and half at the start of the circle.
    head, tail = PRIMER[:14], PRIMER[14:]
    seq = (tail + "ACGTACGTACGTACGT" + head).upper()   # primer wraps: ...head | tail...
    bindings = pl.find_bindings(PRIMER, seq, circular=True)
    assert len(bindings) == 1
    b = bindings[0]
    assert b.strand == 1 and b.footprint == len(PRIMER) and b.start > b.end  # wraps the origin


def test_primers_shorter_than_the_anchor_are_not_searched():
    # A 10-mer would match a plasmid by chance many times: it is not drawn anywhere.
    short = "ACGTACGTAC"
    assert pl.find_bindings(short, linear("TTTT", short, "GGGG"), circular=False) == []


# ---- products and sequencing primers -----------------------------------------------------

def test_a_forward_and_reverse_library_primer_make_a_product():
    fwd = "ACGTTGCATTCGAATGACCTG"
    rev_site = "GGTTACCGATTACGCATTGCA"
    insert = "ATATATATATATATATAT"
    seq = linear(fwd, insert, rc(rev_site))
    primers = [pl.Primer("F", fwd), pl.Primer("R", rev_site)]
    products = pl.amplified_products(primers, seq, circular=False)
    assert products and products[0]["forward"] == "F" and products[0]["reverse"] == "R"
    assert products[0]["size"] == len(seq)


def test_product_over_region_reports_the_pair_that_already_covers_it():
    fwd = "ACGTTGCATTCGAATGACCTG"
    rev_site = "GGTTACCGATTACGCATTGCA"
    middle = "CCCCGGGGAAAATTTTCCCC"
    seq = linear(fwd, middle, rc(rev_site))
    primers = [pl.Primer("KanR_F", fwd), pl.Primer("KanR_R", rev_site)]
    region_start = len(fwd) + 2
    region_end = len(fwd) + len(middle) - 2
    hit = pl.product_over_region(primers, seq, False, region_start, region_end)
    assert hit and hit["forward"] == "KanR_F" and hit["reverse"] == "KanR_R"


def test_sequencing_primers_read_across_the_region_from_the_right_distance():
    fwd = "TTGCATTCGAATGACCTGAT"              # 20 nt
    spacer = "A" * 70                          # ~70 bp upstream of the region (within 50-100)
    region = "GGGGCCCCGGGGCCCCGGGG"
    seq = linear(fwd, spacer + region, "TTTTTTTT")
    primers = [pl.Primer("Seq1", fwd)]
    region_start = len(fwd) + len(spacer)
    region_end = region_start + len(region) - 1
    reads = pl.sequencing_primers(primers, seq, False, region_start, region_end)
    assert len(reads) == 1 and reads[0]["name"] == "Seq1" and reads[0]["strand"] == 1
    assert reads[0]["reads_from"] == region_start - (len(fwd) - 1)   # from the 3' end


# ---- the shared .csv ---------------------------------------------------------------------

def test_csv_with_a_header_round_trips():
    text = "name,sequence,notes,location,added\nM13F,GTAAAACGACGGCCAGT,universal,box 3 A4,2024-01-02\n"
    primers, rejected = pl.parse_csv(text)
    assert rejected == [] and len(primers) == 1
    p = primers[0]
    assert p.name == "M13F" and p.sequence == "GTAAAACGACGGCCAGT" and p.location == "box 3 A4"
    again, _ = pl.parse_csv(pl.merged_text(None, primers)[0])
    assert again == primers


def test_pasted_two_columns_without_a_header_are_read():
    pasted = "M13F\tGTAAAACGACGGCCAGT\nM13R\tCAGGAAACAGCTATGAC\n"   # tab-separated, Excel paste
    primers, rejected = pl.parse_csv(pasted)
    assert rejected == [] and [p.name for p in primers] == ["M13F", "M13R"]


def test_rows_with_non_dna_are_listed_and_left_out():
    text = "name,sequence\nGood,ACGTACGT\nBad,ACGZXY\nEmpty,\n"
    primers, rejected = pl.parse_csv(text)
    assert [p.name for p in primers] == ["Good"]
    assert [r["name"] for r in rejected] == ["Bad", "Empty"]
    assert "not DNA" in rejected[0]["reason"] and "no sequence" in rejected[1]["reason"]


def test_a_repeated_sequence_is_reported_not_duplicated():
    text = "name,sequence\nA,ACGTACGTAC\nB,ACGTACGTAC\n"
    primers, rejected = pl.parse_csv(text)
    assert len(primers) == 1 and len(rejected) == 1 and "same oligo as “A”" in rejected[0]["reason"]


def test_u_is_read_as_t_and_whitespace_is_dropped():
    assert pl.clean_sequence(" acg uac gt ") == "ACGTACGT"
    with pytest.raises(pl.PrimerLibraryError):
        pl.clean_sequence("ACGT123")


# ---- HTTP --------------------------------------------------------------------------------

@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path))
    settings.remember_primer_library_path(str(tmp_path / "lab-primers.csv"))
    app = create_app(TOKEN, static_dir=None)
    return TestClient(app, base_url="http://127.0.0.1")


def headers():
    return {TOKEN_HEADER: TOKEN}


def test_save_then_get_round_trips_through_the_file(client):
    rows = [{"name": "M13F", "sequence": "GTAAAACGACGGCCAGT", "notes": "", "location": "box 1", "added": "2024-01-02"}]
    saved = client.post("/api/primer-library/save", json={"added": rows}, headers=headers())
    assert saved.status_code == 200 and saved.json()["added"] == 1
    got = client.get("/api/primer-library", headers=headers())
    assert got.status_code == 200
    body = got.json()
    assert body["exists"] and body["primers"][0]["name"] == "M13F"
    assert body["primers"][0]["tm"] is not None


def test_bindings_route_finds_a_stored_primer_on_a_sequence(client):
    seq = ("AAAACCCC" + PRIMER + "GGGGTTTT").upper()
    client.post("/api/primer-library/save",
                json={"added": [{"name": "P1", "sequence": PRIMER}]}, headers=headers())
    r = client.post("/api/primer-library/bindings",
                    json={"sequenceData": {"sequence": seq, "circular": False}}, headers=headers())
    assert r.status_code == 200
    primers = r.json()["primers"]
    assert len(primers) == 1 and primers[0]["bindings"][0]["strand"] == 1


def test_parse_route_reports_rejected_rows_without_saving(client):
    r = client.post("/api/primer-library/parse",
                    json={"text": "name,sequence\nGood,ACGTACGT\nBad,ACGXYZ\n"}, headers=headers())
    assert r.status_code == 200
    body = r.json()
    assert len(body["primers"]) == 1 and len(body["rejected"]) == 1 and body["today"]
    # Nothing was written: the library is still empty.
    assert client.get("/api/primer-library", headers=headers()).json()["primers"] == []


def test_bindings_route_uses_the_dialogs_unsaved_list_when_given(client):
    # Nothing saved, but the dialog sends its working list: the new primer binds at once.
    seq = ("AAAACCCC" + PRIMER + "GGGGTTTT").upper()
    r = client.post("/api/primer-library/bindings",
                    json={"sequenceData": {"sequence": seq, "circular": False},
                          "primers": [{"name": "Unsaved", "sequence": PRIMER}]}, headers=headers())
    assert r.status_code == 200
    assert [p["name"] for p in r.json()["primers"]] == ["Unsaved"]


# ---- review findings: coordinates, Excel files, degenerate and repetitive primers ---------

def test_a_reverse_primer_across_the_origin_keeps_the_editors_convention():
    seq = (rc(PRIMER)[10:] + "ACGT" * 10 + rc(PRIMER)[:10]).upper()
    [b] = pl.find_bindings(PRIMER, seq, circular=True)
    assert b.strand == -1 and b.start > b.end          # start > end only when crossing the origin
    assert b.start == len(seq) - 10 and b.end == len(PRIMER) - 10 - 1
    assert b.three_prime == b.start


def test_semicolon_files_from_excel_set_to_italian_are_read():
    text = "Nome;Sequenza;Note\nM13F;GTAAAACGACGGCCAGT;box 3\nM13R;CAGGAAACAGCTATGAC;box 3\n"
    primers, rejected = pl.parse_csv(text)
    assert rejected == [] and [p.name for p in primers] == ["M13F", "M13R"]
    assert primers[0].notes == "box 3"


def test_usual_header_names_are_recognised():
    text = "Primer name,Sequence (5'-3'),Comments,Box\nT7,TAATACGACTCACTATAGGG,promoter,A1\n"
    [p], rejected = pl.parse_csv(text)
    assert rejected == [] and (p.name, p.notes, p.location) == ("T7", "promoter", "A1")


def test_a_comma_in_the_notes_of_pasted_columns_stays_in_the_notes():
    pasted = "M13F\tGTAAAACGACGGCCAGT\tuniversal, -20, box 3, A4\n"
    [p], rejected = pl.parse_csv(pasted)
    assert rejected == [] and p.notes == "universal, -20, box 3, A4"


@pytest.fixture()
def library_file(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path))
    path = tmp_path / "shared.csv"
    settings.remember_primer_library_path(str(path))
    return path


def test_a_windows_code_page_file_is_read(library_file):
    library_file.write_bytes("name;sequence;notes\nM13F;GTAAAACGACGGCCAGT;5 µl, 4 °C\n".encode("cp1252"))
    lib = pl.read_library()
    assert lib["primers"][0]["notes"] == "5 µl, 4 °C"


def test_saving_keeps_unreadable_rows_and_the_files_separator(library_file):
    library_file.write_text(
        "name;sequence;notes\nM13F;GTAAAACGACGGCCAGT;box 3\nInosine;ACGTIACGT;keep me\n"
        "M13F copy;GTAAAACGACGGCCAGT;box 7\n", encoding="utf-8")
    lib = pl.read_library()
    assert [p["name"] for p in lib["primers"]] == ["M13F"] and len(lib["rejected"]) == 2
    new = pl.Primer("KanR_F", "ATGATTGAACAAGATGGATTGCAC")
    result = pl.write_library([new], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "box 3"})])
    assert result["kept"] == 2
    raw = library_file.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf")             # a BOM, so Excel shows µ and accents
    text = raw.decode("utf-8-sig")
    assert "Inosine;ACGTIACGT;keep me" in text and "M13F copy;GTAAAACGACGGCCAGT;box 7" in text
    assert "KanR_F;ATGATTGAACAAGATGGATTGCAC" in text
    assert [p["name"] for p in pl.read_library()["primers"]] == ["M13F", "KanR_F"]


def test_a_degenerate_primer_binds_every_base_it_stands_for():
    degenerate = "GCNTTCGAA" + "TGACCTGRTTAACGT"     # N at 3, R (A/G) inside the anchor
    seq = linear("GGGGGGGG", PRIMER, "AAAAAAAA")
    [b] = pl.find_bindings(degenerate, seq, circular=False)
    assert b.mismatches == 0 and b.tail == 0 and b.footprint == len(degenerate)
    [r] = pl.find_bindings(degenerate, rc(seq), circular=False)
    assert r.strand == -1 and r.mismatches == 0


def test_a_primer_in_a_repeat_counts_once_and_too_many_sites_are_not_drawn():
    oligo_dt = "T" * 20
    seq = "GC" * 20 + "A" * 40 + "GC" * 20
    # 21 overlapping 20-mers fit the 40-A stretch: one site, the stretch is one repeat.
    [site] = pl.find_bindings(oligo_dt, seq, circular=False)
    assert site.positions == 26        # every place its 3′ 15 nt fit the stretch
    repeated = "".join(PRIMER + "ACACACACAC" for _ in range(12))
    [hit] = pl.library_on_sequence([pl.Primer("P", PRIMER)], repeated, circular=False)
    assert hit["sites"] == 12 and hit["bindings"] == []


def test_a_palindromic_primer_is_not_a_product_with_itself():
    palindrome = "GCGAATTCGCGAATTCGC"
    assert palindrome == rc(palindrome)
    seq = linear("ATATATATATATAT", palindrome, "TATATATATATATA")
    assert pl.amplified_products([pl.Primer("P", palindrome)], seq, circular=False) == []


def test_a_region_across_the_origin_is_not_inside_a_product_that_misses_it():
    fwd, rev_site = "ACGTTGCATTCGAATGACCTG", "GGTTACCGATTACGCATTGCA"
    seq = ("T" * 100 + fwd + "C" * 700 + rc(rev_site) + "G" * 1158).upper()   # 2,000 bp
    primers = [pl.Primer("F", fwd), pl.Primer("R", rev_site)]
    [product] = pl.amplified_products(primers, seq, circular=True)
    assert (product["product_start"], product["product_end"]) == (100, 841)
    assert pl.product_over_region(primers, seq, True, 300, 600) is not None
    assert pl.product_over_region(primers, seq, True, 800, 150) is None    # crosses the origin


def test_a_sequencing_primer_a_few_hundred_bases_upstream_is_listed():
    fwd = "TTGCATTCGAATGACCTGAT"
    seq = linear(fwd, "A" * 300 + "GGGGCCCCGGGGCCCC", "T" * 600)
    start = len(fwd) + 300
    [read] = pl.sequencing_primers([pl.Primer("Seq", fwd)], seq, False, start, start + 15)
    assert read["reads_from"] == 301
    # Too close to read (the first ~40 bases after a primer are not readable), or too far.
    assert pl.sequencing_primers([pl.Primer("Seq", fwd)], seq, False, len(fwd) + 10, len(fwd) + 20) == []
    assert pl.sequencing_primers([pl.Primer("Seq", fwd)], seq, False, start, start + 700) == []


# ---- second review: the shared file is kept as the lab wrote it ---------------------------

def test_extra_columns_and_the_header_survive_a_save(library_file):
    original = ("Nome;Sequenza;Tm;Lunghezza;Chi;Note\r\n"
                "M13F;GTAAAACGACGGCCAGT;53;17;Anna;box 3\r\n")
    library_file.write_bytes(original.encode("utf-8"))
    result = pl.write_library([pl.Primer("KanR_F", "ATGATTGAACAAGATGGATTGCAC")],
                              [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "box 4"})])
    text = library_file.read_bytes().decode("utf-8-sig")
    lines = text.split("\r\n")
    assert lines[0] == "Nome;Sequenza;Tm;Lunghezza;Chi;Note"            # header as written
    assert lines[1] == "M13F;GTAAAACGACGGCCAGT;53;17;Anna;box 4"        # Tm, length, who stay
    assert lines[2] == "KanR_F;ATGATTGAACAAGATGGATTGCAC"
    assert (result["added"], result["changed"]) == (1, 1)


def test_a_name_column_wins_over_an_id_column():
    [p], _ = pl.parse_csv("ID,Name,Sequence\n17,M13F,GTAAAACGACGGCCAGT\n")
    assert p.name == "M13F"


def test_rows_stay_in_their_columns_whatever_the_column_order(library_file):
    library_file.write_text("sequence,name,notes\nGTAAAACGACGGCCAGT,M13F,u\nnot-dna,Broken,fix me\n",
                            encoding="utf-8")
    before = pl.read_library()
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "v"})])
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "u"})])
    after = pl.read_library()
    assert [p["name"] for p in after["primers"]] == ["M13F"]   # the name never becomes a sequence
    assert after["rejected"] == before["rejected"]
    assert "not-dna,Broken,fix me" in library_file.read_text(encoding="utf-8-sig")


def test_an_unchanged_file_is_written_back_byte_for_byte(library_file):
    body = ('name;sequence;notes\r\n"M13F";gtaaaacgacggccagt;"5 µl; 4 °C"\r\n'
            "Inosine;ACGTIACGT;keep\r\n")
    library_file.write_bytes(b"\xef\xbb\xbf" + body.encode("utf-8"))
    # The dialog sends every field of every primer as it has them: nothing differs, nothing moves.
    edits = [pl.Change(p["sequence"], {k: p[k] for k in pl.EDITABLE}, {k: p[k] for k in pl.EDITABLE})
             for p in pl.read_library()["primers"]]
    pl.write_library([], edits)
    assert library_file.read_bytes() == b"\xef\xbb\xbf" + body.encode("utf-8")   # lowercase kept


def test_excels_sep_line_names_the_separator_and_stays(library_file):
    library_file.write_text("sep=;\nname;sequence;notes\nM13F;GTAAAACGACGGCCAGT;a, b\n", encoding="utf-8")
    lib = pl.read_library()
    assert lib["delimiter"] == ";" and lib["primers"][0]["notes"] == "a, b"
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "a, b; c"})])
    assert library_file.read_text(encoding="utf-8-sig").startswith("sep=;")


def test_a_tab_inside_a_note_does_not_make_the_file_tab_separated():
    [p], rejected = pl.parse_csv("name,sequence,notes\nM13F,GTAAAACGACGGCCAGT,box\t3\n")
    assert rejected == [] and p.notes == "box\t3"


@pytest.mark.parametrize("data,words", [
    (b"PK\x03\x04rest-of-a-workbook", "Excel workbook"),
    ("name,sequence\nM13F,GTAAAACGACGGCCAGT\n".encode("utf-16"), None),
    (b"LOCUS       pUC19   2686 bp    DNA     circular\nORIGIN\n        1 tcgcgcgttt\n//\n",
     "not a primer library"),
    (b"\x00\x01\x02binary", "not a text file"),
])
def test_files_that_are_not_a_primer_list_are_refused_and_left_alone(library_file, data, words):
    library_file.write_bytes(data)
    if words is None:                      # Excel's "Unicode text" (UTF-16) is read and kept so
        assert pl.read_library()["primers"][0]["name"] == "M13F"
        pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "x"})])
        assert library_file.read_bytes().startswith((b"\xff\xfe", b"\xfe\xff"))
        return
    with pytest.raises(pl.PrimerLibraryError, match=words):
        pl.read_library()
    with pytest.raises(pl.PrimerLibraryError):
        pl.write_library([pl.Primer("M13F", "GTAAAACGACGGCCAGT")])
    assert library_file.read_bytes() == data


def test_old_mac_line_endings_are_read():
    primers, _ = pl.parse_csv("name,sequence\rM13F,GTAAAACGACGGCCAGT\rM13R,CAGGAAACAGCTATGAC\r")
    assert [p.name for p in primers] == ["M13F", "M13R"]


def test_a_library_path_must_be_a_csv(client):
    r = client.post("/api/primer-library/path", json={"path": "/lab/pUC19.gb"}, headers=headers())
    assert r.status_code == 422 and ".csv" in r.json()["detail"]


def test_a_nul_in_pasted_text_is_a_rejected_row_not_a_crash(client):
    r = client.post("/api/primer-library/parse", json={"text": "name,seq\nA,AC\x00GT\n"}, headers=headers())
    assert r.status_code == 200 and len(r.json()["rejected"]) == 1


def test_a_utf8_bom_before_windows_code_page_text_keeps_the_header(library_file):
    library_file.write_bytes(b"\xef\xbb\xbf" + "Nome;Sequenza;Note\nM13F;GTAAAACGACGGCCAGT;5 µl\n".encode("cp1252"))
    [p] = pl.read_library()["primers"]
    assert (p["name"], p["notes"]) == ("M13F", "5 µl")


def test_saving_applies_only_this_windows_changes(client):
    library_file = settings.primer_library_path()
    library_file.write_text("name;sequence;notes\nM13F;GTAAAACGACGGCCAGT;orig\n"
                            "M13R;CAGGAAACAGCTATGAC;\nT7;TAATACGACTCACTATAGGG;\n", encoding="utf-8")
    # B edits M13F's notes and removes T7; A, who loaded the same file, only adds SP6.
    b = client.post("/api/primer-library/save", headers=headers(), json={
        "edited": [{"sequence": "GTAAAACGACGGCCAGT", "was": {"name": "M13F"}, "notes": "B's note"}],
        "removed": [{"sequence": "TAATACGACTCACTATAGGG", "was": {"name": "T7"}}]})
    assert b.status_code == 200
    a = client.post("/api/primer-library/save", headers=headers(), json={
        "added": [{"name": "SP6", "sequence": "ATTTAGGTGACACTATAG"}]})
    body = a.json()
    assert [p["name"] for p in body["primers"]] == ["M13F", "M13R", "SP6"]   # T7 stays deleted
    assert body["primers"][0]["notes"] == "B's note"                          # B's edit stays
    # An edit to a primer someone deleted meanwhile is not saved, and is counted.
    c = client.post("/api/primer-library/save", headers=headers(), json={
        "edited": [{"sequence": "TAATACGACTCACTATAGGG", "was": {"name": "T7"}, "notes": "too late"}]})
    assert c.json()["gone"] == 1 and "T7" not in [p["name"] for p in c.json()["primers"]]



def test_a_dinucleotide_repeat_is_one_site_and_long_n_primers_stay_fast():
    import random
    import time
    rnd = random.Random(5)
    seq = "".join(rnd.choice("ACGT") for _ in range(7000)) + "CA" * 60 + \
        "".join(rnd.choice("ACGT") for _ in range(8000))
    sites = pl.find_bindings("CA" * 10, seq, circular=True)
    assert [b.strand for b in sites] == [1]           # one site, not 41...
    assert sites[0].positions > 10                    # ...where it primes at many positions
    [hit] = pl.library_on_sequence([pl.Primer("CA", "CA" * 10)], seq, circular=True)
    assert hit["sites"] > pl.MAX_SITES and hit["bindings"] == []    # a repeat, not drawn
    t = time.perf_counter()
    assert len(pl.find_bindings("N" * 18, seq, circular=True)) == 2    # all one big "repeat"
    assert time.perf_counter() - t < 2


def test_a_repeat_across_the_origin_is_one_site():
    seq = "A" * 30 + "GCGTCGATCGTAGCTAGCATCGAT" * 10 + "A" * 30
    assert len([b for b in pl.find_bindings("T" * 20, seq, circular=True) if b.strand == -1]) == 1


def test_primers_with_overlapping_footprints_still_make_a_product():
    rnd_left = "GATCGATTGCAAGCTTGGCACTGGCCGTCGTTTTACAACG"
    insert = "ACGTTGCATTCGAATGACCTGATTAACGTAGCTAGCATCGATCGGATCCAATTGC"   # 55 bp
    seq = linear(rnd_left, insert, "TTTTCCCAGTCACGACGTTGTAAAACGACG")
    f = insert[:30]                   # forward footprint 0..29 of the insert
    r = rc(insert[15:])               # reverse footprint 15..54: overlaps the forward one
    [product] = pl.amplified_products([pl.Primer("F", f), pl.Primer("R", r)], seq, circular=False)
    assert product["size"] == 55


def test_the_library_pair_must_fit_the_product_size_and_flank_the_region():
    fwd, rev_site = "ACGTTGCATTCGAATGACCTG", "GGTTACCGATTACGCATTGCA"
    seq = ("T" * 100 + fwd + "C" * 700 + rc(rev_site) + "G" * 1158).upper()   # product 742 bp
    primers = [pl.Primer("F", fwd), pl.Primer("R", rev_site)]
    assert pl.product_over_region(primers, seq, False, 300, 600) is not None
    assert pl.product_over_region(primers, seq, False, 300, 600, 100, 500) is None   # too big
    # A region under the forward primer: the product carries the primer there, not the template.
    assert pl.product_over_region(primers, seq, False, 105, 400) is None


def test_a_sequencing_read_reaches_exactly_800_bases():
    fwd = "TTGCATTCGAATGACCTGAT"
    seq = linear(fwd, "A" * 1000, "")
    three = len(fwd) - 1
    assert pl.sequencing_primers([pl.Primer("S", fwd)], seq, False, three + 40, three + 800)
    assert not pl.sequencing_primers([pl.Primer("S", fwd)], seq, False, three + 40, three + 801)


def test_a_reversed_region_on_a_linear_sequence_is_refused(client):
    r = client.post("/api/primer-library/product", headers=headers(), json={
        "sequenceData": {"sequence": "ACGT" * 100, "circular": False}, "start": 300, "end": 10})
    assert r.status_code == 422


# ---- third review: the merge, odd files, repeats -----------------------------------------

def test_a_new_column_goes_past_cells_the_header_does_not_name(library_file):
    library_file.write_text("Name;Sequence\nM13F;GTAAAACGACGGCCAGT;box 3\nM13R;CAGGAAACAGCTATGAC;box 4\n",
                            encoding="utf-8")
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "NEW"})])
    text = library_file.read_text(encoding="utf-8-sig")
    assert "M13F;GTAAAACGACGGCCAGT;box 3;NEW" in text and "M13R;CAGGAAACAGCTATGAC;box 4\n" in text
    assert text.startswith("Name;Sequence;;notes")


def test_a_quote_never_closed_is_refused_not_swallowing_the_rows_below(library_file):
    data = b'name;sequence;notes\nM13F;GTAAAACGACGGCCAGT;"box 3\nM13R;CAGGAAACAGCTATGAC;x\n'
    library_file.write_bytes(data)
    with pytest.raises(pl.PrimerLibraryError, match="quote"):
        pl.read_library()
    with pytest.raises(pl.PrimerLibraryError):
        pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "NEW"})])
    assert library_file.read_bytes() == data


def test_title_rows_above_the_header_are_kept(library_file):
    body = "Lab primers 2024;;;\n;;;\nName;Tm;Sequence;Notes\nM13F;53;GTAAAACGACGGCCAGT;box 3\n"
    library_file.write_text(body, encoding="utf-8")
    [p] = pl.read_library()["primers"]
    assert (p["name"], p["notes"]) == ("M13F", "box 3")
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "box 4"})])
    assert library_file.read_text(encoding="utf-8-sig") == body.replace("box 3", "box 4")


@pytest.mark.parametrize("header,row", [
    ("Sequence name;Sequence;Notes", "M13F;GTAAAACGACGGCCAGT;box 3"),
    ("Sequence ID;Name;Sequence (5'-3');Notes", "17;M13F;GTAAAACGACGGCCAGT;box 3"),
    ("Name;Sequence length;Sequence;Notes", "M13F;17;GTAAAACGACGGCCAGT;box 3"),
])
def test_sequence_like_headers_that_are_not_the_sequence(header, row):
    [p], rejected = pl.parse_csv(f"{header}\n{row}\n")
    assert rejected == [] and (p.name, p.sequence, p.notes) == ("M13F", "GTAAAACGACGGCCAGT", "box 3")


def test_a_no_op_save_leaves_odd_rows_exactly_as_written(library_file):
    body = ('name;sequence;notes\r\n;GTAAAACGACGGCCAGT;nameless\r\n M13R ;CAGGAAACAGCTATGAC\r\n'
            'T7;TAATACGACTCACTATAGGG;"line1\nline2"\r\n')
    library_file.write_bytes(body.encode("utf-8"))
    edits = [pl.Change(p["sequence"], {k: p[k] for k in pl.EDITABLE}, {k: p[k] for k in pl.EDITABLE})
             for p in pl.read_library()["primers"]]
    pl.write_library([], edits)
    assert library_file.read_bytes().decode("utf-8-sig") == body


def test_a_headerless_file_is_not_padded_on_save(library_file):
    library_file.write_text("M13F,GTAAAACGACGGCCAGT\nM13R,CAGGAAACAGCTATGAC\n", encoding="utf-8")
    pl.write_library([pl.Primer("T7", "TAATACGACTCACTATAGGG")])
    assert library_file.read_text(encoding="utf-8-sig") == (
        "M13F,GTAAAACGACGGCCAGT\nM13R,CAGGAAACAGCTATGAC\nT7,TAATACGACTCACTATAGGG\n")


def test_a_utf8_file_with_one_windows_line_keeps_both(library_file):
    library_file.write_bytes("name;sequence;notes\nM13F;GTAAAACGACGGCCAGT;5 µM\n".encode()
                             + "M13R;CAGGAAACAGCTATGAC;4 °C\n".encode("cp1252"))
    notes = [p["notes"] for p in pl.read_library()["primers"]]
    assert notes == ["5 µM", "4 °C"]


def test_a_refused_file_does_not_become_the_library(client, tmp_path):
    good = settings.primer_library_path()
    bad = tmp_path / "workbook.csv"
    bad.write_bytes(b"PK\x03\x04xlsx")
    r = client.post("/api/primer-library/path", json={"path": str(bad)}, headers=headers())
    assert r.status_code == 422 and settings.primer_library_path() == good


@pytest.mark.parametrize("pasted", [
    "GTAAAACGACGGCCAGT\nCAGGAAACAGCTATGAC\n",                 # one column of sequences
    "GTAAAACGACGGCCAGT\tM13F\nCAGGAAACAGCTATGAC\tM13R\n",     # sequence, then name
])
def test_pastes_with_the_sequence_first(pasted):
    primers, rejected = pl.parse_csv(pasted)
    assert rejected == [] and [p.sequence for p in primers] == ["GTAAAACGACGGCCAGT", "CAGGAAACAGCTATGAC"]


def test_a_primer_in_a_short_repeat_is_not_offered_as_a_product():
    fwd = "ACGTTGCATTCGAATGACCTG"
    rep = "ACACACACACACACACACAC"
    seq = ("G" * 50 + fwd + "T" * 300 + rep + "ACAC" + "G" * 50).upper()
    rev = rc(rep)
    [site] = [b for b in pl.find_bindings(rev, seq, False) if b.strand == -1]
    assert site.positions > 1
    assert pl.amplified_products([pl.Primer("F", fwd), pl.Primer("R", rev)], seq, False) == []


# ---- fourth review: two tubes of one oligo, look-alike headers, nearby sites -------------

TUBES = "name;sequence;location\nM13F-tube2;GTAAAACGACGGCCAGT;box7\nM13F;GTAAAACGACGGCCAGT;box1\n"
LOADED_M13F = {"name": "M13F", "notes": "", "location": "box1", "added": ""}


def test_an_edit_lands_on_the_tube_it_was_made_on(library_file):
    # The dialog loaded the file with M13F first; a colleague then put tube 2 above it.
    library_file.write_text(TUBES, encoding="utf-8")
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", LOADED_M13F, {"location": "box1 A4"})])
    text = library_file.read_text(encoding="utf-8-sig")
    assert "M13F-tube2;GTAAAACGACGGCCAGT;box7" in text and "M13F;GTAAAACGACGGCCAGT;box1 A4" in text


def test_removing_a_tube_removes_that_tube_only(library_file):
    library_file.write_text(TUBES, encoding="utf-8")
    pl.write_library([], removed=[pl.Change("GTAAAACGACGGCCAGT", LOADED_M13F)])
    assert library_file.read_text(encoding="utf-8-sig") == (
        "name;sequence;location\nM13F-tube2;GTAAAACGACGGCCAGT;box7\n")
    # The colleague had already removed that tube: the other tube is left alone.
    result = pl.write_library([], removed=[pl.Change("GTAAAACGACGGCCAGT", LOADED_M13F)])
    assert result["gone"] == 1 and [p["name"] for p in result["primers"]] == ["M13F-tube2"]


def test_a_primer_removed_and_added_again_comes_back_even_if_deleted_meanwhile(library_file):
    library_file.write_text("name;sequence\nM13R;CAGGAAACAGCTATGAC\n", encoding="utf-8")
    result = pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"},
                                             {"name": "M13F", "notes": "again"}, readd=True)])
    assert [p["name"] for p in result["primers"]] == ["M13R", "M13F"] and result["gone"] == 0


def test_a_data_row_starting_with_seq_is_not_a_header():
    text = "M13F;GTAAAACGACGGCCAGT\nSeq1;CAGGAAACAGCTATGAC;primer\nT7;TAATACGACTCACTATAGGG\n"
    primers, rejected = pl.parse_csv(text)
    assert rejected == [] and [p.name for p in primers] == ["M13F", "Seq1", "T7"]


def test_gene_names_are_not_taken_for_sequences():
    primers, _ = pl.parse_csv("ACTB;CATGTACGTTGCTATCCAGGC;box1\nMYC;to order;\nHBB;?;\n")
    assert [(p.name, p.sequence) for p in primers] == [("ACTB", "CATGTACGTTGCTATCCAGGC")]


@pytest.mark.parametrize("row", ['A;GTAAAACGACGGCCAGT;"GFP" fusion', '"A" old;GTAAAACGACGGCCAGT;x',
                                 'A;GTAAAACGACGGCCAGT;"Box 3", shelf'])
def test_hand_typed_quotes_are_read_as_excel_reads_them(row):
    primers, rejected = pl.parse_csv(f"name;sequence;notes\n{row}\nB;CAGGAAACAGCTATGAC;y\n")
    assert rejected == [] and len(primers) == 2


def test_a_primer_with_a_second_site_nearby_is_not_offered_for_sequencing_or_a_product():
    fwd = "TTGCATTCGAATGACCTGAT"
    seq = linear("G" * 200 + fwd + "CCCCC" + fwd, "A" * 300, "T" * 200)
    assert pl.sequencing_primers([pl.Primer("S", fwd)], seq, False, 330, 360) == []
    rev = "GGTTACCGATTACGCATTGC"
    seq2 = seq + rc(rev) + "G" * 50
    primers = [pl.Primer("F", fwd), pl.Primer("R", rev)]
    assert len(pl.amplified_products(primers, seq2, False)) == 2      # two bands...
    assert pl.product_over_region(primers, seq2, False, 300, 400) is None   # ...so not offered


# ---- fifth review ------------------------------------------------------------------------

SAME_NAME = "name;sequence;location\nM13F;GTAAAACGACGGCCAGT;box1\nM13F;GTAAAACGACGGCCAGT;box2\n"


def test_with_two_tubes_of_one_name_a_change_never_falls_on_the_other(library_file):
    # Loaded with two tubes; the colleague then deleted the box1 tube.
    library_file.write_text("name;sequence;location\nM13F;GTAAAACGACGGCCAGT;box2\n", encoding="utf-8")
    was = {"name": "M13F", "notes": "", "location": "box1", "added": ""}
    r = pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", was, {"location": "box1 A4"}, twins=2)],
                         [pl.Change("GTAAAACGACGGCCAGT", was, twins=2)])
    assert r["gone"] == 2
    assert library_file.read_text(encoding="utf-8-sig").endswith("M13F;GTAAAACGACGGCCAGT;box2\n")


def test_an_edit_follows_a_single_tube_a_colleague_renamed(library_file):
    library_file.write_text("name;sequence;notes\nM13F-new;GTAAAACGACGGCCAGT;\n", encoding="utf-8")
    r = pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F"}, {"notes": "mine"})])
    assert r["changed"] == 1 and r["primers"][0]["name"] == "M13F-new"


def test_adding_back_a_primer_still_in_the_file_keeps_its_date_and_notes(library_file):
    library_file.write_text("name;sequence;notes;added\nM13F;GTAAAACGACGGCCAGT;theirs;2020-01-01\n",
                            encoding="utf-8")
    pl.write_library([], [pl.Change("GTAAAACGACGGCCAGT", {"name": "M13F", "notes": ""},
                                    {"name": "M13F", "notes": "", "added": "2026-10-21"}, readd=True)])
    assert "M13F;GTAAAACGACGGCCAGT;theirs;2020-01-01" in library_file.read_text(encoding="utf-8-sig")


@pytest.mark.parametrize("body", [
    'Name;Sequence;Notes\nA;ACGTACGTACGTACGTAA;"odd\nB;ACGTACGTAC;x\nC;/5Phos/ACGTACGTACG;y\nD;to order;z\n',
    'Sequence;Name\nACGTACGTACGTACGTAA;"P1\nGTAAAACGACGGCCAGT;P2\nCAGGAAACAGCTATGAC;P3\n',
])
def test_an_unclosed_quote_anywhere_is_refused(library_file, body):
    library_file.write_text(body, encoding="utf-8")
    with pytest.raises(pl.PrimerLibraryError, match="never closed"):
        pl.read_library()


def test_a_plain_header_is_kept_even_when_its_sequences_carry_modifications():
    text = "Name;Target site;Sequence;Notes\nqF1;GGCCAATTGGCCAATTGG;/5Phos/ACGTACGTACGTACGTAA;x\n"
    primers, rejected = pl.parse_csv(text)
    assert primers == [] and len(rejected) == 1        # read under its header, not by position


def test_a_pair_is_not_offered_when_one_primer_also_pairs_with_itself():
    f, r = "ACGTTGCATTCGAATGACCTG", "GGTTACCGATTACGCATTGCA"
    seq = ("T" * 200 + f + "C" * 680 + rc(r) + "A" * 1100 + rc(f) + "G" * 500).upper()
    primers = [pl.Primer("F", f), pl.Primer("R", r)]
    assert pl.product_over_region(primers, seq, False, 300, 600) is None   # F with F: 2nd band
