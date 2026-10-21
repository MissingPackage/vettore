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
    again, _ = pl.parse_csv(pl.to_csv(primers))
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
    saved = client.post("/api/primer-library/save", json={"primers": rows}, headers=headers())
    assert saved.status_code == 200 and saved.json()["count"] == 1
    got = client.get("/api/primer-library", headers=headers())
    assert got.status_code == 200
    body = got.json()
    assert body["exists"] and body["primers"][0]["name"] == "M13F"
    assert body["primers"][0]["tm"] is not None


def test_bindings_route_finds_a_stored_primer_on_a_sequence(client):
    seq = ("AAAACCCC" + PRIMER + "GGGGTTTT").upper()
    client.post("/api/primer-library/save",
                json={"primers": [{"name": "P1", "sequence": PRIMER}]}, headers=headers())
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
    result = pl.write_library([pl.Primer("M13F", "GTAAAACGACGGCCAGT", "box 3"), new])
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
    # 21 overlapping 20-mers fit the 40-A stretch; only non-overlapping sites are kept (two).
    assert len(pl.find_bindings(oligo_dt, seq, circular=False)) == 2
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
