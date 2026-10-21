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
    assert rc(seq[b.end:b.start + 1]) == PRIMER


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


def test_short_primer_requires_the_whole_oligo_to_match():
    short = "ACGTACGTAC"      # 10 nt, shorter than the 15-nt anchor
    seq = linear("TTTT", short, "GGGG")
    [b] = pl.find_bindings(short, seq, circular=False)
    assert b.footprint == len(short) and b.mismatches == 0


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
    assert pl.SEQ_PRIMER_MIN <= reads[0]["reads_from"] <= pl.SEQ_PRIMER_MAX


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
    assert len(rejected) == 1 and rejected[0]["name"] == "Bad"


def test_a_repeated_sequence_is_reported_not_duplicated():
    text = "name,sequence\nA,ACGTACGTAC\nB,ACGTACGTAC\n"
    primers, rejected = pl.parse_csv(text)
    assert len(primers) == 1 and len(rejected) == 1 and "repeats" in rejected[0]["reason"]


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
