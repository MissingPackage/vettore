"""PCR primer design with primer3 (``vettore.primer_design``).

Every returned pair is checked against the template by hand: the forward primer is the top
strand at its position, the reverse primer the reverse complement of the bottom strand, both
lie outside the target, and the product size is the distance between their 5' ends.
"""

import math
import sys
from itertools import pairwise
from pathlib import Path

import pytest
from Bio.SeqUtils import MeltingTemp as mt
from cloning_helpers import clean_dna, doc, rc
from fastapi.testclient import TestClient

from vettore import fileio, primer_design, tm
from vettore.primer_design import DesignRequest, PrimerDesignError, design
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
TOKEN = "t0ken"
TEMPLATE = clean_dna(1500, 411)


def site(template: str, start: int, end: int) -> str:
    """Top-strand bases start..end (0-based inclusive), across the origin if end < start."""
    return template[start:end + 1] if start <= end else template[start:] + template[:end + 1]


def check_pair(pair: dict, template: str, target_start: int, target_end: int, circular: bool,
               body: DesignRequest, product_range: list[int]) -> None:
    length = len(template)
    fwd, rev = pair["forward"], pair["reverse"]
    assert fwd["sequence"] == site(template, fwd["start"], fwd["end"])
    assert rev["sequence"] == rc(site(template, rev["start"], rev["end"]))
    # Product: from the forward primer's 5' end to the reverse primer's 5' end.
    assert pair["product_size"] == (rev["end"] - fwd["start"]) % length + 1
    assert product_range[0] <= pair["product_size"] <= product_range[1]
    # Both primers lie outside the target, which the product contains.
    def offset(p: int) -> int:          # position along the product
        return (p - fwd["start"]) % length

    assert offset(fwd["end"]) < offset(target_start)
    assert offset(target_end) < offset(rev["start"])
    for primer in (fwd, rev):
        assert body.size_min <= primer["length"] <= body.size_max
        assert body.gc_min <= primer["gc"] <= body.gc_max
        # Displayed Tm is Vettore's Standard preset, and primer3 optimised the same number.
        assert primer["tm"] == round(tm.primer_tm(primer["sequence"]), 1)
        assert abs(primer["tm"] - primer["tm_primer3"]) <= 0.5
        assert body.tm_min - 0.5 <= primer["tm"] <= body.tm_max + 0.5
    assert abs(fwd["tm_primer3"] - rev["tm_primer3"]) <= body.max_tm_difference + 0.1
    if not circular:
        assert fwd["start"] <= fwd["end"] < rev["start"] <= rev["end"]


def test_pairs_flank_the_target_and_match_the_template():
    body = DesignRequest(template=doc("t", TEMPLATE, False), start=600, end=799)
    out = design(body)
    assert len(out["pairs"]) == primer_design.NUM_RETURN
    # Default product range: the smallest product that holds the region and two primers,
    # up to 300 bp more.
    assert out["product_range"] == [200 + 2 * 18, 200 + 2 * 18 + 300]
    for pair in out["pairs"]:
        check_pair(pair, TEMPLATE, 600, 799, False, body, out["product_range"])
    assert [p["rank"] for p in out["pairs"]] == [1, 2, 3, 4, 5]
    assert "SantaLucia" in out["tm_method"] and "250 nM" in out["tm_method"]


def test_user_constraints_are_respected():
    body = DesignRequest(template=doc("t", TEMPLATE, False), start=500, end=699,
                         product_min=400, product_max=450, size_min=20, size_opt=22,
                         size_max=24, tm_min=58, tm_opt=60, tm_max=62, gc_min=45, gc_max=60,
                         max_tm_difference=2)
    out = design(body)
    assert out["pairs"]
    for pair in out["pairs"]:
        check_pair(pair, TEMPLATE, 500, 699, False, body, [400, 450])


def plasmid() -> dict:
    return fileio.read_documents("MT612431.1.gb", (DATA / "MT612431.1.gb").read_bytes())[0]


def test_a_target_across_the_origin_of_a_circular_plasmid():
    data = plasmid()
    sequence = data["sequence"].upper()
    body = DesignRequest(template=data, start=7100, end=100)
    out = design(body)
    assert out["target"]["length"] == len(sequence) - 7100 + 101
    for pair in out["pairs"]:
        check_pair(pair, sequence, 7100, 100, True, body, out["product_range"])
    # The product crosses the origin: it starts before 7100 and ends after 100.
    assert all(p["product_start"] > p["product_end"] for p in out["pairs"])


def test_pairs_without_warnings_are_listed_first():
    data = plasmid()
    cds = next(f for f in data["features"] if f["name"] == "Gephyrin.FingR-GFP-CCR5TC")
    out = design(DesignRequest(template=data, start=cds["start"], end=cds["end"]))
    flags = [bool(p["warnings"]) for p in out["pairs"]]
    assert flags == sorted(flags)                       # clean pairs first
    clean = [p["penalty"] for p in out["pairs"] if not p["warnings"]]
    assert clean == sorted(clean)


def test_structure_warnings_name_the_primer_and_the_dg():
    # The 20-mer is its own reverse complement, so it pairs with itself along its length.
    primer = "GCGAATTCGCGCGAATTCGC"
    assert primer == rc(primer)
    other = "ATGGTGAGCAAGGGCGAGGAG"
    structures = primer_design._structures(primer_design.primer3_bindings(), primer, other)
    assert structures["forward_self_dimer"] < primer_design.DG_WARNING
    warnings = primer_design.pair_warnings(structures, primer, other, TEMPLATE, False, "t")
    assert any(w.startswith("The forward primer can pair with itself (self-dimer) (ΔG -")
               and "kcal/mol at 37 °C" in w for w in warnings)
    assert not any("reverse primer can pair with itself" in w for w in warnings)


def nearest_neighbour_dg(seq: str, sodium_molar: float, celsius: float) -> float:
    """ΔG (kcal/mol) of a perfect duplex worked out from the SantaLucia & Hicks 2004
    nearest-neighbour table (Biopython's ``DNA_NN4``) with the SantaLucia 1998 entropy salt
    correction: initiation + terminal A·T penalties + the sum over the 19 stacks."""
    table = mt.DNA_NN4
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    dh, ds = table["init"]
    for end in (seq[0], seq[-1]):
        h, s = table["init_A/T" if end in "AT" else "init_G/C"]
        dh, ds = dh + h, ds + s
    for a, b in pairwise(seq):
        key = f"{a}{b}/{complement[a]}{complement[b]}"
        h, s = table[key] if key in table else table[f"{complement[b]}{complement[a]}/{b}{a}"]
        dh, ds = dh + h, ds + s
    ds += 0.368 * (len(seq) - 1) * math.log(sodium_molar)
    return dh - (celsius + 273.15) * ds / 1000.0


@pytest.mark.parametrize("primer", ["ATGGTGAGCAAGGGCGAGGAG", "GCTAGTTATTGCTCAGCGG"])
def test_dimer_dg_agrees_with_the_nearest_neighbour_sum(primer):
    bindings = primer_design.primer3_bindings()
    structures = primer_design._structures(bindings, primer, rc(primer))
    expected = nearest_neighbour_dg(primer, 0.05, 37.0)
    assert structures["hetero_dimer"] == pytest.approx(expected, abs=0.1)


def test_primer3_tm_matches_the_standard_preset():
    bindings = primer_design.primer3_bindings()
    for primer in ("AGCTTGCATGCCTGCAGGTCG", "GTAAAACGACGGCCAGT", "TAATACGACTCACTATAGGG",
                   "TTACTTGTACAGCTCGTCCATGCC"):
        value = bindings.calc_tm(primer, mv_conc=50, dv_conc=0, dntp_conc=0,
                                 dna_conc=primer_design.PRIMER3_DNA_CONC)
        assert value == pytest.approx(tm.primer_tm(primer), abs=0.2)


def test_a_repeated_three_prime_end_is_reported():
    repeat = TEMPLATE[560:600]           # 40 bases where primer3 finds forward primers
    template = repeat + TEMPLATE[600:1300] + repeat + TEMPLATE[1300:]
    # The target starts right after the first copy, so every forward primer lies inside the
    # repeat and its 3' end also matches the second copy.
    out = design(DesignRequest(template=doc("t", template, False), start=40, end=239))
    assert out["pairs"]
    for pair in out["pairs"]:
        assert pair["forward"]["end"] < 40
        assert any("forward primer also match 1 other place" in w for w in pair["warnings"])


@pytest.mark.parametrize("kwargs, message", [
    ({"start": 0, "end": 50}, "too close to the start or the end"),
    ({"start": 900, "end": 800}, "start is after its end"),
    ({"start": 600, "end": 799, "product_max": 230}, "at least 236 bp"),
    ({"start": 600, "end": 799, "product_min": 500, "product_max": 400}, "minimum is larger"),
    ({"start": 600, "end": 799, "tm_min": 75, "tm_opt": 78, "tm_max": 80}, "Tm too low"),
    ({"start": 600, "end": 799, "size_min": 25, "size_opt": 20, "size_max": 30}, "length"),
    ({"start": 600, "end": 5000}, "between 1 and 1,500"),
])
def test_impossible_requests_explain_what_to_change(kwargs, message):
    with pytest.raises(PrimerDesignError) as error:
        design(DesignRequest(template=doc("t", TEMPLATE, False), **kwargs))
    assert message in str(error.value)


def test_ambiguous_bases_are_avoided():
    template = TEMPLATE[:560] + "N" * 30 + TEMPLATE[590:]
    out = design(DesignRequest(template=doc("t", template, False), start=600, end=799))
    for pair in out["pairs"]:
        assert "N" not in pair["forward"]["sequence"] + pair["reverse"]["sequence"]


# ---- primer3 availability (packaged builds) -----------------------------------------------

def test_status_reports_the_parameter_files():
    out = primer_design.status()
    assert out["ready"] is True and out["missing"] == []
    assert out["parameters"].endswith("primer3_config")


def test_missing_parameter_files_give_an_actionable_message(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "primer3.bindings", None)    # import now fails
    with pytest.raises(PrimerDesignError) as error:
        design(DesignRequest(template=doc("t", TEMPLATE, False), start=600, end=799))
    assert "Unzip a fresh copy of Vettore" in str(error.value)
    monkeypatch.setattr(primer_design, "parameters_dir", lambda: tmp_path)
    out = primer_design.status()
    assert out["ready"] is False and "stack.ds" in out["missing"]


def test_a_non_ascii_install_folder_is_named_on_windows(monkeypatch):
    monkeypatch.setattr(primer_design.os, "name", "nt")
    message = primer_design._unavailable_message(Path("C:/Users/José/Vettore/primer3_config"))
    assert "accented" in message and "C:\\Vettore" in message


def test_http_routes():
    client = TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")
    headers = {TOKEN_HEADER: TOKEN}
    assert client.get("/api/primer-design/status", headers=headers).json()["ready"] is True
    ok = client.post("/api/primer-design/design", headers=headers,
                     json={"template": doc("t", TEMPLATE, False), "start": 600, "end": 799})
    assert ok.status_code == 200 and len(ok.json()["pairs"]) == 5
    bad = client.post("/api/primer-design/design", headers=headers,
                      json={"template": doc("t", TEMPLATE, False), "start": 0, "end": 10})
    assert bad.status_code == 422 and "Primers need" in bad.json()["detail"]
