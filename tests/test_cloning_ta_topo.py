"""TA, GC and TOPO cloning. Expected products are written out by hand from the ends each method
relies on (sources in vettore/cloning/gateway_ta.py): a T-vector opened between X and Y gives
X + T + insert + A + Y, which is the familiar ...GATT[insert]AATC... junction of a vector
opened at EcoRV (GAT^ATC) and T-tailed."""

import pytest
from cloning_helpers import clean_dna, doc, feature, rc
from fastapi import FastAPI
from fastapi.testclient import TestClient

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.gateway_ta import ta_topo
from vettore.cloning.pcr import PrimerIn

X = clean_dna(1200, 201) + "GAT"
Y = "ATC" + clean_dna(800, 202)
VECTOR = SequenceDoc(**doc("pT", X + Y, True, features=[feature("lacZa", len(X) - 50, len(X) + 49),
                                                         feature("AmpR", 100, 899)]))
INSERT_SEQ = clean_dna(400, 203)
INSERT = SequenceDoc(**doc("amplicon", INSERT_SEQ, False, features=[feature("ORF", 30, 329)]))
OPEN = len(X)


def seqs(result):
    return [p.doc["sequence"].upper() for p in result.products]


def test_ta_cloning_gives_both_orientations_with_t_and_a_at_the_junctions():
    result = ta_topo("TA", VECTOR, OPEN, insert=INSERT)
    assert seqs(result) == [X + "T" + INSERT_SEQ + "A" + Y, X + "T" + rc(INSERT_SEQ) + "A" + Y]
    forward, backward = result.products
    assert forward.label.startswith("Insert forward") and backward.label.startswith("Insert reversed")
    assert forward.operation == "TA cloning"
    assert forward.parents == ["pT", "amplicon"]
    assert forward.doc["name"] == "pT + amplicon"
    assert backward.doc["name"] == "pT + amplicon, product 2"
    assert any("either orientation" in w for w in result.warnings)
    assert any("Taq" in w for w in result.warnings)


def test_ta_product_keeps_vector_numbering_and_features():
    [product, _] = ta_topo("TA", VECTOR, OPEN, insert=INSERT).products
    feats = {f["name"]: (f["start"], f["end"]) for f in product.doc["features"]}
    assert feats["AmpR"] == (100, 899)
    assert feats["ORF"] == (OPEN + 1 + 30, OPEN + 1 + 329)   # after the vector's T
    assert "lacZa" not in feats       # disrupted by the insert


def test_gc_cloning_uses_c_and_g_overhangs():
    result = ta_topo("GC", VECTOR, OPEN, insert=INSERT)
    assert seqs(result) == [X + "C" + INSERT_SEQ + "G" + Y, X + "C" + rc(INSERT_SEQ) + "G" + Y]
    assert result.products[0].operation == "GC cloning"


def test_topo_ta_is_ta_cloning_into_a_topoisomerase_vector():
    result = ta_topo("TOPO-TA", VECTOR, OPEN, insert=INSERT)
    assert seqs(result) == seqs(ta_topo("TA", VECTOR, OPEN, insert=INSERT))
    assert result.products[0].operation == "TOPO TA cloning"


def test_blunt_topo_both_orientations_without_extra_bases():
    result = ta_topo("blunt-TOPO", VECTOR, OPEN, insert=INSERT)
    assert seqs(result) == [X + INSERT_SEQ + Y, X + rc(INSERT_SEQ) + Y]
    assert any("proofreading" in w for w in result.warnings)


def test_directional_topo_takes_the_cacc_end_one_way_only():
    cacc = SequenceDoc(**doc("ORF PCR", "CACC" + INSERT_SEQ, False))
    result = ta_topo("directional-TOPO", VECTOR, OPEN, insert=cacc)
    # The vector's GTGG overhang pairs with CACC, which is written once in the clone.
    assert seqs(result) == [X + "CACC" + INSERT_SEQ + Y]
    assert not any("either orientation" in w for w in result.warnings)
    assert result.products[0].operation == "Directional TOPO cloning"


def test_directional_topo_without_cacc_is_an_error_that_says_what_to_do():
    with pytest.raises(CloningError) as err:
        ta_topo("directional-TOPO", VECTOR, OPEN, insert=INSERT)
    assert "does not start with CACC" in str(err.value)
    assert "Add CACC to the 5′ end of the forward primer" in str(err.value)


def test_directional_topo_with_cacc_on_the_reverse_primer_goes_in_reversed():
    backwards = SequenceDoc(**doc("ORF PCR", INSERT_SEQ + rc("CACC"), False))
    result = ta_topo("directional-TOPO", VECTOR, OPEN, insert=backwards)
    assert seqs(result) == [X + "CACC" + rc(INSERT_SEQ) + Y]
    assert result.products[0].label.startswith("Insert reversed")
    assert any("reverse primer" in w for w in result.warnings)


def test_directional_topo_with_cacc_on_both_primers_warns_and_gives_both():
    both = SequenceDoc(**doc("ORF PCR", "CACC" + INSERT_SEQ + rc("CACC"), False))
    result = ta_topo("directional-TOPO", VECTOR, OPEN, insert=both)
    assert len(result.products) == 2
    assert any("Both ends of the PCR product start with CACC" in w for w in result.warnings)


def test_pcr_from_template_and_primers_with_the_cacc_tail():
    template = SequenceDoc(**doc("cDNA", clean_dna(100, 204) + INSERT_SEQ + clean_dna(100, 205),
                                 True))
    forward = PrimerIn(name="ORF_F", sequence="CACC" + INSERT_SEQ[:22])
    reverse = PrimerIn(name="ORF_R", sequence=rc(INSERT_SEQ[-22:]))
    result = ta_topo("directional-TOPO", VECTOR, OPEN, template=template, forward=forward,
                     reverse=reverse)
    assert seqs(result) == [X + "CACC" + INSERT_SEQ + Y]
    assert [o.name for o in result.oligos] == ["ORF_F", "ORF_R"]
    assert result.products[0].parents == ["pT", "cDNA"]


def test_directional_forward_primer_without_cacc_names_the_primer_to_order():
    template = SequenceDoc(**doc("cDNA", INSERT_SEQ, False))
    forward = PrimerIn(name="F", sequence=INSERT_SEQ[:22])
    reverse = PrimerIn(name="R", sequence=rc(INSERT_SEQ[-22:]))
    with pytest.raises(CloningError) as err:
        ta_topo("directional-TOPO", VECTOR, OPEN, template=template, forward=forward,
                reverse=reverse)
    assert f"Order the forward primer as CACC{INSERT_SEQ[:22]}" in str(err.value)


def test_linearised_vector_document_is_used_as_supplied():
    linear = SequenceDoc(**doc("pT linear", Y + X, False))
    result = ta_topo("TA", linear, insert=INSERT)
    assert seqs(result)[0] == Y + X + "T" + INSERT_SEQ + "A"


def test_circular_vector_needs_the_opening_position():
    with pytest.raises(CloningError, match="enter the position where the vector is opened"):
        ta_topo("TA", VECTOR, None, insert=INSERT)


def test_circular_insert_is_refused():
    with pytest.raises(CloningError, match="choose the linear PCR product"):
        ta_topo("TA", VECTOR, OPEN, insert=VECTOR)


def test_http_route():
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    res = client.post("/api/cloning/ta-topo", json={
        "kind": "TA", "vector": doc("pT", X + Y, True), "position": OPEN,
        "insert": doc("amplicon", INSERT_SEQ, False)})
    assert res.status_code == 200
    assert len(res.json()["products"]) == 2
