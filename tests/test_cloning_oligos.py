import pytest
from cloning_helpers import clean_dna, doc, rc, same_circular
from fastapi import FastAPI
from fastapi.testclient import TestClient

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.oligos import anneal, design

BACKBONE = clean_dna(2500, 401)
DROPOUT = clean_dna(60, 402)

# BbsI is GAAGAC(2/6): it cuts two bases past the site on the top strand and six on the
# bottom, leaving a four-base 5' overhang. Placing one site on each side of a dropout, facing
# it, leaves the backbone with the overhangs the gRNA protocol expects: GTTT at its 5' end and
# CACC at its 3' end (so the insert needs CACC on top and AAAC on the bottom strand).
VECTOR_SEQ = (BACKBONE + "CACC" + "TT" + "GTCTTC" + DROPOUT + "GAAGAC" + "TT" + "GTTT")
VECTOR = SequenceDoc(**doc("pGuide", VECTOR_SEQ, True))
N20 = "GTCGCCCTCGAACTTCACCT"
CORE = "G" + N20


def test_designed_oligos_follow_the_protocol_shape():
    result = design(CORE, VECTOR, ["BbsI"], name="gRNA1")
    top, bottom = result.oligos
    assert top.name == "gRNA1_top" and bottom.name == "gRNA1_bottom"
    assert top.sequence == "CACC" + CORE
    assert bottom.sequence == "AAAC" + rc(CORE)


def test_the_designed_product_is_the_vector_with_the_insert_in_place():
    [product] = design(CORE, VECTOR, ["BbsI"]).products
    expected = "GTTT" + BACKBONE + "CACC" + CORE
    assert product.doc["circular"] is True
    assert same_circular(product.doc["sequence"], expected)
    assert CORE in product.doc["sequence"]
    assert DROPOUT not in product.doc["sequence"]      # the filler is gone


def test_feeding_the_designed_oligos_back_gives_the_same_product():
    designed = design(CORE, VECTOR, ["BbsI"])
    top, bottom = designed.oligos
    annealed = anneal("top", top.sequence, "bottom", bottom.sequence, VECTOR, ["BbsI"])
    assert annealed.products[0].doc["sequence"] == designed.products[0].doc["sequence"]


def test_annealing_without_a_vector_describes_the_duplex():
    [product] = anneal("o1", "CACC" + CORE, "o2", "AAAC" + rc(CORE)).products
    assert "5' CACC / 5' AAAC" in product.details
    assert product.doc["circular"] is False


def test_oligos_that_do_not_anneal_are_rejected():
    with pytest.raises(CloningError, match="do not anneal"):
        anneal("o1", "ACGTACGTACGTACGTAC", "o2", "TTTTTTTTTTTTTTTTTT")


def test_overhangs_that_do_not_match_the_vector_name_what_is_needed():
    with pytest.raises(CloningError) as err:
        anneal("o1", "AATT" + CORE, "o2", "AATT" + rc(CORE), VECTOR, ["BbsI"])
    assert "5′ CACC overhang on the top strand" in str(err.value)
    assert "5′ AAAC overhang on the bottom strand" in str(err.value)


def test_phosphorylation_is_always_reminded():
    result = design(CORE, VECTOR, ["BbsI"])
    assert any("T4 PNK" in w for w in result.warnings)
    assert any("Anneal the two oligos" in w for w in result.warnings)


def test_a_blunt_cut_vector_is_refused():
    blunt = SequenceDoc(**doc("pBlunt", clean_dna(1500, 403) + "GATATC" + clean_dna(800, 404), True))
    with pytest.raises(CloningError, match="blunt ends"):
        design(CORE, blunt, ["EcoRV"])


def test_a_vector_cut_more_than_twice_is_refused():
    busy = SequenceDoc(**doc("pBusy", VECTOR_SEQ + "GAAGAC" + clean_dna(40, 405), True))
    with pytest.raises(CloningError, match="BbsI cuts the vector 3 times; the vector must be cut once or twice"):
        design(CORE, busy, ["BbsI"])


def test_http_route():
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    res = client.post("/api/cloning/anneal-oligos",
                      json={"mode": "design", "insert": CORE, "insert_name": "gRNA1",
                            "vector": doc("pGuide", VECTOR_SEQ, True), "enzymes": ["BbsI"]})
    assert res.status_code == 200
    assert [o["sequence"] for o in res.json()["oligos"]] == ["CACC" + CORE, "AAAC" + rc(CORE)]
