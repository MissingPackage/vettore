import pytest
from cloning_helpers import clean_dna, doc, feature, rc, same_circular
from fastapi import FastAPI
from fastapi.testclient import TestClient

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.overlap import InsertIn, Linearize, VectorIn, overlap_assembly

BB = clean_dna(2400, 41)
VECTOR_SEQ = BB[:1200] + "GAATTC" + BB[1200:]
VECTOR = SequenceDoc(**doc("pAAV-MCS", VECTOR_SEQ, True, features=[feature("AmpR", 50, 649)]))
GFP = clean_dna(720, 42)
INSERT = SequenceDoc(**doc("GFP", GFP, False, features=[feature("GFP CDS", 0, 719)]))
WPRE = clean_dna(400, 43)
INSERT2 = SequenceDoc(**doc("WPRE", WPRE, False))


def ecori_backbone() -> str:
    """The EcoRI-linearised vector at full length, as the kits repair it: the plasmid rotated
    to the top-strand cut (G^AATTC) with the 4 nt overhang present on both ends."""
    cut = VECTOR_SEQ.index("GAATTC") + 1
    return VECTOR_SEQ[cut:] + VECTOR_SEQ[:cut] + "AATT"


def test_gibson_one_insert_into_an_ecori_cut_vector():
    result = overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                              [InsertIn(doc=INSERT)])
    [product] = result.products
    expected = ecori_backbone() + GFP
    assert len(product.doc["sequence"]) == len(expected)
    assert same_circular(product.doc["sequence"], expected)
    assert product.doc["circular"] is True
    names = {f["name"] for f in product.doc["features"]}
    assert names == {"AmpR", "GFP CDS"}
    assert product.operation == "Gibson assembly"


def test_gibson_primers_carry_the_homology_arms():
    result = overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                              [InsertIn(doc=INSERT)], overlap=25)
    forward, reverse = result.oligos
    backbone = ecori_backbone()
    assert forward.sequence[:25] == backbone[-25:]
    assert forward.sequence[25:] == GFP[: forward.annealing_length]
    assert reverse.sequence[:25] == rc(backbone[:25])
    assert reverse.sequence[25:] == rc(GFP[-reverse.annealing_length :])


def test_two_inserts_keep_the_requested_order():
    result = overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                              [InsertIn(doc=INSERT), InsertIn(doc=INSERT2)])
    expected = ecori_backbone() + GFP + WPRE
    assert same_circular(result.products[0].doc["sequence"], expected)
    assert len(result.oligos) == 4


def test_in_fusion_uses_15_bp_arms_by_default():
    result = overlap_assembly("infusion", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                              [InsertIn(doc=INSERT)])
    assert result.oligos[0].sequence[:15] == ecori_backbone()[-15:]
    assert "15 bp homology arms" in result.products[0].details
    assert result.products[0].operation == "In-Fusion cloning"


def test_vector_linearised_by_pcr_gets_its_own_primers():
    result = overlap_assembly("hifi", VectorIn(doc=VECTOR, linearize=Linearize(pcr={"start": 0, "end": 1999})),
                              [InsertIn(doc=INSERT)])
    expected = VECTOR_SEQ[:2000] + GFP
    assert same_circular(result.products[0].doc["sequence"], expected)
    assert [o.name for o in result.oligos][:2] == ["pAAV-MCS_F", "pAAV-MCS_R"]
    assert len(result.oligos) == 4


def test_insert_region_of_a_circular_donor_is_required():
    donor = SequenceDoc(**doc("pDONOR", clean_dna(1500, 44), True))
    with pytest.raises(CloningError, match="choose the region to amplify"):
        overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                         [InsertIn(doc=donor)])


def test_insert_region_of_a_circular_donor_is_amplified():
    core = clean_dna(1500, 45)
    donor = SequenceDoc(**doc("pDONOR", core, True))
    result = overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                              [InsertIn(doc=donor, start=100, end=699)])
    assert same_circular(result.products[0].doc["sequence"], ecori_backbone() + core[100:700])


def test_too_short_overlap_is_rejected():
    with pytest.raises(CloningError, match="at least 10 bp"):
        overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["EcoRI"])),
                         [InsertIn(doc=INSERT)], overlap=5)


def test_enzyme_that_does_not_cut_is_reported():
    with pytest.raises(CloningError, match="does not cut"):
        overlap_assembly("gibson", VectorIn(doc=VECTOR, linearize=Linearize(enzymes=["NotI"])),
                         [InsertIn(doc=INSERT)])


def test_http_route():
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    body = {"method": "gibson",
            "vector": {"doc": doc("pAAV-MCS", VECTOR_SEQ, True), "linearize": {"enzymes": ["EcoRI"]}},
            "inserts": [{"doc": doc("GFP", GFP, False)}]}
    res = client.post("/api/cloning/overlap-assembly", json=body)
    assert res.status_code == 200
    assert res.json()["products"][0]["doc"]["circular"] is True
