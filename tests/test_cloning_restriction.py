import pytest
from Bio import Restriction
from Bio.Seq import Seq
from cloning_helpers import clean_dna, doc, feature, rc, same_circular
from fastapi import FastAPI
from fastapi.testclient import TestClient

from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.restriction import delete_fragment, digest, insert_fragment, router

A = clean_dna(2000, 1)
B = clean_dna(1000, 2)
MCS = "GAATTC" "GGTACC" "GGATCC" "AGATCT" "TCTAGA" "CTGCAG" "AAGCTT"
VECTOR_SEQ = A + MCS + B
CORE = clean_dna(700, 3)
INSERT_SEQ = "CCC" + "GAATTC" + CORE + "GGATCC" + "CCC"


def vector_doc():
    return SequenceDoc(**doc("pVEC", VECTOR_SEQ, True,
                             features=[feature("AmpR", 100, 899), feature("ori", 2200, 2799)]))


def insert_doc():
    return SequenceDoc(**doc("GFP PCR", INSERT_SEQ, False,
                             features=[feature("GFP", 9 + 50, 9 + 649)]))


def test_directional_ecori_bamhi_insertion_gives_exact_product():
    result = insert_fragment(vector_doc(), ["EcoRI", "BamHI"], insert_doc(), ["EcoRI", "BamHI"])
    [product] = result.products
    expected = A + "GAATTC" + CORE + "GGATCC" + "AGATCT" "TCTAGA" "CTGCAG" "AAGCTT" + B
    seq = product.doc["sequence"]
    assert product.doc["circular"] is True
    assert seq.upper() == expected  # vector numbering is kept: same origin as the vector
    assert same_circular(seq, expected)
    assert product.doc["name"] == "pVEC + GFP PCR"
    assert product.operation == "Restriction cloning"
    assert product.parents == ["pVEC", "GFP PCR"]
    assert "EcoRI–BamHI fragment of GFP PCR" in product.details
    assert not any("both orientations" in w for w in result.warnings)


def test_features_of_both_parents_are_carried_with_shifted_coordinates():
    [product] = insert_fragment(vector_doc(), ["EcoRI", "BamHI"], insert_doc(),
                                ["EcoRI", "BamHI"]).products
    feats = {f["name"]: (f["start"], f["end"]) for f in product.doc["features"]}
    assert feats["AmpR"] == (100, 899)
    gfp_start = len(A) + 6 + 50
    assert feats["GFP"] == (gfp_start, gfp_start + 599)
    shift = 700 + 6 - (6 + 6)  # insert (GAATTC+CORE+G..) replaces GAATTCGGTACCG
    assert feats["ori"] == (2200 + shift, 2799 + shift)


def test_wrong_orientation_for_directional_cloning_is_an_error_with_hint():
    with pytest.raises(CloningError) as err:
        insert_fragment(vector_doc(), ["EcoRI", "BamHI"], insert_doc(), ["EcoRI", "BamHI"],
                        orientation="reverse")
    assert "compatible in the other orientation" in str(err.value)


def test_single_enzyme_cloning_is_non_directional_and_both_orientations_work():
    ins = SequenceDoc(**doc("frag", "CC" + "GAATTC" + CORE + "GAATTC" + "CC", False))
    fwd = insert_fragment(vector_doc(), ["EcoRI"], ins, ["EcoRI"])
    rev = insert_fragment(vector_doc(), ["EcoRI"], ins, ["EcoRI"], orientation="reverse")
    rest = "GGTACCGGATCCAGATCTTCTAGACTGCAGAAGCTT" + B
    assert fwd.products[0].doc["sequence"].upper() == A + "GAATTC" + CORE + "GAATTC" + rest
    assert rev.products[0].doc["sequence"].upper() == A + "GAATTC" + rc(CORE) + "GAATTC" + rest
    assert any("both orientations" in w for w in fwd.warnings)
    assert any("re-ligate without insert" in w for w in fwd.warnings)


def test_incompatible_ends_are_reported():
    ins = SequenceDoc(**doc("frag", "CC" + "GAATTC" + CORE + "CTGCAG" + "CC", False))
    with pytest.raises(CloningError) as err:
        insert_fragment(vector_doc(), ["EcoRI", "BamHI"], ins, ["EcoRI", "PstI"])
    assert "do not match" in str(err.value)


def test_enzyme_that_does_not_cut_is_an_error():
    with pytest.raises(CloningError, match="NotI does not cut the vector"):
        insert_fragment(vector_doc(), ["NotI"], insert_doc(), ["EcoRI"])


def test_unknown_enzyme_is_an_error():
    with pytest.raises(CloningError, match="Unknown or unsupported enzyme"):
        digest(vector_doc(), ["NotAnEnzyme"])


def test_double_cutter_in_vector_is_warned():
    vec = SequenceDoc(**doc("pDUP", A + MCS + B[:500] + "GGATCC" + B[500:], True))
    result = insert_fragment(vec, ["EcoRI", "BamHI"], insert_doc(), ["EcoRI", "BamHI"])
    assert any("BamHI cuts the vector 2 times" in w for w in result.warnings)


def test_insert_excised_from_a_donor_plasmid_uses_the_smaller_fragment():
    donor_backbone = clean_dna(2500, 4)
    donor = SequenceDoc(**doc("pDONOR", donor_backbone + "GAATTC" + CORE + "GGATCC", True))
    [product] = insert_fragment(vector_doc(), ["EcoRI", "BamHI"], donor, ["EcoRI", "BamHI"]).products
    assert product.doc["sequence"].upper() == (A + "GAATTC" + CORE + "GGATCC"
                                               + "AGATCTTCTAGACTGCAGAAGCTT" + B)


def test_digest_fragments_match_biopython_and_report_origin_spanning_fragment():
    frags = digest(vector_doc(), ["EcoRI", "HindIII"])
    seq = Seq(VECTOR_SEQ)
    eco = Restriction.EcoRI.search(seq, linear=False)[0] - 1   # 0-based index after the cut
    hind = Restriction.HindIII.search(seq, linear=False)[0] - 1
    lengths = sorted(f["length"] for f in frags)
    assert lengths == sorted([hind - eco, len(VECTOR_SEQ) - (hind - eco)])
    wrapping = next(f for f in frags if f["start"] > f["end"])
    assert (wrapping["start"], wrapping["end"]) == (hind, eco - 1)
    assert wrapping["left_enzyme"] == "HindIII" and wrapping["right_enzyme"] == "EcoRI"
    assert wrapping["left_end"] == "5' overhang AGCT"


def test_single_enzyme_fragment_lengths_match_biopython_catalyse():
    seq = Seq(INSERT_SEQ)
    expected = [len(f) for f in Restriction.EcoRI.catalyse(seq, linear=True)]
    assert [f["length"] for f in digest(insert_doc(), ["EcoRI"])] == expected


def test_linear_digest_reports_molecule_ends():
    frags = digest(insert_doc(), ["EcoRI", "BamHI"])
    # CCC G|AATTC <700> G|GATCC CCC: top-strand cuts at 4 and 710 of 718 bp.
    assert [f["length"] for f in frags] == [4, 706, 8]
    assert frags[0]["left_end"] == "end of molecule"
    assert frags[1]["label"] == "EcoRI–BamHI (706 bp)"


def test_delete_fragment_between_compatible_bamhi_and_bglii_sites():
    [product] = delete_fragment(vector_doc(), ["BamHI", "BglII"]).products
    i = VECTOR_SEQ.index("GGATCCAGATCT")
    expected = VECTOR_SEQ[: i + 1] + VECTOR_SEQ[i + 7:]  # drop GATCCA, joining G|GATCT
    assert product.doc["sequence"].upper() == expected
    assert "GGATCT" in expected[i - 2: i + 8]


def test_delete_fragment_with_incompatible_ends_is_an_error():
    with pytest.raises(CloningError, match="blunting is not supported"):
        delete_fragment(vector_doc(), ["KpnI", "XbaI"])


def test_http_routes():
    app = FastAPI()
    app.include_router(router, prefix="/api/cloning")
    client = TestClient(app)
    body = {"vector": doc("pVEC", VECTOR_SEQ, True), "vector_enzymes": ["EcoRI", "BamHI"],
            "insert": doc("GFP PCR", INSERT_SEQ, False), "insert_enzymes": ["ecori", "bamhi"]}
    res = client.post("/api/cloning/restriction-insert", json=body)
    assert res.status_code == 200
    assert len(res.json()["products"][0]["doc"]["sequence"]) == len(A) + 6 + 700 + 6 + 24 + len(B)
    res = client.post("/api/cloning/digest", json={"doc": doc("x", VECTOR_SEQ, True),
                                                   "enzymes": ["EcoRI"]})
    assert res.json()["fragments"][0]["length"] == len(VECTOR_SEQ)
