"""Multi-fragment and linear ligation. Expected sequences are written out by hand from the
enzymes' cut positions (EcoRI G^AATTC, BamHI G^GATCC, XhoI C^TCGAG, EcoRV GAT^ATC, all from
REBASE) and cross-checked against pydna's own ligation_assembly."""

import pytest
from Bio.Restriction import BamHI, EcoRI, XhoI
from cloning_helpers import clean_dna, doc, feature, rc, same_circular
from fastapi import FastAPI
from fastapi.testclient import TestClient
from pydna.assembly2 import ligation_assembly
from pydna.dseqrecord import Dseqrecord

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.ligation import FragmentIn, insert_multiple, ligate_fragments

A = clean_dna(1500, 101)
B = clean_dna(500, 102)
C1 = clean_dna(300, 103)
C2 = clean_dna(200, 104)
STUFFER = "AAAACCCC"
VECTOR_SEQ = A + "GAATTC" + STUFFER + "CTCGAG" + B
INS1_SEQ = "CC" + "GAATTC" + C1 + "GGATCC" + "CC"
INS2_SEQ = "CC" + "GGATCC" + C2 + "CTCGAG" + "CC"


def vector():
    return SequenceDoc(**doc("pVEC", VECTOR_SEQ, True, features=[feature("AmpR", 100, 899)]))


def ins1():
    return SequenceDoc(**doc("promoter", INS1_SEQ, False, features=[feature("P", 8 + 10, 8 + 209)]))


def ins2():
    return SequenceDoc(**doc("CDS", INS2_SEQ, False))


def frag(d, *enzymes, **kw):
    return FragmentIn(doc=d, enzymes=list(enzymes), **kw)


# ---- Insert Multiple Fragments -------------------------------------------------------------

def test_vector_and_two_directional_inserts_give_one_exact_product():
    result = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                             [frag(ins1(), "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    [product] = result.products
    # G^AATTC ... G^GATCC ... C^TCGAG: the stuffer leaves, the two inserts take its place.
    expected = A + "GAATTC" + C1 + "GGATCC" + C2 + "CTCGAG" + B
    assert product.doc["sequence"].upper() == expected   # vector numbering kept
    assert product.doc["circular"] is True
    assert product.doc["name"] == "pVEC + promoter + CDS"
    assert product.parents == ["pVEC", "promoter", "CDS"]
    assert product.operation == "Insert multiple fragments"
    assert product.label.startswith("pVEC → promoter → CDS")
    assert "Order: pVEC → promoter → CDS" in product.details
    assert not any("different products" in w for w in result.warnings)
    assert not any("close without an insert" in w for w in result.warnings)


def test_features_of_vector_and_inserts_are_carried():
    [product] = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                                [frag(ins1(), "EcoRI", "BamHI"),
                                 frag(ins2(), "BamHI", "XhoI")]).products
    feats = {f["name"]: (f["start"], f["end"]) for f in product.doc["features"]}
    assert feats["AmpR"] == (100, 899)
    start = len(A) + 6 + 10          # the insert's own coordinates, shifted behind GAATTC
    assert feats["P"] == (start, start + 199)


def test_result_matches_pydna_ligation_assembly():
    result = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                             [frag(ins1(), "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    backbone = max(Dseqrecord(VECTOR_SEQ, circular=True).cut(EcoRI, XhoI), key=len)
    piece1 = Dseqrecord(INS1_SEQ).cut(EcoRI, BamHI)[1]
    piece2 = Dseqrecord(INS2_SEQ).cut(BamHI, XhoI)[1]
    [reference] = ligation_assembly([backbone, piece1, piece2], circular_only=True)
    assert same_circular(result.products[0].doc["sequence"], str(reference.seq))


def test_inserts_listed_out_of_order_still_ligate_with_a_warning():
    result = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                             [frag(ins2(), "BamHI", "XhoI"), frag(ins1(), "EcoRI", "BamHI")])
    [product] = result.products
    assert same_circular(product.doc["sequence"], A + "GAATTC" + C1 + "GGATCC" + C2 + "CTCGAG" + B)
    assert "pVEC → promoter → CDS" in product.label
    assert "not ligate in the order listed" in result.warnings[0]


def test_one_enzyme_everywhere_gives_every_order_and_orientation():
    vec = SequenceDoc(**doc("pE", A + "GAATTC" + B, True))
    x = SequenceDoc(**doc("x", "CC" + "GAATTC" + C1 + "GAATTC" + "CC", False))
    y = SequenceDoc(**doc("y", "CC" + "GAATTC" + C2 + "GAATTC" + "CC", False))
    result = insert_multiple(frag(vec, "EcoRI"), [frag(x, "EcoRI"), frag(y, "EcoRI")])
    # vector fixed forward; two orders of the inserts, each insert either way round
    assert len(result.products) == 8
    seqs = [p.doc["sequence"].upper() for p in result.products]
    assert seqs[0] == A + "GAATTC" + C1 + "GAATTC" + C2 + "GAATTC" + B
    assert any(same_circular(s, A + "GAATTC" + rc(C2) + "GAATTC" + C1 + "GAATTC" + B) for s in seqs)
    assert any("8 different products" in w for w in result.warnings)
    assert any("close without an insert" in w and "Dephosphorylate" in w for w in result.warnings)
    assert result.products[1].doc["name"] == "pE + x + y, product 2"
    assert len({p.label for p in result.products}) == 8


def test_pcr_product_with_the_same_site_at_both_ends_is_not_warned_about():
    """EcoRI at each end of a PCR product cuts it twice, one cut per end: nothing to check.
    A third EcoRI site inside it is worth a warning."""
    vec = SequenceDoc(**doc("pE", A + "GAATTC" + B, True))
    pcr = SequenceDoc(**doc("pcr", "CC" + "GAATTC" + C1 + "GAATTC" + "CC", False))
    result = insert_multiple(frag(vec, "EcoRI"), [frag(pcr, "EcoRI")])
    assert result.products[0].doc["sequence"].upper() == A + "GAATTC" + C1 + "GAATTC" + B
    assert not any("cuts" in w and "times" in w for w in result.warnings)

    inner = SequenceDoc(**doc("pcr", "CC" + "GAATTC" + C1 + "GAATTC" + C2 + "GAATTC" + "CC",
                              False))
    result = insert_multiple(frag(vec, "EcoRI"), [frag(inner, "EcoRI")])
    assert any(w.startswith("EcoRI cuts “pcr” 3 times; check that the fragment used")
               for w in result.warnings)


def test_details_count_inserts_and_fragments_in_plain_words():
    vec = SequenceDoc(**doc("pE", A + "GAATTC" + B, True))
    pcr = SequenceDoc(**doc("pcr", "CC" + "GAATTC" + C1 + "GAATTC" + "CC", False))
    one = insert_multiple(frag(vec, "EcoRI"), [frag(pcr, "EcoRI")]).products[0]
    assert one.details.startswith("Ligated 1 insert into pE (cut with EcoRI): ")
    two = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                          [frag(ins1(), "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    assert two.products[0].details.startswith("Ligated 2 inserts into pVEC")
    assert "(s)" not in one.details + two.products[0].details


def test_directional_insert_with_an_extra_site_is_warned_about():
    """EcoRI–BamHI from a PCR product: one EcoRI end is expected; a second EcoRI site is not
    (the old count of distinct enzymes let it through)."""
    extra = SequenceDoc(**doc("promoter", "CC" + "GAATTC" + C2 + "GAATTC" + C1 + "GGATCC" + "CC",
                              False))
    result = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                             [frag(extra, "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    assert any(w.startswith("EcoRI cuts “promoter” 2 times") for w in result.warnings)
    clean = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                            [frag(ins1(), "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    assert not any("times;" in w for w in clean.warnings)


def test_blunt_pcr_insert_into_ecorv_cut_vector_both_orientations():
    vec = SequenceDoc(**doc("pBlunt", A + "GATATC" + B, True))
    pcr = SequenceDoc(**doc("amplicon", C1, False))
    result = insert_multiple(frag(vec, "EcoRV"), [FragmentIn(doc=pcr)])
    seqs = [p.doc["sequence"].upper() for p in result.products]
    # GAT^ATC is blunt: the PCR product sits between GAT and ATC, either way round.
    assert seqs == [A + "GAT" + C1 + "ATC" + B, A + "GAT" + rc(C1) + "ATC" + B]
    assert "amplicon (reversed)" in result.products[1].label
    text = " ".join(result.warnings)
    assert "T4 polynucleotide kinase" in text
    assert "Blunt ends ligate less efficiently" in text
    assert "close without an insert" in text


def test_vector_made_by_inverse_pcr_is_used_whole():
    linear_vector = SequenceDoc(**doc("pPCR", B + A, False))
    result = insert_multiple(FragmentIn(doc=linear_vector), [FragmentIn(doc=SequenceDoc(
        **doc("amplicon", C1, False)))])
    assert result.products[0].doc["sequence"].upper() == B + A + C1
    assert any("phosphorylate only the insert" in w for w in result.warnings)


def test_region_of_a_document_is_used_as_a_pcr_product():
    vec = SequenceDoc(**doc("pBlunt", A + "GATATC" + B, True))
    template = SequenceDoc(**doc("template", C2 + C1 + C2, True))
    start, end = len(C2), len(C2) + len(C1) - 1
    result = insert_multiple(frag(vec, "EcoRV"), [FragmentIn(doc=template, start=start, end=end)])
    assert result.products[0].doc["sequence"].upper() == A + "GAT" + C1 + "ATC" + B
    assert "template 201..500 (blunt PCR product, 300 bp)" in result.products[0].details


def test_incompatible_ends_are_an_error_naming_the_junction():
    bad = SequenceDoc(**doc("bad", "CC" + "GAATTC" + C1 + "CTGCAG" + "CC", False))
    with pytest.raises(CloningError) as err:
        insert_multiple(frag(vector(), "EcoRI", "XhoI"), [frag(bad, "EcoRI", "PstI")])
    message = str(err.value)
    assert "cannot ligate into a plasmid" in message
    assert "“bad” ends with 3' overhang TGCA" in message
    assert "“pVEC” starts with" in message


def test_explicit_fragment_choice_and_out_of_range_index():
    donor = SequenceDoc(**doc("pDonor", clean_dna(2500, 105) + "GAATTC" + C1 + "GGATCC", True))
    # The digest of the donor lists the 2500 bp backbone and the insert; pick the insert by index.
    result = insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                             [frag(donor, "EcoRI", "BamHI", fragment_index=0),
                              frag(ins2(), "BamHI", "XhoI")])
    assert result.products[0].doc["sequence"].upper() == (
        A + "GAATTC" + C1 + "GGATCC" + C2 + "CTCGAG" + B)
    with pytest.raises(CloningError, match="Fragment 5 of “pDonor” does not exist"):
        insert_multiple(frag(vector(), "EcoRI", "XhoI"),
                        [frag(donor, "EcoRI", "BamHI", fragment_index=4)])


def test_circular_vector_needs_enzymes_and_circular_insert_needs_a_region():
    with pytest.raises(CloningError, match="Choose the enzyme"):
        insert_multiple(FragmentIn(doc=vector()), [frag(ins1(), "EcoRI", "BamHI")])
    with pytest.raises(CloningError, match="is circular: choose enzymes"):
        insert_multiple(frag(vector(), "EcoRI", "XhoI"), [FragmentIn(doc=vector())])


def test_enzyme_that_does_not_cut_an_insert_is_an_error():
    with pytest.raises(CloningError, match="XhoI does not cut “promoter”"):
        insert_multiple(frag(vector(), "EcoRI", "XhoI"), [frag(ins1(), "EcoRI", "XhoI")])


# ---- Ligate Fragments ----------------------------------------------------------------------

def test_two_sticky_fragments_join_into_one_linear_product():
    result = ligate_fragments([frag(ins1(), "EcoRI", "BamHI"), frag(ins2(), "BamHI", "XhoI")])
    [product] = result.products
    # Each piece keeps its single-stranded ends: AATT on the left, TCGA on the right.
    assert product.doc["sequence"].upper() == "AATTC" + C1 + "GGATCC" + C2 + "CTCGA"
    assert product.doc["circular"] is False
    assert product.operation == "Ligate fragments"
    assert product.parents == ["promoter", "CDS"]
    assert "Ends: 5' overhang AATT / 5' overhang TCGA" in product.details
    assert not any("close into a circle" in w for w in result.warnings)


def test_two_blunt_pcr_products_give_four_linear_molecules():
    p1 = SequenceDoc(**doc("left", C1, False))
    p2 = SequenceDoc(**doc("right", C2, False))
    result = ligate_fragments([FragmentIn(doc=p1), FragmentIn(doc=p2)])
    seqs = [p.doc["sequence"].upper() for p in result.products]
    assert seqs[0] == C1 + C2
    assert len(seqs) == 4
    molecules = {min(s, rc(s)) for s in seqs}
    expected = {C1 + C2, C2 + C1, C1 + rc(C2), rc(C1) + C2}
    assert molecules == {min(s, rc(s)) for s in expected}
    text = " ".join(result.warnings)
    assert "4 different products" in text
    assert "close into a circle" in text
    assert "T4 polynucleotide kinase" in text


def test_a_single_fragment_can_be_closed_into_a_circle():
    piece = SequenceDoc(**doc("frag", "CC" + "GAATTC" + C1 + "GAATTC" + "CC", False))
    result = ligate_fragments([frag(piece, "EcoRI")], circular=True)
    [product] = result.products
    assert product.doc["circular"] is True
    assert same_circular(product.doc["sequence"], "GAATTC" + C1)


def test_linear_ligation_needs_two_fragments():
    with pytest.raises(CloningError, match="Add a second fragment"):
        ligate_fragments([frag(ins1(), "EcoRI", "BamHI")])


def test_fragments_that_cannot_join_are_an_error():
    bad = SequenceDoc(**doc("bad", "CC" + "GAATTC" + C1 + "CTGCAG" + "CC", False))
    with pytest.raises(CloningError) as err:
        ligate_fragments([frag(ins2(), "BamHI", "XhoI"), frag(bad, "EcoRI", "PstI")])
    assert "cannot ligate into a linear molecule" in str(err.value)
    assert "“CDS” ends with 5' overhang TCGA but “bad” starts with 5' overhang AATT" in str(err.value)


def test_two_copies_of_a_fragment_join_through_matching_ends():
    # EcoRI–BamHI twice: the BamHI ends pair, giving a head-to-head dimer.
    result = ligate_fragments([frag(ins1(), "EcoRI", "BamHI"), frag(ins1(), "EcoRI", "BamHI")])
    seqs = {p.doc["sequence"].upper() for p in result.products}
    head_to_head = "AATTC" + C1 + "GGATCC" + rc(C1) + "GAATT"
    tail_to_tail = "GATCC" + rc(C1) + "GAATTC" + C1 + "G" + "GATC"
    assert seqs <= {head_to_head, rc(head_to_head), tail_to_tail, rc(tail_to_tail)}
    assert {min(s, rc(s)) for s in seqs} == {min(head_to_head, rc(head_to_head)),
                                             min(tail_to_tail, rc(tail_to_tail))}
    assert result.products[0].parents == ["promoter"]


def test_too_many_blunt_fragments_are_refused_with_advice():
    pieces = [FragmentIn(doc=SequenceDoc(**doc(f"p{i}", clean_dna(60, 110 + i), False)))
              for i in range(5)]
    with pytest.raises(CloningError, match="more than 50 different ways"):
        ligate_fragments(pieces)


# ---- HTTP ----------------------------------------------------------------------------------

def test_http_routes():
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    body = {"vector": {"doc": doc("pVEC", VECTOR_SEQ, True), "enzymes": ["EcoRI", "XhoI"]},
            "inserts": [{"doc": doc("promoter", INS1_SEQ, False), "enzymes": ["EcoRI", "BamHI"]},
                        {"doc": doc("CDS", INS2_SEQ, False), "enzymes": ["BamHI", "XhoI"]}]}
    res = client.post("/api/cloning/insert-multiple", json=body)
    assert res.status_code == 200
    [product] = res.json()["products"]
    assert product["label"].startswith("pVEC → promoter → CDS")
    res = client.post("/api/cloning/ligate", json={"fragments": body["inserts"]})
    assert res.status_code == 200
    assert res.json()["products"][0]["doc"]["sequence"].upper().startswith("AATTC" + C1[:20])
