"""MultiSite Gateway: two to four entry clones in one LR reaction, BP into a P4-P1R donor.

Test data — public GenBank records fetched from NCBI with E-utilities, saved unchanged:

- ``tests/data/cloning-multisite-5prime.gb`` — LT726877.1, "Vector pEntry-L4
  rtTA-IRES-Puro-pA-Ins/Ins-TRE R1", 6006 bp, attL4 and attR1 annotated (BCCM/LMBP plasmid
  collection, Ghent).
- ``tests/data/cloning-gateway-entry.gb`` — LT726670.1, entry clone pENTR207-CATMA5a57070,
  3561 bp, attL1 and attL2 (the middle entry clone; also used by ``test_cloning_gateway.py``).
- ``tests/data/cloning-multisite-3prime.gb`` — LT726875.1, "Vector
  pEntry-attR2-IRES-eGFP-luc+-pA-attL3", 6273 bp, attR2 and attL3 annotated (BCCM/LMBP).
- ``tests/data/cloning-multisite-destination.gb`` — LT726807.1, "Mammalian vector
  pdcDNA-R4R3-dhCMV", 6125 bp, attR4 and attR3 around a ccdB cassette (BCCM/LMBP).
- ``tests/data/cloning-multisite-destination-r4r2.gb`` — AP018631.1, "Gateway vector
  R4pMpGWB495", 11750 bp, attR4 and attR2: a destination vector for two entry clones (Mano
  et al., "Novel gateway binary vectors for rapid tripartite DNA assembly").
- ``tests/data/cloning-multisite-entry-l1r5.gb`` — PP281278.1, "Genome editing vector
  pGN1101CK", 4521 bp, attL1 and attR5 annotated (Wu et al.).
- ``tests/data/cloning-multisite-entry-l5l4.gb`` — PP281279.1, "Genome editing vector
  pGN1102CK", 4430 bp, attL5 and attL4 annotated; used here as the published source of an
  attL5 site.
- ``tests/data/cloning-multisite-expression.gb`` — KC800697.1, "Expression vector
  mTol2-B-P-Gal4VP16-CG3", 6522 bp: a **published three-fragment MultiSite product**, with its
  attB4, attB1, attB2 and attB3 sites annotated (Otsuna et al., zebrafish Gal4 enhancer trap).
- ``tests/data/cloning-multisite-donor.gb`` — HM623915.1, "Marker vector pA-hyg-OSCAR",
  4833 bp, attP4 and attP1r around a hygromycin cassette (Paz et al., "One Step Construction
  of Agrobacterium-Recombination-ready-plasmids").
- ``tests/data/cloning-multisite-entry-l4r1.gb`` — LC477681.1, "Gateway adapted binary vector
  pCas9_P4P1r", 8581 bp: a published clone carrying the att site pair a P4-P1R donor leaves
  behind (the record names them attR4/attR1; by the site definitions the first one is an
  attL4, which is what the detection reports).

No published record gives a MultiSite expression clone together with all of its parents, so
the expected products here are **built by hand** from the parent records and the published att
core sequences (GTACAAA, GTACAAG, GTATAAT, GTATAGA, GTATACA for sites 1 to 5, each inside the
13 bp recombination region the two partners share): every molecule contributes the arc between
its two att cores that carries a named feature of that parent (the insert of an entry clone,
the backbone of a vector), and the arcs are joined at the cores. That expectation is
independent of the assembler.

What a published record does check is the junctions: the attB4, attB1, attB2 and attB3 sites
of the simulated three-entry product are the same sites, base for base, as the ones annotated
in the published three-fragment clone KC800697.1.
"""

import re
from pathlib import Path

import pytest
from Bio import SeqIO
from cloning_helpers import clean_dna, doc, rc
from fastapi import FastAPI
from fastapi.testclient import TestClient
from pydna.dseqrecord import Dseqrecord

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.gateway import gateway, site_names
from vettore.fileio import read_documents

DATA = Path(__file__).parent / "data"

# The 7 bp core of each Gateway att site, as it reads in the attB sites of the published
# records (see the attB4 and attB3 annotations of KC800697.1).
CORES = {1: "GTACAAA", 2: "GTACAAG", 3: "GTATAAT", 4: "GTATAGA", 5: "GTATACA"}

# The attB sites annotated in the published three-fragment clone KC800697.1.
PUBLISHED_SITES = {"attB4": "CAACTTTGTATAGAAAAGTTG",
                   "attB1": "CAAGTTTGTACAAAAAAGCAGGCT",
                   "attB2": "ACCCAGCTTTCTTGTACAAAGTGG",
                   "attB3": "CAACTTTATTATACATAGTTG"}


def load(stem: str) -> SequenceDoc:
    path = DATA / f"{stem}.gb"
    [document] = read_documents(path.name, path.read_bytes())
    document["name"] = stem
    return SequenceDoc(**document)


def record(stem: str) -> Dseqrecord:
    return Dseqrecord(SeqIO.read(DATA / f"{stem}.gb", "genbank"), circular=True)


def sequence(stem: str) -> str:
    return str(record(stem).seq).upper()


def marker(stem: str, note: str) -> str:
    """60 bp of a feature annotated in a record, to say which arc of it the reaction keeps."""
    annotated = SeqIO.read(DATA / f"{stem}.gb", "genbank")
    for feature in annotated.features:
        text = " ".join(value for key in ("note", "label", "gene", "product")
                        for value in feature.qualifiers.get(key, []))
        if note.lower() in text.lower():
            return str(feature.extract(annotated.seq)).upper()[:60]
    raise AssertionError(f"no feature “{note}” in {stem}")


def recombination_region(seq: str, number: int) -> tuple[int, int]:
    """The single 13 bp region (T W T + core + A A A) where att{number} recombines."""
    core = CORES[number]
    found = []
    for pattern in (re.compile(f"T[AT]T{core}AAA"), re.compile(f"TTT{rc(core)}A[AT]A")):
        found += [(match.start(), match.end()) for match in pattern.finditer(seq + seq)
                  if match.start() < len(seq)]
    assert len(found) == 1, f"{len(found)} att{number} recombination regions"
    return found[0]


def contributed(seq: str, first: int, second: int, note: str) -> str:
    """What a molecule puts into the product: the arc from the recombination region of
    att{first} to the one of att{second} (not included) that carries the given feature."""
    n = len(seq)
    start, start_end = recombination_region(seq, first)
    second_start, second_end = recombination_region(seq, second)
    forward = (seq * 2)[start: start + (second_start - start) % n]
    backward = rc((seq * 2)[second_end: second_end + (start_end - second_end) % n])
    arcs = [arc for arc in (forward, backward) if note in arc or rc(note) in arc]
    assert len(arcs) == 1, "the feature must be in one of the two arcs"
    return arcs[0]


def joined(*arcs: str) -> str:
    return "".join(arcs)


def is_circular_match(product, expected: str) -> bool:
    found = product.doc["sequence"].upper()
    return (product.doc["circular"] and len(found) == len(expected)
            and (found in expected * 2 or rc(found) in expected * 2))


def entry_l5_l2() -> SequenceDoc:
    """An attL5 - insert - attL2 entry clone, built from two public records: the published
    middle entry clone LT726670.1 with its attL1 site replaced by the attL5 of PP281279.1.
    No deposited record pairs an attL5 with an attL2, and a four-fragment reaction needs one;
    both att sites and the backbone are published sequence, only the join is ours."""
    middle = SeqIO.read(DATA / "cloning-gateway-entry.gb", "genbank")
    donor = SeqIO.read(DATA / "cloning-multisite-entry-l5l4.gb", "genbank")
    [attl1] = [f for f in middle.features if "attL1" in " ".join(f.qualifiers.get("note", []))]
    [attl5] = [f for f in donor.features if " ".join(f.qualifiers.get("note", [])) == "attL5"]
    text = str(middle.seq).upper()
    built = (text[:int(attl1.location.start)] + str(attl5.extract(donor.seq)).upper()
             + text[int(attl1.location.end):])
    return SequenceDoc(name="entry-l5-l2", sequence=built, circular=True)


# ---- LR with two, three and four entry clones ----------------------------------------------

def test_three_entry_lr_builds_the_hand_built_expression_clone():
    """attR4 destination x (attL4-attR1) x (attL1-attL2) x (attR2-attL3): the product must be
    the four arcs joined at the att cores, base for base."""
    expected = joined(
        contributed(sequence("cloning-multisite-5prime"), 4, 1,
                    marker("cloning-multisite-5prime", "rtTA")),
        contributed(sequence("cloning-gateway-entry"), 1, 2,
                    marker("cloning-gateway-entry", "CATMA5a57070")),
        contributed(sequence("cloning-multisite-3prime"), 2, 3,
                    marker("cloning-multisite-3prime", "EGFP")),
        contributed(sequence("cloning-multisite-destination"), 3, 4,
                    marker("cloning-multisite-destination", "ampicillin")))
    result = gateway("LR", load("cloning-multisite-5prime"), load("cloning-gateway-entry"),
                     load("cloning-multisite-3prime"), load("cloning-multisite-destination"))
    clone = result.products[0]
    assert clone.label == f"Expression clone — {len(expected)} bp circular"
    assert len(expected) == 11974
    assert is_circular_match(clone, expected)
    assert clone.parents == ["cloning-multisite-5prime", "cloning-gateway-entry",
                             "cloning-multisite-3prime", "cloning-multisite-destination"]
    assert "3 entry clones" in clone.details
    assert [p.label.split(" — ")[0] for p in result.products] == ["Expression clone",
                                                                  "attP by-product"]


def test_three_entry_product_has_the_four_published_attb_junctions():
    """Every junction the simulation builds must be the att site a published three-fragment
    clone carries (KC800697.1), base for base and once each."""
    result = gateway("LR", load("cloning-multisite-5prime"), load("cloning-gateway-entry"),
                     load("cloning-multisite-3prime"), load("cloning-multisite-destination"))
    clone = Dseqrecord(result.products[0].doc["sequence"], circular=True)
    assert site_names(clone) == ["attB1", "attB2", "attB3", "attB4"]
    published = sequence("cloning-multisite-expression")
    text = str(clone.seq).upper()
    for name, site in PUBLISHED_SITES.items():
        assert published.count(site) + published.count(rc(site)) == 1, name
        assert (text * 2).count(site) + (text * 2).count(rc(site)) in (1, 2), name
        assert site in text * 2 or rc(site) in text * 2, name


def test_two_entry_lr_into_a_destination_for_two_fragments():
    """A 5' entry clone and a middle entry clone into an attR4/attR2 destination vector."""
    expected = joined(
        contributed(sequence("cloning-multisite-5prime"), 4, 1,
                    marker("cloning-multisite-5prime", "rtTA")),
        contributed(sequence("cloning-gateway-entry"), 1, 2,
                    marker("cloning-gateway-entry", "CATMA5a57070")),
        contributed(sequence("cloning-multisite-destination-r4r2"), 2, 4,
                    marker("cloning-multisite-destination-r4r2", "nptII")))
    result = gateway("LR", load("cloning-multisite-5prime"), load("cloning-gateway-entry"),
                     load("cloning-multisite-destination-r4r2"))
    clone = result.products[0]
    assert clone.label == f"Expression clone — {len(expected)} bp circular"
    assert is_circular_match(clone, expected)
    assert site_names(Dseqrecord(clone.doc["sequence"], circular=True)) == [
        "attB1", "attB2", "attB4"]
    assert clone.doc["name"] == ("cloning-multisite-destination-r4r2 + "
                                 "cloning-multisite-5prime + cloning-gateway-entry")


def test_four_entry_lr_through_the_att5_sites():
    """Four entry clones chained 4 - 1 - 5 - 2 into an attR4/attR3 destination vector."""
    built = entry_l5_l2()
    expected = joined(
        contributed(sequence("cloning-multisite-5prime"), 4, 1,
                    marker("cloning-multisite-5prime", "rtTA")),
        contributed(sequence("cloning-multisite-entry-l1r5"), 1, 5,
                    marker("cloning-multisite-entry-l1r5", "gRNA scaffold")),
        contributed(built.sequence, 5, 2, marker("cloning-gateway-entry", "CATMA5a57070")),
        contributed(sequence("cloning-multisite-3prime"), 2, 3,
                    marker("cloning-multisite-3prime", "EGFP")),
        contributed(sequence("cloning-multisite-destination"), 3, 4,
                    marker("cloning-multisite-destination", "ampicillin")))
    result = gateway("LR", load("cloning-multisite-5prime"),
                     load("cloning-multisite-entry-l1r5"), built,
                     load("cloning-multisite-3prime"), load("cloning-multisite-destination"))
    clone = result.products[0]
    assert clone.label == f"Expression clone — {len(expected)} bp circular"
    assert is_circular_match(clone, expected)
    assert site_names(Dseqrecord(clone.doc["sequence"], circular=True)) == [
        "attB1", "attB2", "attB3", "attB4", "attB5"]
    assert "4 entry clones" in clone.details


def test_the_ccdb_of_an_entry_clone_is_reported_by_name():
    """pGN1101CK carries a ccdB stuffer between its att sites, so the product carries it too:
    the warning must say which molecule it came from, not blame the destination vector."""
    result = gateway("LR", load("cloning-multisite-5prime"),
                     load("cloning-multisite-entry-l1r5"), entry_l5_l2(),
                     load("cloning-multisite-3prime"), load("cloning-multisite-destination"))
    [warning] = [w for w in result.warnings if "ccdB" in w]
    assert "“cloning-multisite-entry-l1r5”" in warning
    assert "clone your insert into that entry clone first" in warning


def test_a_five_entry_reaction_is_refused():
    entries = [load("cloning-multisite-5prime")] * 5
    with pytest.raises(CloningError) as err:
        gateway("LR", *entries, load("cloning-multisite-destination"))
    assert "at most 4 entry clones" in str(err.value)


# ---- messages when the sites do not chain --------------------------------------------------

def test_a_missing_middle_entry_clone_names_the_sites_that_cannot_pair():
    """5' entry (attL4, attR1) + 3' entry (attR2, attL3) with no attL1/attL2 clone between
    them: the message must name the sites of each molecule and the junction that is missing."""
    with pytest.raises(CloningError) as err:
        gateway("LR", load("cloning-multisite-5prime"), load("cloning-multisite-3prime"),
                load("cloning-multisite-destination"))
    message = str(err.value)
    assert "do not pair" in message
    assert "“cloning-multisite-5prime”: attL4, attR1" in message
    assert "“cloning-multisite-3prime”: attL3, attR2" in message
    assert "attR1 in “cloning-multisite-5prime” has no attL1 in the other molecules" in message
    assert "attR2 in “cloning-multisite-3prime” has no attL2 in the other molecules" in message


def test_entry_clones_in_the_wrong_order_still_react():
    """The order of the entry clones in the list does not change the biology: the reaction
    finds the chain. Only the destination vector must come last."""
    straight = gateway("LR", load("cloning-multisite-5prime"), load("cloning-gateway-entry"),
                       load("cloning-multisite-3prime"),
                       load("cloning-multisite-destination")).products[0]
    shuffled = gateway("LR", load("cloning-multisite-3prime"), load("cloning-multisite-5prime"),
                       load("cloning-gateway-entry"),
                       load("cloning-multisite-destination")).products[0]
    assert is_circular_match(shuffled, straight.doc["sequence"].upper())


def test_a_destination_vector_used_as_an_entry_clone_is_reported():
    with pytest.raises(CloningError) as err:
        gateway("LR", load("cloning-multisite-destination"), load("cloning-multisite-5prime"),
                load("cloning-multisite-destination"))
    message = str(err.value)
    assert "entry clone 1 (attL or attR sites)" in message
    assert "has no attL site" in message
    assert "entry clones and attR sites in the destination vector" in message


# ---- BP into a MultiSite donor vector ------------------------------------------------------

def test_bp_into_a_p4_p1r_donor_gives_an_attl4_attr1_entry_clone():
    """An attB4 - insert - attB1 product recombined with a donor whose attP4 and attP1r sites
    read the same way round: the entry clone keeps the piece between the donor's att sites
    (here the hygromycin cassette of pA-hyg-OSCAR) and comes out with attL4 and attR1, the
    pair the published P4-P1R clone LC477681.1 carries."""
    cassette = attb4_attb1_product()
    expected = joined(
        contributed(cassette.sequence, 4, 1, marker("cloning-multisite-5prime", "rtTA")),
        contributed(sequence("cloning-multisite-donor"), 1, 4,
                    marker("cloning-multisite-donor", "hph")))
    result = gateway("BP", cassette, load("cloning-multisite-donor"))
    clone = result.products[0]
    assert clone.label == f"Entry clone — {len(expected)} bp circular"
    assert is_circular_match(clone, expected)
    assert site_names(Dseqrecord(clone.doc["sequence"], circular=True)) == ["attL4", "attR1"]
    assert site_names(record("cloning-multisite-entry-l4r1")) == ["attL4", "attR1"]


def attb4_attb1_product() -> SequenceDoc:
    """The attB4 - 5' element - attB1 piece of the hand-built three-entry expression clone,
    with 10 bp of primer tail on each side: what a PCR with attB4 and attB1r primers gives."""
    clone = joined(
        contributed(sequence("cloning-multisite-5prime"), 4, 1,
                    marker("cloning-multisite-5prime", "rtTA")),
        contributed(sequence("cloning-gateway-entry"), 1, 2,
                    marker("cloning-gateway-entry", "CATMA5a57070")),
        contributed(sequence("cloning-multisite-3prime"), 2, 3,
                    marker("cloning-multisite-3prime", "EGFP")),
        contributed(sequence("cloning-multisite-destination"), 3, 4,
                    marker("cloning-multisite-destination", "ampicillin")))
    doubled = clone + clone
    start = doubled.index(PUBLISHED_SITES["attB4"])
    end = doubled.index(PUBLISHED_SITES["attB1"], start) + len(PUBLISHED_SITES["attB1"])
    return SequenceDoc(name="attB4 PCR product", sequence=doubled[start - 10: end + 10],
                       circular=False)


def test_a_bp_with_three_molecules_is_refused():
    with pytest.raises(CloningError) as err:
        gateway("BP", attb4_attb1_product(), attb4_attb1_product(),
                load("cloning-multisite-donor"))
    assert "one attB-flanked molecule and one donor vector" in str(err.value)


def test_the_att_pair_chosen_for_the_donor_is_checked():
    """The BP form offers the att pairs of the donor; a pair the donor does not carry is a
    mistake the user can correct."""
    with pytest.raises(CloningError) as err:
        gateway("BP", attb4_attb1_product(), load("cloning-multisite-donor"),
                donor_sites=[1, 2])
    message = str(err.value)
    assert "does not carry the att pair you chose (attP1 and attP2)" in message
    assert "“cloning-multisite-donor”: attP1, attP4" in message
    result = gateway("BP", attb4_attb1_product(), load("cloning-multisite-donor"),
                     donor_sites=[4, 1])
    assert result.products[0].label.startswith("Entry clone")


def test_a_donor_vector_reads_as_attp_sites_only():
    """attL and attR are halves of an attP: a donor vector must read attP1, attP4, not six
    site names, or the dialog offers it as an entry clone."""
    assert site_names(record("cloning-multisite-donor")) == ["attP1", "attP4"]


# ---- HTTP ----------------------------------------------------------------------------------

def client() -> TestClient:
    app = FastAPI()
    app.include_router(router)
    return TestClient(app)


def test_http_route_takes_three_entry_clones():
    res = client().post("/api/cloning/gateway", json={
        "reaction": "LR",
        "inputs": [load(stem).model_dump() for stem in (
            "cloning-multisite-5prime", "cloning-gateway-entry", "cloning-multisite-3prime",
            "cloning-multisite-destination")]})
    assert res.status_code == 200
    payload = res.json()
    assert payload["products"][0]["label"] == "Expression clone — 11974 bp circular"


def test_http_route_rejects_a_single_molecule():
    with pytest.raises(CloningError) as err:      # the app turns this into HTTP 422
        client().post("/api/cloning/gateway", json={
            "reaction": "LR", "inputs": [load("cloning-gateway-entry").model_dump()]})
    assert "Choose both molecules" in str(err.value)


def test_gateway_sites_route_reports_donor_and_phic31_sites():
    res = client().post("/api/cloning/gateway-sites", json={"docs": [
        load("cloning-multisite-5prime").model_dump(),
        load("cloning-multisite-donor").model_dump(),
        load("cloning-phic31-attb").model_dump(),
        doc("plain", clean_dna(500, 7), True)]})
    assert res.status_code == 200
    assert res.json() == {
        "sites": [["attL4", "attR1"], ["attP1", "attP4"], [], []],
        "other": [[], [], ["phiC31 attB"], []]}
