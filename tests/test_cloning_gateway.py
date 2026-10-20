"""Gateway BP and LR, validated against public records that carry annotated att sites.

Test data (fetched from NCBI with E-utilities, saved unchanged):

- ``tests/data/cloning-gateway-entry.gb`` — GenBank LT726670.1, "Gateway entry vector
  pENTR207-CATMA5a57070", 3561 bp, attL1 and attL2 annotated (BCCM/LMBP collection; the GST
  library is described in Hilson et al. 2004, Genome Res 14:2176).
- ``tests/data/cloning-gateway-destination.gb`` — GenBank AY568055.1, "Cloning vector
  pAGRIKOLA", 10459 bp, two attR1/attR2 cassettes annotated (Hilson et al. 2004).
- ``tests/data/cloning-gateway-expression.gb`` — GenBank LT725632.1, "Gateway expression vector
  pAGRIKOLA-CATMA5a57070", 9443 bp: the published expression clone of that entry clone, with
  its attB1 and attB2 sites annotated.

The expression clone was built from the same entry clone and a pAGRIKOLA preparation whose
backbone differs from the AY568055 deposit (different submitters, different annotation), so the
whole plasmids are not expected to be identical; what the simulation must reproduce is the
recombined cassette, attB1 to attB2, which it does base for base. The BP reaction run back on
that product regenerates the published entry clone exactly.
"""

from pathlib import Path

import pytest
from Bio import SeqIO
from cloning_helpers import clean_dna, doc, rc, same_circular
from fastapi import FastAPI
from fastapi.testclient import TestClient
from pydna.dseqrecord import Dseqrecord
from pydna.gateway import find_gateway_sites

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.gateway import gateway
from vettore.fileio import read_documents

DATA = Path(__file__).parent / "data"


def load(stem: str) -> SequenceDoc:
    path = DATA / f"cloning-gateway-{stem}.gb"
    [document] = read_documents(path.name, path.read_bytes())
    document["name"] = stem
    return SequenceDoc(**document)


def record(stem: str) -> Dseqrecord:
    path = DATA / f"cloning-gateway-{stem}.gb"
    return Dseqrecord(SeqIO.read(path, "genbank"), circular=True)


def cassette(record: Dseqrecord, first: str, second: str) -> str:
    """The recombined cassette: first site to second site, on the strand of the first site."""
    sites = find_gateway_sites(record, False)
    sequence = str(record.seq).upper()
    n = len(sequence)
    found = []
    for a in sites.get(first, []):
        if int(a.end) - int(a.start) > 300:
            continue
        for b in sites.get(second, []):
            if int(b.end) - int(b.start) > 300:
                continue
            if a.strand == 1 and b.strand == -1:
                start, end = int(a.start), int(b.end)
            elif a.strand == -1 and b.strand == 1:
                start, end = int(b.start), int(a.end)
            else:
                continue
            length = (end - start) % n
            if 0 < length < 2000:
                piece = (sequence * 2)[start: start + length]
                found.append(piece if a.strand == 1
                             else str(Dseqrecord(piece).reverse_complement().seq).upper())
    assert found, f"no {first}..{second} cassette in {record.name}"
    return max(set(found), key=found.count)


def occurrences(record: Dseqrecord, probe: str) -> int:
    """How many times a piece appears in a (circular) molecule, either strand."""
    sequence = str(record.seq).upper()
    window = sequence + sequence[: len(probe) - 1] if record.circular else sequence
    reverse = str(Dseqrecord(probe).reverse_complement().seq).upper()
    return window.count(probe) + window.count(reverse)


def product_record(product) -> Dseqrecord:
    return Dseqrecord(product.doc["sequence"], circular=product.doc["circular"])


def by_product(result) -> SequenceDoc:
    """The attP molecule the LR reaction releases, as an input for the BP reaction back."""
    [found] = [p for p in result.products if p.label.startswith("attP by-product")
               and find_gateway_sites(product_record(p), False).get("attP1")][:1]
    return SequenceDoc(**{**found.doc, "name": "attP by-product"})


# ---- LR: entry clone x destination vector --------------------------------------------------

def test_lr_reproduces_the_published_expression_clone_cassette():
    """The recombined attB1 - tag - attB2 cassette must match the published clone exactly."""
    result = gateway("LR", load("entry"), load("destination"))
    expression = result.products[0]
    assert expression.operation == "Gateway LR reaction"
    assert expression.label.startswith("Expression clone")
    assert expression.parents == ["entry", "destination"]
    assert expression.doc["circular"] is True
    published = cassette(record("expression"), "attB1", "attB2")
    assert len(published) == 239     # attB1 + the gene-specific tag + attB2
    assert occurrences(product_record(expression), published) == 1


def test_a_second_lr_fills_the_other_cassette_and_gives_the_published_hairpin():
    """pAGRIKOLA takes the tag twice, in opposite orientations; the published clone carries the
    cassette twice, and so does the product of two rounds of the reaction."""
    first = gateway("LR", load("entry"), load("destination")).products[0]
    product = SequenceDoc(**{**first.doc, "name": "round 1"})
    second = gateway("LR", load("entry"), product).products[0]
    published = cassette(record("expression"), "attB1", "attB2")
    assert occurrences(record("expression"), published) == 2
    assert occurrences(product_record(second), published) == 2
    assert not any("ccdB" in w for w in gateway("LR", load("entry"), product).warnings)


def test_lr_also_reports_the_attp_by_product():
    result = gateway("LR", load("entry"), load("destination"))
    released = [p for p in result.products if p.label.startswith("attP by-product")]
    assert released, [p.label for p in result.products]
    assert all("by-product" in p.doc["name"] for p in released)
    assert find_gateway_sites(product_record(released[0]), False).get("attP1")
    # Every product opens as its own tab, so no two may share a name.
    assert len({p.doc["name"] for p in result.products}) == len(result.products)


def test_lr_into_a_two_cassette_vector_reports_both_clones_and_the_leftover_ccdb():
    """pAGRIKOLA carries two Gateway cassettes (it builds a hairpin), so one reaction fills
    one of them and the product still carries the ccdB of the other."""
    result = gateway("LR", load("entry"), load("destination"))
    clones = [p for p in result.products if p.label.startswith("Expression clone")]
    assert len(clones) == 2
    assert any("one per Gateway cassette" in w for w in result.warnings)
    assert any("ccdB" in w for w in result.warnings)


def test_lr_warns_that_single_site_products_are_not_shown():
    result = gateway("LR", load("entry"), load("destination"))
    assert any("co-integrate" in w for w in result.warnings)


# ---- BP: the reverse reaction --------------------------------------------------------------

def test_bp_regenerates_the_published_entry_clone_base_for_base():
    """attB expression clone x attP by-product -> the entry clone it came from."""
    lr = gateway("LR", load("entry"), load("destination"))
    result = gateway("BP", load("expression"), by_product(lr))
    entry_clone = result.products[0]
    assert entry_clone.operation == "Gateway BP reaction"
    assert entry_clone.label.startswith("Entry clone")
    assert product_record(entry_clone).seguid() == record("entry").seguid()


def test_bp_product_carries_the_attl_sites():
    lr = gateway("LR", load("entry"), load("destination"))
    [entry_clone, *_] = gateway("BP", load("expression"), by_product(lr)).products
    sites = find_gateway_sites(product_record(entry_clone), False)
    assert sites.get("attL1") and sites.get("attL2")


def test_bp_with_a_linear_attb_pcr_product_gives_the_entry_clone():
    """The textbook BP reaction: a linear PCR product with attB tails. Its only piece between
    the two sites is the insert; nothing may be read through its ends as if it were circular.
    attB1-tag-attB2 cut from the published expression clone, with 10 bp of primer overhang on
    each side, recombined into the attP donor, must give back the published entry clone."""
    amplicon = "GGGGCCCCGG" + cassette(record("expression"), "attB1", "attB2") + "CCGGGGCCCC"
    pcr = SequenceDoc(name="attB PCR product", sequence=amplicon, circular=False)
    lr = gateway("LR", load("entry"), load("destination"))
    result = gateway("BP", pcr, by_product(lr))
    entry_clone = result.products[0]
    assert entry_clone.label == "Entry clone — 3561 bp circular"
    assert entry_clone.doc["name"] == "attP by-product + attB PCR product"
    assert product_record(entry_clone).seguid() == record("entry").seguid()
    assert not any("right way round" in w for w in result.warnings)


def test_lr_with_an_insert_longer_than_the_entry_backbone():
    """A 5 kb ORF in a 3.2 kb entry backbone: the insert is the piece the att sites face, not
    the shorter piece. Expected product: the expression clone of the published entry clone
    with the same 5 kb inserted at the same place in its insert."""
    entry = load("entry")
    sequence = entry.sequence.upper()
    insert = cassette(record("entry"), "attL1", "attL2")[96:-96]   # between attL1 and attL2
    assert len(insert) == 197
    filler = clean_dna(5000, 2024)
    at = sequence.index(insert) + 100
    big = SequenceDoc(name="entry-big", sequence=sequence[:at] + filler + sequence[at:],
                      circular=True)
    big_insert = insert[:100] + filler + insert[100:]

    normal = product_record(gateway("LR", entry, load("destination")).products[0])
    text = str(normal.seq).upper()
    assert text.count(insert) + text.count(rc(insert)) == 1
    expected = (text.replace(insert, big_insert) if insert in text
                else text.replace(rc(insert), rc(big_insert)))

    result = gateway("LR", big, load("destination"))
    clone = result.products[0]
    assert clone.label == f"Expression clone — {len(expected)} bp circular"
    assert clone.doc["name"] == "destination + entry-big"
    assert same_circular(clone.doc["sequence"], expected)
    assert not any("right way round" in w for w in result.warnings)
    assert any("ccdB" in w for w in result.warnings)     # the second cassette is still there
    assert [p.label.split(" — ")[0] for p in result.products] == [
        "Expression clone", "Expression clone", "attP by-product", "attP by-product"]


# ---- errors a user can act on --------------------------------------------------------------

def test_molecule_without_att_sites_names_what_is_missing():
    plain = SequenceDoc(**doc("pUC-like", clean_dna(2000, 301), True))
    with pytest.raises(CloningError) as err:
        gateway("LR", plain, load("destination"))
    assert "“pUC-like” has no Gateway att site" in str(err.value)


def test_a_bp_reaction_on_an_lr_pair_is_refused_with_the_sites_found():
    with pytest.raises(CloningError) as err:
        gateway("BP", load("entry"), load("destination"))   # attL x attR is an LR reaction
    message = str(err.value)
    assert "“entry” has no attB site" in message
    assert "which recombine in an LR reaction" in message
    assert "attR1" in message and "attL1" in message


def test_wrong_reaction_between_two_entry_clones_is_an_error():
    with pytest.raises(CloningError) as err:
        gateway("LR", load("entry"), load("entry"))
    assert "no attR site" in str(err.value)


# ---- HTTP ----------------------------------------------------------------------------------

def test_http_route():
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    res = client.post("/api/cloning/gateway", json={
        "reaction": "LR",
        "inputs": [load("entry").model_dump(), load("destination").model_dump()]})
    assert res.status_code == 200
    payload = res.json()
    assert payload["products"][0]["label"].startswith("Expression clone")
    assert payload["products"][0]["doc"]["name"] == "destination + entry"


def test_gateway_sites_route_names_the_sites_of_each_document():
    """The Gateway dialogs use this to pick their default inputs. attB is the core of attL and
    attR, so an entry clone reads attL1, attL2 only."""
    app = FastAPI()
    app.include_router(router)
    client = TestClient(app)
    amplicon = "GGGGCCCCGG" + cassette(record("expression"), "attB1", "attB2") + "CCGGGGCCCC"
    res = client.post("/api/cloning/gateway-sites", json={"docs": [
        load("entry").model_dump(), load("destination").model_dump(),
        load("expression").model_dump(),
        {"name": "attB PCR product", "sequence": amplicon, "circular": False},
        doc("pUC-like", clean_dna(2000, 301), True),
        {"name": "empty", "sequence": ""}, {"name": "no sequence"}]})
    assert res.status_code == 200
    assert res.json() == {
        "sites": [["attL1", "attL2"], ["attR1", "attR2"], ["attB1", "attB2"],
                  ["attB1", "attB2"], [], [], []],
        "other": [[], [], [], [], [], [], []]}


def test_lr_messages_use_the_right_article():
    plain = SequenceDoc(**doc("pUC-like", clean_dna(2000, 301), True))
    with pytest.raises(CloningError) as err:
        gateway("LR", plain, load("destination"))
    assert "An LR reaction needs attL sites" in str(err.value)
    with pytest.raises(CloningError) as err:
        gateway("BP", plain, load("destination"))
    assert "A BP reaction needs attB sites" in str(err.value)
