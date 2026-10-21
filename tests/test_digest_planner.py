"""Diagnostic digest planner (S5.1), checked against digests worked out by hand on pBR322.

pBR322 (ENA J01749, 4,361 bp, circular) cut sites from the record's own numbering, as the
position of the last base before the top-strand cut: HindIII A^AGCTT at 29, BamHI G^GATCC at
375, SalI G^TCGAC at 651, PstI CTGCA^G at 3611, EcoRI G^AATTC at 4359.
"""

from pathlib import Path

import pytest
from Bio import SeqIO
from fastapi.testclient import TestClient

from vettore import digest_planner as dp
from vettore import gel
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.restriction import insert_fragment
from vettore.server import TOKEN_HEADER, create_app

DATA = Path(__file__).parent / "data"
PBR322 = str(SeqIO.read(DATA / "exchange-pBR322.embl", "embl").seq).upper()
PACYC184 = str(list(SeqIO.parse(DATA / "exchange-plasmids.fasta", "fasta"))[2].seq).upper()


def mol(sequence, circular=True, name="m"):
    m = dp.Molecule(name, sequence, circular)
    dp._cut_positions([m], dp.commercial_enzymes())
    return m


@pytest.mark.parametrize(("enzymes", "sizes"), [
    (("BamHI",), [4361]),                       # one cut opens the circle
    (("BamHI", "SalI"), [4085, 276]),           # 651 - 375 = 276
    (("BamHI", "PstI"), [3236, 1125]),          # 3611 - 375 = 3236; 4361 - 3236
    (("EcoRI", "HindIII"), [4330, 31]),         # 29 + (4361 - 4359) = 31
    (("HindIII", "BamHI", "SalI"), [3739, 346, 276]),  # 375 - 29; 4361 - 346 - 276
])
def test_fragment_sizes_match_the_hand_worked_pbr322_digests(enzymes, sizes):
    assert dp.fragment_sizes(mol(PBR322), enzymes) == sizes


def test_linear_molecule_keeps_its_ends():
    # Linear pBR322 cut by BamHI after base 375: 375 bp and 4,361 - 375 = 3,986 bp.
    assert dp.fragment_sizes(mol(PBR322, circular=False), ("BamHI",)) == [3986, 375]
    # No cut: the whole molecule.
    assert dp.fragment_sizes(mol(PBR322, circular=False), ("NotI",)) == [4361]


def bamhi_pcr_insert() -> SequenceDoc:
    """A 1,000 bp PCR product of pACYC184 with a BamHI site added at each end."""
    core = PACYC184[2000:3000]
    assert "GGATCC" not in core
    return SequenceDoc(name="cat PCR", sequence="CGC" + "GGATCC" + core + "GGATCC" + "CGC",
                       circular=False)


def cloning():
    vector = SequenceDoc(name="pBR322", sequence=PBR322, circular=True)
    return insert_fragment(vector, ["BamHI"], bamhi_pcr_insert(), ["BamHI"])


def test_restriction_cloning_lists_the_wrong_clones_to_screen_against():
    result = cloning()
    labels = [p.label for p in result.screen_against]
    # One BamHI site: the vector closed again without insert is the parent vector itself.
    assert labels == ["insert in the other orientation", "parent vector"]
    product = result.products[0].doc["sequence"]
    reverse, parent = (p.doc["sequence"] for p in result.screen_against)
    assert len(reverse) == len(product) == 4361 + 1006
    assert reverse.upper() != product.upper()
    assert parent.upper() == PBR322


def test_two_sites_add_the_re_ligated_empty_vector():
    vector = SequenceDoc(name="pBR322", sequence=PBR322, circular=True)
    insert = SequenceDoc(name="cat PCR", sequence="CGCGGATCC" + PACYC184[2000:3000]
                         + "GTCGACCGC", circular=False)
    result = insert_fragment(vector, ["BamHI", "SalI"], insert, ["BamHI", "SalI"])
    # BamHI and SalI ends do not pair: directional, no re-ligation, only the parent vector.
    assert [p.label for p in result.screen_against] == ["parent vector"]
    # BamHI and BglII ends pair (GATC): with the SalI site of pBR322 changed into a BglII site,
    # the cut vector can close on itself and lose the 276 bp between the two sites.
    salI = PBR322.index("GTCGAC")
    bglii_vector = SequenceDoc(name="pBR322-BglII", circular=True,
                               sequence=PBR322[:salI] + "AGATCT" + PBR322[salI + 6:])
    insert = SequenceDoc(name="cat PCR", sequence="CGCGGATCC" + PACYC184[2000:3000]
                         + "AGATCTCGC", circular=False)
    result = insert_fragment(bglii_vector, ["BamHI", "BglII"], insert, ["BamHI", "BglII"])
    labels = [p.label for p in result.screen_against]
    assert labels == ["insert in the other orientation", "empty vector, re-ligated",
                      "parent vector"]
    assert len(result.screen_against[1].doc["sequence"]) == 4361 - 276


def plan_for_cloning(**kwargs):
    result = cloning()
    expected = ("Expected clone", result.products[0].doc)
    wrong = [(p.label, p.doc) for p in result.screen_against]
    args = {"enzymes": None, "excluded": [], "percent": 1.0, "ladder": "neb-1kb-plus"} | kwargs
    return dp.plan(expected, wrong, **args)


def test_planner_tells_the_orientation_and_the_empty_vector_apart():
    out = plan_for_cloning()
    assert out["reason"] is None
    assert 0 < len(out["results"]) <= dp.TOP
    best = out["results"][0]
    assert [lane["name"] for lane in best["lanes"]] == [
        "Expected clone", "insert in the other orientation", "parent vector"]
    assert best["min_separation"] >= dp.CLEAR_FACTOR * gel.RESOLUTION
    assert best["sentence"].startswith("Expected clone: ")
    # A single enzyme cannot see the orientation of a fragment it cuts out symmetrically:
    # nothing listed is BamHI alone.
    assert ["BamHI"] not in [r["enzymes"] for r in out["results"]]


def test_sizes_agree_with_the_gel_simulator():
    """What the planner lists is what Simulate Agarose Gel draws (pydna's digest)."""
    result = cloning()
    docs = [("Expected clone", result.products[0].doc)] + [
        (p.label, p.doc) for p in result.screen_against]
    out = plan_for_cloning()
    for res in out["results"]:
        lanes = gel.simulate([(n, d, res["enzymes"]) for n, d in docs], "neb-1kb-plus", 1.0)
        assert [lane["sizes"] for lane in lanes["lanes"]] == [
            lane["sizes"] for lane in res["lanes"]], res["enzymes"]


def test_ranking_is_deterministic():
    assert plan_for_cloning() == plan_for_cloning()


def test_only_my_enzymes_are_tried_and_methylation_blocked_ones_never():
    mine = ["EcoRI", "HindIII", "PstI", "XbaI", "ClaI", "NcoI", "SalI", "NdeI"]
    out = plan_for_cloning(enzymes=mine, excluded=["ClaI"])
    assert out["enzymes_tried"] == 7
    used = {e for r in out["results"]
            for combo in [r["enzymes"], *r["same_pattern"]] for e in combo}
    assert used <= set(mine) - {"ClaI"}


def test_same_bands_are_one_result():
    out = plan_for_cloning()
    patterns = [tuple(tuple(lane["sizes"]) for lane in r["lanes"]) for r in out["results"]]
    assert len(patterns) == len(set(patterns))


def test_a_point_change_cannot_be_seen_and_the_dialog_says_why():
    changed = PBR322[:1000] + ("A" if PBR322[1000] != "A" else "C") + PBR322[1001:]
    out = dp.plan(("Expected clone", {"sequence": changed, "circular": True}),
                  [("parent", {"sequence": PBR322, "circular": True})],
                  ["EcoRI", "BamHI"], [], 1.0, "neb-1kb-plus")
    assert out["results"] == []
    assert "differs from the expected clone by only 1 base" in out["reason"]
    assert "colony PCR or sequencing" in out["reason"]


def test_the_same_molecule_is_named():
    rotated = PBR322[2000:] + PBR322[:2000]
    out = dp.plan(("Expected clone", {"sequence": PBR322, "circular": True}),
                  [("copy", {"sequence": rotated, "circular": True})],
                  ["EcoRI"], [], 1.0, "neb-1kb-plus")
    assert "is the same molecule as the expected clone" in out["reason"]


def test_small_bands_do_not_count():
    # EcoRI + HindIII cut a 31 bp piece out of pBR322: too small to see.
    m = mol(PBR322, name="pBR322")
    lane = dp._lane(m, ("EcoRI", "HindIII"), 1.0)
    assert [b.size for b in lane.visible] == [4330]
    assert dp.lane_text(lane) == "4.3\u00a0kb (and 31\u00a0bp, too small to see)"


def test_nothing_to_tell_apart_is_an_error():
    with pytest.raises(CloningError, match="at least one wrong clone"):
        dp.plan(("Expected clone", {"sequence": PBR322, "circular": True}), [], None, [], 1.0,
                "neb-1kb-plus")


def test_route(tmp_path, monkeypatch):
    monkeypatch.setenv("VETTORE_HOME", str(tmp_path))
    client = TestClient(create_app("t", None), base_url="http://127.0.0.1")
    result = cloning()
    body = {"expected": {"name": "Expected clone", "sequenceData": result.products[0].doc},
            "alternatives": [{"name": p.label, "sequenceData": p.doc}
                             for p in result.screen_against],
            "enzymes": ["EcoRI", "PstI", "SphI"], "percent": 1.0}
    res = client.post("/api/diagnostic-digest", json=body, headers={TOKEN_HEADER: "t"})
    assert res.status_code == 200
    assert res.json()["results"]
    body["alternatives"] = []
    res = client.post("/api/diagnostic-digest", json=body, headers={TOKEN_HEADER: "t"})
    assert res.status_code == 422
    assert "wrong clone" in res.json()["detail"]
