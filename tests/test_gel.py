import pytest
from fastapi.testclient import TestClient

from cloning_helpers import clean_dna, doc
from vettore.gel import (DEFAULT_LADDER, LADDERS, lane_bands, ladder_bands, migration,
                         resolving_range)
from vettore.server import TOKEN_HEADER, create_app

TOKEN = "t0ken"


def test_migration_is_monotonic_in_size():
    sizes = [100, 250, 500, 1000, 3000, 10000]
    rfs = [migration(s, 1.0) for s in sizes]
    assert rfs == sorted(rfs, reverse=True)
    assert all(0 < rf < 1 for rf in rfs)


def test_bands_outside_the_resolving_range_stay_on_the_gel():
    assert 0 < migration(20, 0.7) < 1
    assert 0 < migration(100000, 2.0) < 1


def test_a_higher_agarose_percentage_resolves_smaller_fragments():
    low, high = resolving_range(2.0)
    assert low < resolving_range(0.7)[0] and high < resolving_range(0.7)[1]
    # A 200 bp fragment runs further from the well on a 0.7 % gel than on a 2 % gel.
    assert migration(200, 0.7) > migration(200, 2.0)


def test_lane_intensity_follows_fragment_mass():
    bands = lane_bands([4000, 500], 1.0)
    big = next(b for b in bands if b.size == 4000)
    small = next(b for b in bands if b.size == 500)
    assert big.intensity == 1.0
    assert small.intensity < big.intensity


def test_unresolvable_fragments_merge_into_one_band():
    bands = lane_bands([1000, 1010, 300], 1.0)
    assert len(bands) == 2
    merged = next(b for b in bands if b.rf < max(x.rf for x in bands))
    assert "1,000 bp" in merged.label and "1,010 bp" in merged.label


def test_supercoiled_plasmid_runs_ahead_of_its_linear_length():
    [linear] = lane_bands([5000], 1.0)
    [supercoiled] = lane_bands([5000], 1.0, supercoiled=True)
    assert supercoiled.rf > linear.rf
    assert supercoiled.apparent_size < 5000
    assert "supercoiled" in supercoiled.label


def test_every_ladder_has_ordered_bands_and_reference_bands_present():
    for key, entry in LADDERS.items():
        assert entry["bands"] == sorted(entry["bands"], reverse=True), key
        assert set(entry["reference"]) <= set(entry["bands"]), key
        assert len(ladder_bands(key, 1.0)) == len(entry["bands"])


def test_gel_endpoint_digests_each_lane():
    sequence = clean_dna(1500, 61) + "GAATTC" + clean_dna(900, 62) + "GGATCC" + clean_dna(600, 63)
    client = TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")
    body = {"lanes": [{"name": "pTest", "sequenceData": doc("pTest", sequence, True),
                       "enzymes": ["EcoRI", "BamHI"]},
                      {"name": "pTest", "sequenceData": doc("pTest", sequence, True)}],
            "ladder": DEFAULT_LADDER, "percent": 1.0}
    res = client.post("/api/gel", headers={TOKEN_HEADER: TOKEN}, json=body)
    assert res.status_code == 200
    data = res.json()
    assert data["lanes"][0]["sizes"] == [2106, 906]
    assert data["lanes"][0]["name"] == "pTest / EcoRI + BamHI"
    assert data["lanes"][1]["bands"][0]["label"].endswith("(supercoiled)")
    assert len(data["ladder"]) == len(LADDERS[DEFAULT_LADDER]["bands"])


def test_ladders_endpoint():
    client = TestClient(create_app(TOKEN, None), base_url="http://127.0.0.1")
    data = client.get("/api/gel/ladders", headers={TOKEN_HEADER: TOKEN}).json()
    assert data["default_ladder"] == DEFAULT_LADDER
    assert {entry["id"] for entry in data["ladders"]} == set(LADDERS)
