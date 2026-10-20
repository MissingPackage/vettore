import pytest
from cloning_helpers import clean_dna, doc, feature, rc
from fastapi import FastAPI
from fastapi.testclient import TestClient

from vettore.cloning.api import router
from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.pcr import PrimerIn, design_primers, run_pcr
from vettore.tm import primer_tm

T = clean_dna(1200, 21)
PLASMID = SequenceDoc(**doc("pTEST", T, True, features=[feature("GFP", 200, 499)]))


def test_pcr_product_includes_tails_and_carries_template_features():
    result = run_pcr(PLASMID, PrimerIn(name="F", sequence="GAATTC" + T[150:172]),
                     PrimerIn(name="R", sequence="GGATCC" + rc(T[550:572])))
    [product] = result.products
    assert product.doc["sequence"] == "GAATTC" + T[150:572] + rc("GGATCC")
    assert product.doc["circular"] is False
    assert [(f["name"], f["start"], f["end"]) for f in product.doc["features"]] == [("GFP", 56, 355)]
    names = {(p["name"], p["strand"]) for p in product.doc["primers"]}
    assert names == {("F", 1), ("R", -1)}
    assert [o.annealing_length for o in result.oligos] == [22, 22]
    assert result.oligos[0].tm == pytest.approx(primer_tm(T[150:172]), abs=0.01)


def test_amplicon_spanning_the_origin_of_a_circular_template():
    result = run_pcr(PLASMID, PrimerIn(name="F", sequence=T[1100:1122]),
                     PrimerIn(name="R", sequence=rc(T[:122])[:22]))
    assert result.products[0].doc["sequence"] == T[1100:] + T[:122]


def test_primer_that_does_not_bind_is_an_error():
    with pytest.raises(CloningError, match="Forward primer “F” does not bind"):
        run_pcr(PLASMID, PrimerIn(name="F", sequence="ACGT" * 6),
                PrimerIn(name="R", sequence=rc(T[550:572])))


def test_primer_that_does_not_bind_says_what_to_check():
    with pytest.raises(CloningError) as err:
        run_pcr(PLASMID, PrimerIn(name="F", sequence=T[150:172]),
                PrimerIn(name="R", sequence="ACGT" * 6))
    assert str(err.value) == ("Reverse primer “R” does not bind “pTEST”. Check that it is "
                              "written 5′→3′ and that the right template is selected.")


def test_a_reverse_primer_copied_from_the_top_strand_gets_its_reverse_complement():
    """The commonest slip: the reverse primer typed as it reads on the map (top strand)."""
    with pytest.raises(CloningError) as err:
        run_pcr(PLASMID, PrimerIn(name="F", sequence=T[150:172]),
                PrimerIn(name="R", sequence=T[550:572]))
    message = str(err.value)
    assert message.startswith("Reverse primer “R” matches the top strand of “pTEST”, so it "
                              "points the same way as the forward primer.")
    assert message.endswith(f"use its reverse complement, 5′-{rc(T[550:572])}-3′.")
    # ... and with that reverse complement the PCR works.
    [product] = run_pcr(PLASMID, PrimerIn(name="F", sequence=T[150:172]),
                        PrimerIn(name="R", sequence=rc(T[550:572]))).products
    assert product.doc["sequence"] == T[150:572]


def test_a_forward_primer_written_on_the_bottom_strand_gets_its_reverse_complement():
    with pytest.raises(CloningError) as err:
        run_pcr(PLASMID, PrimerIn(name="F", sequence=rc(T[150:172])),
                PrimerIn(name="R", sequence=rc(T[550:572])))
    message = str(err.value)
    assert message.startswith("Forward primer “F” matches the bottom strand of “pTEST”")
    assert message.endswith(f"use its reverse complement, 5′-{T[150:172]}-3′.")


def test_swapped_primers_give_the_same_product_with_each_primer_where_it_binds():
    result = run_pcr(PLASMID, PrimerIn(name="F", sequence="GGATCC" + rc(T[550:572])),
                     PrimerIn(name="R", sequence="GAATTC" + T[150:172]))
    [product] = result.products
    assert product.doc["sequence"] == "GAATTC" + T[150:572] + rc("GGATCC")
    primers = {p["name"]: (p["start"], p["end"], p["strand"]) for p in product.doc["primers"]}
    length = len(product.doc["sequence"])
    assert primers == {"R": (0, 27, 1), "F": (length - 28, length - 1, -1)}
    assert result.warnings[0].startswith("“R” binds as the forward primer and “F” as the "
                                         "reverse primer")


def test_primers_facing_away_on_a_linear_template_is_an_error():
    linear = SequenceDoc(**doc("frag", T[:600], False))
    with pytest.raises(CloningError, match="face away"):
        run_pcr(linear, PrimerIn(name="F", sequence=T[400:422]),
                PrimerIn(name="R", sequence=rc(T[100:122])))


def test_tm_difference_is_warned():
    result = run_pcr(PLASMID, PrimerIn(name="F", sequence=T[150:190]),
                     PrimerIn(name="R", sequence=rc(T[550:566])))
    assert any("differ by" in w for w in result.warnings)


def test_designed_primers_reach_the_target_tm_and_carry_the_tails():
    result = design_primers(PLASMID, 200, 700, 62.0, "GAATTC", "AAGCTT")
    forward, reverse = result.oligos
    assert forward.sequence.startswith("GAATTC") and reverse.sequence.startswith("AAGCTT")
    assert primer_tm(forward.sequence[6:]) >= 62.0
    assert primer_tm(reverse.sequence[6:]) >= 62.0
    assert forward.sequence[6:] == T[200 : 200 + forward.annealing_length]
    assert reverse.sequence[6:] == rc(T[701 - reverse.annealing_length : 701])
    assert result.products[0].doc["sequence"] == "GAATTC" + T[200:701] + rc("AAGCTT")


def test_designed_primers_on_a_region_spanning_the_origin():
    result = design_primers(PLASMID, 1100, 199, 60.0)
    assert result.products[0].doc["sequence"] == T[1100:] + T[:200]


def test_region_outside_the_sequence_is_rejected():
    with pytest.raises(CloningError) as err:
        design_primers(PLASMID, 10, 5000)
    assert str(err.value) == ("The region 11..5,001 is not on the template “pTEST”, which runs "
                              "from 1 to 1,200. Choose a region inside it")


def test_http_routes():
    client = TestClient(_app())
    body = {"template": doc("pTEST", T, True),
            "forward": {"name": "F", "sequence": T[150:172]},
            "reverse": {"name": "R", "sequence": rc(T[550:572])}}
    res = client.post("/api/cloning/pcr", json=body)
    assert res.status_code == 200
    assert len(res.json()["products"][0]["doc"]["sequence"]) == 422
    res = client.post("/api/cloning/design-pcr-primers",
                      json={"template": doc("pTEST", T, True), "start": 100, "end": 600})
    assert len(res.json()["oligos"]) == 2


def _app() -> FastAPI:
    app = FastAPI()
    app.include_router(router)
    return app
