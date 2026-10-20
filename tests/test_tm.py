import pytest

from vettore.tm import PRESETS, all_presets, primer_tm

SEQ = "ACGTAGCTAGCTAGGCTATCGGATCCTA"


def test_q5_preset_matches_the_benchmarked_value():
    # Benchmarked against NEB's Tm API when the preset was calibrated: 70.53 °C for this 28-mer.
    assert primer_tm(SEQ, "neb_q5") == pytest.approx(70.53, abs=0.01)


def test_q5_preset_tracks_neb_calculator_example():
    # NEB's Tm API returns 66 °C (Q5, 0.5 µM) for this primer.
    assert primer_tm("AGCGGATAACAATTTCACACAGGA", "neb_q5") == pytest.approx(66, abs=0.5)


def test_default_preset_is_the_standard_conditions():
    assert primer_tm(SEQ) == primer_tm(SEQ, "standard")
    assert "250 nM" in PRESETS["standard"].description


def test_tm_is_case_insensitive_and_accepts_u():
    assert primer_tm("acgu" + SEQ[4:].lower()) == pytest.approx(primer_tm(SEQ))


@pytest.mark.parametrize("bad", ["", "ACGTN", "A" * 101])
def test_tm_rejects_invalid_input(bad):
    with pytest.raises(ValueError):
        primer_tm(bad)


def test_all_presets_report_label_description_and_value():
    rows = all_presets(SEQ)
    assert [r["key"] for r in rows] == ["standard", "neb_q5", "neb_taq"]
    assert all(r["description"] and r["tm"] > 50 for r in rows)
