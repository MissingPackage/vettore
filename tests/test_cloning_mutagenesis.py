import pytest
from cloning_helpers import clean_dna, doc, feature, rc

from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.mutagenesis import mutagenesis
from vettore.tm import primer_tm

T = clean_dna(3000, 31)
PLASMID = SequenceDoc(**doc("pAAV", T, True,
                            features=[feature("GFP", 100, 799), feature("ori", 2000, 2499)]))


def test_substitution_edits_the_plasmid_and_keeps_features():
    result = mutagenesis(PLASMID, "substitution", 500, 502, "GCA")
    [product] = result.products
    assert product.doc["sequence"] == T[:500] + "GCA" + T[503:]
    assert product.doc["circular"] is True
    assert [(f["name"], f["start"], f["end"]) for f in product.doc["features"]] == [
        ("GFP", 100, 799), ("ori", 2000, 2499)]
    assert product.details == f"Substituted 3 bp at 501..503 ({T[500:503]} → GCA)"


def test_substitution_primers_are_back_to_back_with_the_change_on_the_forward_tail():
    result = mutagenesis(PLASMID, "substitution", 500, 502, "GCA")
    forward, reverse = result.oligos
    assert forward.sequence.startswith("GCA")
    assert forward.sequence[3:] == T[503 : 503 + forward.annealing_length]
    assert reverse.sequence == rc(T[500 - reverse.annealing_length : 500])
    assert primer_tm(forward.sequence[3:]) >= 60.0
    assert primer_tm(reverse.sequence) >= 60.0


def test_deletion_removes_the_region_without_tails():
    result = mutagenesis(PLASMID, "deletion", 1000, 1099)
    assert result.products[0].doc["sequence"] == T[:1000] + T[1100:]
    forward, reverse = result.oligos
    assert forward.sequence == T[1100 : 1100 + forward.annealing_length]
    assert reverse.sequence == rc(T[1000 - reverse.annealing_length : 1000])


def test_short_insertion_sits_on_the_forward_primer():
    result = mutagenesis(PLASMID, "insertion", 1500, 1500, "GGATCC")
    assert result.products[0].doc["sequence"] == T[:1501] + "GGATCC" + T[1501:]
    assert result.oligos[0].sequence.startswith("GGATCC")
    assert result.oligos[1].sequence == rc(T[1501 - result.oligos[1].annealing_length : 1501])


def test_long_insertion_is_split_between_both_primers():
    flag = "ATGGATTACAAGGATGACGACGATAAG"  # 27 bp FLAG-like tag
    result = mutagenesis(PLASMID, "insertion", 1500, 1500, flag)
    assert result.products[0].doc["sequence"] == T[:1501] + flag + T[1501:]
    forward, reverse = result.oligos
    assert forward.sequence.startswith(flag[13:])
    assert reverse.sequence.startswith(rc(flag[:13]))
    assert any("split between both primers" in w for w in result.warnings)


def test_kld_reminder_is_always_given():
    result = mutagenesis(PLASMID, "deletion", 10, 20)
    assert any("KLD" in w for w in result.warnings)


def test_linear_input_is_rejected():
    linear = SequenceDoc(**doc("frag", T[:500], False))
    with pytest.raises(CloningError, match="linear"):
        mutagenesis(linear, "deletion", 10, 20)


def test_substitution_without_new_bases_is_rejected():
    with pytest.raises(CloningError, match="Enter the new bases"):
        mutagenesis(PLASMID, "substitution", 10, 20)
