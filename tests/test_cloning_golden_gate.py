import pytest
from cloning_helpers import clean_dna, doc, same_circular

from vettore.cloning.base import CloningError, SequenceDoc
from vettore.cloning.golden_gate import golden_gate

# BsaI: GGTCTC(N1) cuts leaving a 4 nt 5' overhang starting one base after the site.
# A part is  <filler> GGTCTC N <overhang> ...core... <overhang> N GAGACC <filler>
BSA_F, BSA_R = "GGTCTC", "GAGACC"


def part(core: str, left: str, right: str, seed: int) -> str:
    pad = clean_dna(20, seed)
    return pad + BSA_F + "A" + left + core + right + "T" + BSA_R + pad


def vector(backbone: str, dropout: str, first: str, last: str) -> str:
    """A Golden Gate destination vector, written from the backbone outwards: cutting releases
    the dropout and leaves a backbone that starts with `last` and ends with `first`."""
    return last + backbone + first + "T" + BSA_R + dropout + BSA_F + "A"


A = clean_dna(300, 71)
B = clean_dna(250, 72)
C = clean_dna(200, 73)
DROPOUT = clean_dna(150, 74)
BACKBONE = clean_dna(1200, 70)

VECTOR = SequenceDoc(**doc("pGG-dest", vector(BACKBONE, DROPOUT, "AATG", "CGCT"), True))
PART1 = SequenceDoc(**doc("promoter", part(A, "AATG", "GCTT", 75), False))
PART2 = SequenceDoc(**doc("CDS", part(B, "GCTT", "TACA", 76), False))
PART3 = SequenceDoc(**doc("terminator", part(C, "TACA", "CGCT", 77), False))


def expected_product() -> str:
    """Each junction overhang appears once in the circle: CGCT closes it."""
    return "CGCT" + BACKBONE + "AATG" + A + "GCTT" + B + "TACA" + C


def test_three_part_assembly_gives_the_expected_plasmid():
    result = golden_gate(VECTOR, [PART1, PART2, PART3])
    [product] = result.products
    assert product.doc["circular"] is True
    assert same_circular(product.doc["sequence"], expected_product())
    assert "junction overhangs AATG, GCTT, TACA, CGCT" in product.details
    assert product.parents == ["pGG-dest", "promoter", "CDS", "terminator"]


def test_dropout_fragment_is_reported():
    result = golden_gate(VECTOR, [PART1, PART2, PART3])
    assert any("dropout fragment" in w for w in result.warnings)


def test_duplicate_overhangs_are_warned():
    same = SequenceDoc(**doc("CDS2", part(B, "GCTT", "GCTT", 78), False))
    part3 = SequenceDoc(**doc("terminator", part(C, "GCTT", "CGCT", 79), False))
    result = golden_gate(VECTOR, [PART1, same, part3])
    assert any("Duplicate junction overhangs" in w for w in result.warnings)


def test_palindromic_overhang_is_warned():
    first = SequenceDoc(**doc("promoter", part(A, "AATG", "GATC", 80), False))
    second = SequenceDoc(**doc("CDS", part(B, "GATC", "CGCT", 81), False))
    result = golden_gate(VECTOR, [first, second])
    assert any("Palindromic overhangs" in w for w in result.warnings)


def test_mismatched_junction_is_an_error_naming_the_parts():
    wrong = SequenceDoc(**doc("CDS", part(B, "TTTT", "TACA", 82), False))
    with pytest.raises(CloningError) as err:
        golden_gate(VECTOR, [PART1, wrong, PART3])
    assert "promoter and CDS" in str(err.value)
    assert "GCTT" in str(err.value) and "TTTT" in str(err.value)


def test_internal_site_needs_domestication():
    core = A[:100] + BSA_F + "ACGTACGT" + A[100:]
    dirty = SequenceDoc(**doc("promoter", part(core, "AATG", "GCTT", 83), False))
    with pytest.raises(CloningError, match="domesticated"):
        golden_gate(VECTOR, [dirty, PART2, PART3])


def test_part_without_sites_is_an_error():
    plain = SequenceDoc(**doc("plain", A, False))
    with pytest.raises(CloningError, match="does not cut the part"):
        golden_gate(VECTOR, [plain])


def test_bsmbi_variant_works():
    pad = clean_dna(20, 91)
    backbone = clean_dna(900, 90)
    esp_part = pad + "CGTCTC" + "A" + "AATG" + A + "CGCT" + "T" + "GAGACG" + pad
    vec = "CGCT" + backbone + "AATG" + "T" + "GAGACG" + DROPOUT + "CGTCTC" + "A"
    result = golden_gate(SequenceDoc(**doc("pGG2", vec, True)),
                         [SequenceDoc(**doc("part", esp_part, False))], enzyme_name="BsmBI")
    [product] = result.products
    assert product.doc["circular"] is True
    assert same_circular(product.doc["sequence"], "CGCT" + backbone + "AATG" + A)
