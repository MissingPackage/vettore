import json
import subprocess
import sys
from pathlib import Path

from vettore.enzymes import MODIFICATION_DEPENDENT, build_database

ROOT = Path(__file__).resolve().parents[1]
BUNDLED = ROOT / "frontend" / "src" / "data" / "enzymes.json"
FIXTURE = ROOT / "frontend" / "src" / "__tests__" / "fixtures" / "cuts.json"


def test_bundled_database_is_up_to_date():
    bundled = json.loads(BUNDLED.read_text())["enzymes"]
    assert bundled == build_database(), (
        "Regenerate: python -m vettore.enzymes > frontend/src/data/enzymes.json")


def test_cross_engine_fixture_is_up_to_date():
    out = subprocess.run([sys.executable, "-I", str(ROOT / "scripts" / "make_enzyme_fixture.py")],
                         capture_output=True, text=True, check=True).stdout
    assert json.loads(out) == json.loads(FIXTURE.read_text()), (
        "Regenerate: python -I scripts/make_enzyme_fixture.py > " + str(FIXTURE.relative_to(ROOT)))


def test_classic_enzymes_have_expected_cut_offsets():
    db = build_database()
    expected = {  # name: (site, top offset, bottom offset)
        "EcoRI": ("gaattc", 1, 5),   # G^AATTC
        "BamHI": ("ggatcc", 1, 5),   # G^GATCC
        "PstI": ("ctgcag", 5, 1),    # CTGCA^G, 3' overhang
        "SmaI": ("cccggg", 3, 3),    # CCC^GGG, blunt
        "XmaI": ("cccggg", 1, 5),    # C^CCGGG
        "NotI": ("gcggccgc", 2, 6),  # GC^GGCCGC
        "BsaI": ("ggtctc", 7, 11),   # GGTCTC(1/5)
        "BsmBI": ("cgtctc", 7, 11),  # CGTCTC(1/5)
    }
    for name, (site, top, bottom) in expected.items():
        e = db[name.lower()]
        assert (e["name"], e["site"], e["topSnipOffset"], e["bottomSnipOffset"]) == (
            name, site, top, bottom)


def test_type_iis_and_golden_gate_flags():
    db = build_database()
    assert db["bsai"]["isType2S"] and db["bsai"]["goldenGate"]
    assert not db["ecori"]["isType2S"] and not db["ecori"]["goldenGate"]


def test_one_enzyme_per_specificity_and_no_modification_dependent_enzymes():
    db = build_database()
    keys = [(e["site"], e["topSnipOffset"], e["bottomSnipOffset"]) for e in db.values()]
    assert len(keys) == len(set(keys))
    assert not {e["name"] for e in db.values()} & MODIFICATION_DEPENDENT
    assert "Esp3I" in db["bsmbi"]["isoschizomers"]
