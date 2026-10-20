"""Vettore's own sources never name another sequence-editing product.

The names are stored only as BLAKE2b digests, so this file does not spell them. Every run of one
to three words in a line (case-insensitive; camelCase, spaces, hyphens and underscores all
count as word breaks, and the words are joined without separators) is hashed and compared.

The one allowed occurrence is the identifier Biopython gives the binary .dna format
(``DNA_FILE_FORMAT`` in ``src/vettore/fileio.py``), which that library fixes; its digest is
taken from that constant. The vendored editor library (``frontend/node_modules``) and the built
bundle (``frontend/dist`` and its packaging copy ``src/vettore/static``) are third-party code and
are not scanned.
"""

import hashlib
import re
from pathlib import Path

from vettore import fileio

ROOT = Path(__file__).resolve().parents[1]

# Commercial sequence editors, their suites and vendors, normalised as lower-case words joined
# without separators.
DIGESTS = frozenset({
    "26570cf85c0602eb", "32b348a5eadbe04a", "3cfd85f316ad6385", "42b392b65f1e9699",
    "4a8f8890a6b9f833", "5a2281774cb31d59", "8962055b7dfe9f3d", "ae51a1115c5fb498",
    "aef4b51b8e02acf7", "afabb90ef1b3f654", "bd1196ee3066722d", "c7e2e61a9032ab9b",
    "cf15a7b389a5ee1d", "f2b0ccda01fbe167", "f73039e6c380f264",
})

SCANNED = ["src", "frontend/src", "frontend/e2e", "docs", "README.md"]
TEXT_SUFFIXES = {".py", ".ts", ".tsx", ".js", ".mjs", ".css", ".md", ".html", ".json", ".txt",
                 ".svg", ".toml", ".yml", ".yaml"}
ALLOWED = ("src/vettore/fileio.py", "DNA_FILE_FORMAT = ")


def digest(text: str) -> str:
    return hashlib.blake2b(text.encode(), digest_size=8).hexdigest()


FORBIDDEN = DIGESTS | {digest(fileio.DNA_FILE_FORMAT.lower())}

_CAMEL = re.compile(r"(?<=[a-z0-9])(?=[A-Z])")
_WORD = re.compile(r"[A-Za-z0-9]+")


def keys(line: str) -> set[str]:
    """Every run of one to three words of a line, lower case, joined without separators; a
    camelCase identifier counts both whole and split into its words."""
    tokens = _WORD.findall(line)
    words = [w.lower() for token in tokens for w in _CAMEL.sub(" ", token).split()]
    found = {token.lower() for token in tokens}
    for size in (1, 2, 3):
        for i in range(len(words) - size + 1):
            found.add("".join(words[i:i + size]))
    return found


def sources() -> list[Path]:
    files = [ROOT / "tests" / p.name for p in (ROOT / "tests").glob("*.py")]
    for entry in SCANNED:
        path = ROOT / entry
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files += [p for p in path.rglob("*")
                      if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES
                      and "node_modules" not in p.parts
                      # the packaging copy of the built bundle (third-party code, git-ignored)
                      and not p.is_relative_to(ROOT / "src" / "vettore" / "static")]
    return sorted(set(files))


def hits(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return []
    relative = path.relative_to(ROOT).as_posix()
    out = []
    for number, line in enumerate(text.splitlines(), 1):
        if (relative, line.lstrip()[:len(ALLOWED[1])]) == ALLOWED:
            continue
        if any(digest(key) in FORBIDDEN for key in keys(line)):
            out.append(f"{relative}:{number}")
    return out


def test_no_source_file_names_another_sequence_editor():
    found = [hit for path in sources() for hit in hits(path)]
    assert not found, ("These lines name another sequence-editing product; describe the "
                       "format or the behaviour instead: " + ", ".join(found))


def test_the_scan_sees_the_allowed_identifier_and_the_word_breaks():
    """Guards the guard: the scan reads fileio.py, and spelling variants are caught."""
    assert ROOT / "src" / "vettore" / "fileio.py" in sources()
    name = fileio.DNA_FILE_FORMAT.lower()
    half = len(name) // 2
    variants = [name, name.upper(), f"{name[:half]} {name[half:]}",
                f"{name[:half].capitalize()}{name[half:].capitalize()}",
                f"read_{name}_file", f"{name[:half]}-{name[half:]}"]
    for variant in variants:
        assert any(digest(key) in FORBIDDEN for key in keys(f"x {variant} y")), variant
    assert hits(ROOT / "src" / "vettore" / "fileio.py") == []
