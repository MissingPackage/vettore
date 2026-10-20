"""Windows reads and writes text in its code page (cp1252) unless told otherwise, so every
text file the engine opens names its encoding. Binary I/O (read_bytes, mode "rb") needs none."""

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = sorted((ROOT / "src" / "vettore").rglob("*.py")) + [ROOT / "packaging" / "smoke_test.py"]
TEXT_CALLS = {"read_text", "write_text", "open", "FileHandler", "TextIOWrapper"}
SUBPROCESS = {"run", "check_output", "Popen"}


def _keywords(call: ast.Call) -> set[str]:
    return {k.arg for k in call.keywords if k.arg}


def _binary_mode(call: ast.Call) -> bool:
    modes = [a for a in call.args[1:2]] + [k.value for k in call.keywords if k.arg == "mode"]
    return any(isinstance(m, ast.Constant) and "b" in str(m.value) for m in modes)


def problems(path: Path) -> list[str]:
    out = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        shown = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.name
        where = f"{shown}:{node.lineno}"
        keywords = _keywords(node)
        if name in TEXT_CALLS and name != "open" or (name == "open" and isinstance(func, ast.Name)):
            # read_text("utf-8") passes the encoding positionally
            positional = name == "read_text" and node.args
            if not ({"encoding"} & keywords or positional or _binary_mode(node)):
                out.append(f"{where} {name}() without an encoding")
        if name in SUBPROCESS and ({"text", "universal_newlines"} & keywords
                                   and "encoding" not in keywords):
            out.append(f"{where} {name}(text=True) without an encoding")
    return out


def test_every_text_file_is_opened_with_an_explicit_encoding():
    found = [problem for path in SOURCES for problem in problems(path)]
    assert not found, found


def test_the_check_catches_a_missing_encoding(tmp_path):
    sample = tmp_path / "sample.py"
    sample.write_text("from pathlib import Path\nPath('x').read_text()\nopen('y', 'w')\n"
                      "open('z', 'rb')\nPath('x').read_text('utf-8')\n", encoding="utf-8")
    found = problems(sample)
    assert [p.split(" ", 1)[1] for p in found] == ["read_text() without an encoding",
                                                   "open() without an encoding"]
