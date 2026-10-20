"""Build the methylation-sensitivity table from REBASE.

REBASE publishes, per enzyme, whether methylation inside its recognition site blocks it. That
is the difference between a useful warning and a nagging one: BamHI's GGATCC always overlaps a
Dam site and is cut anyway, while ClaI's ATCGAT is blocked when a Dam site overlaps it.

Run (needs network):
    python -I scripts/build_methylation_table.py > frontend/src/data/methylation.json

Source: REBASE, Roberts RJ, Vincze T, Posfai J, Macelis D, Nucleic Acids Res.
<http://rebase.neb.com>. REBASE may be redistributed free of charge with attribution.
"""

from __future__ import annotations

import json
import sys
import urllib.request

URL = "ftp://ftp.neb.com/pub/rebase/msbuffmore.txt"
#: The methylations a plasmid picks up in ordinary cloning strains.
KINDS = ("Dam", "Dcm", "CpG")
#: Statuses REBASE uses. "cut" means methylation does not stop the enzyme; the interface
#: reports the others, and says so when REBASE records nothing for an enzyme at all.
KNOWN = {"blocked", "some blocked", "impaired", "variable", "cut"}


def main() -> None:
    with urllib.request.urlopen(URL, timeout=120) as response:
        text = response.read().decode("latin-1")

    table: dict[str, dict[str, str]] = {}
    for line in text.splitlines():
        parts = line.split(";")
        if len(parts) < 5 or not parts[0] or " " in parts[0]:
            continue
        name, sensitivity = parts[0].strip(), parts[4].strip()
        if not sensitivity:
            continue
        statuses: dict[str, str] = {}
        for item in sensitivity.split(","):
            kind, _, status = item.partition(":")
            kind, status = kind.strip(), status.strip()
            if kind in KINDS and status in KNOWN:
                statuses[kind] = status
        if statuses:
            table[name.lower()] = statuses

    if len(table) < 200:
        raise SystemExit(f"only {len(table)} enzymes parsed; the REBASE format may have changed")
    json.dump({"source": f"REBASE {URL} (http://rebase.neb.com)",
               "statuses": dict(sorted(table.items()))}, sys.stdout, indent=1)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
