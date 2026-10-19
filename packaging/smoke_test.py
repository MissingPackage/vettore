"""Smoke-test a packaged Vettore build: start it headless, open a .dna file through the API.

Usage: python packaging/smoke_test.py <path-to-executable>
"""

import base64
import json
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

PORT, TOKEN = 8799, "smoke-token"
DATA = Path(__file__).resolve().parents[1] / "tests" / "data" / "pFA-KanMX4.dna"


def call(path: str, body: dict | None = None) -> dict:
    req = urllib.request.Request(
        f"http://127.0.0.1:{PORT}/api/{path}",
        data=None if body is None else json.dumps(body).encode(),
        headers={"X-Vettore-Token": TOKEN, "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=10) as res:
        return json.loads(res.read())


def main() -> int:
    exe = sys.argv[1]
    proc = subprocess.Popen([exe, "--no-window", "--port", str(PORT), "--token", TOKEN])
    try:
        deadline = time.monotonic() + 90
        while True:
            try:
                assert call("health") == {"status": "ok"}
                break
            except OSError:
                if proc.poll() is not None:
                    print(f"executable exited early with code {proc.returncode}")
                    return 1
                if time.monotonic() > deadline:
                    print("server did not start within 90 s")
                    return 1
                time.sleep(1)
        result = call("open-content", {
            "filename": DATA.name,
            "content_b64": base64.b64encode(DATA.read_bytes()).decode(),
        })
        [doc] = result["documents"]
        assert len(doc["sequence"]) == 3941, len(doc["sequence"])
        assert len(doc["features"]) == 9, len(doc["features"])
        with urllib.request.urlopen(f"http://127.0.0.1:{PORT}/index.html", timeout=10) as res:
            assert b'id="root"' in res.read()
        print("smoke test passed: server up, .dna parsed (3941 bp, 9 features), UI served")
        return 0
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    sys.exit(main())
