"""Parity harness for the /convex sibling app.

Replays every record in bench/data/convex-probes.jsonl (verbatim site
capture, see bench/probe_convex.py) against our Flask app via test_client
and compares status code + response body byte-for-byte.

Record schema: ``form`` is the intended POST form dict;
``encode`` records how the site request was sent (form/json/raw/multipart)
— we always replay urlencoded since the server only reads request.form and
the non-form records all carry empty form dicts anyway. ``raw`` is the
verbatim response body. Last record per tag wins.

Usage:
    uv run python bench/parity_convex.py [--out bench/out/parity_convex.jsonl]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.server import app

DATA = Path(__file__).resolve().parent / "data"


def replay(client, rec: dict) -> tuple[int, str]:
    """Replay one probe record; return (status_code, response_body_text)."""
    path, method = rec["endpoint"], rec["method"]
    form = rec.get("form") or {}
    if method in ("GET", "OPTIONS", "HEAD"):
        resp = client.open(path, method=method)
    else:
        resp = client.open(path, method=method, data=form)
    return resp.status_code, resp.get_data(as_text=True)


def check(client, rec: dict) -> dict:
    """Compare one record byte-for-byte (status + raw body)."""
    out = {"tag": rec["tag"], "site_status": int(rec["http"])}
    status, body = replay(client, rec)
    out["ours_status"] = status
    out["status_match"] = status == int(rec["http"])
    site_body = rec["raw"] or ""
    out["body_match"] = body == site_body
    if not out["body_match"]:
        out["ours_body"] = body
        out["site_body"] = site_body
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path,
                   default=Path("bench/out/parity_convex.jsonl"))
    p.add_argument("--probes", type=Path, default=DATA / "convex-probes.jsonl")
    args = p.parse_args()

    latest = {}
    for line in args.probes.open():
        if line.strip():
            r = json.loads(line)
            latest[r["tag"]] = r

    client = app.test_client()
    results = [check(client, r) for r in latest.values()]
    bad = [r for r in results if not (r["status_match"] and r["body_match"])]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    print(f"records: {len(latest)} unique tags")
    print(f"compared: {len(results)}; mismatches: {len(bad)}")
    for r in bad[:30]:
        print(f"  FAIL {r['tag']}: site={r['site_status']} ours={r['ours_status']}")
        if not r["body_match"]:
            print(f"    site: {(r['site_body'] or '')[:200]}")
            print(f"    ours: {(r['ours_body'] or '')[:200]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
