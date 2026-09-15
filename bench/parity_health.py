"""Parity harness for the /health sibling app.

Replays every record in bench/data/health-probes.jsonl (verbatim site capture,
see bench/probe_health.py) against our Flask app via test_client and compares
status code + response body byte-for-byte.

Caveats:
- 429 records (site rate limiter) are skipped — we deliberately do not
  reproduce it (a real limiter would break replay); see docs/health-notes.md.
- ``record_id`` is normalized to 0 in both bodies before comparison — the
  site's counter is global/stateful, ours is a local JSON file (HEALTH_DB).
- The last record per tag wins (429 retries leave duplicate tags in the log).

Usage:
    uv run python bench/parity_health.py [--out bench/out/parity_health.jsonl]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from pathlib import Path

# 计数器写到临时文件，别污染仓库（HEALTH_DB 语义见 health.next_record_id）
os.environ["HEALTH_DB"] = os.path.join(tempfile.mkdtemp(), "records.json")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.server import app

DATA = Path(__file__).resolve().parent / "data"
RID_RE = re.compile(r'"record_id":\d+')


def replay(client, rec: dict) -> tuple[int, str]:
    """Replay one probe record; return (status_code, response_body_text)."""
    path, method = rec["path"], rec["method"]
    if method in ("GET", "OPTIONS", "HEAD"):
        resp = client.open(path, method=method)
    elif rec.get("raw") is not None:
        resp = client.open(path, method=method, data=rec["raw"],
                           content_type=rec.get("content_type"))
    else:
        resp = client.open(path, method=method, json=rec.get("json"))
    return resp.status_code, resp.get_data(as_text=True)


def check(client, rec: dict) -> dict:
    """Compare one record; normalize record_id before byte-compare."""
    out = {"tag": rec["tag"], "site_status": rec["http"]}
    status, body = replay(client, rec)
    out["ours_status"] = status
    out["status_match"] = status == rec["http"]
    site_body = RID_RE.sub('"record_id":0', rec["response_text"] or "")
    ours_body = RID_RE.sub('"record_id":0', body)
    out["body_match"] = ours_body == site_body
    if not out["body_match"]:
        out["ours_body"] = ours_body
        out["site_body"] = site_body
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path,
                   default=Path("bench/out/parity_health.jsonl"))
    p.add_argument("--probes", type=Path, default=DATA / "health-probes.jsonl")
    args = p.parse_args()

    latest = {}
    for line in args.probes.open():
        if line.strip():
            r = json.loads(line)
            latest[r["tag"]] = r
    skipped_429 = sum(1 for r in latest.values() if r["http"] == 429)

    client = app.test_client()
    results = [check(client, r) for r in latest.values() if r["http"] != 429]
    bad = [r for r in results if not (r["status_match"] and r["body_match"])]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")

    print(f"records: {len(latest)} unique tags ({skipped_429} still-429 skipped)")
    print(f"compared: {len(results)}; mismatches: {len(bad)}")
    for r in bad[:30]:
        print(f"  FAIL {r['tag']}: site={r['site_status']} ours={r['ours_status']}")
        if not r["body_match"]:
            print(f"    site: {(r['site_body'] or '')[:160]}")
            print(f"    ours: {(r['ours_body'] or '')[:160]}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
