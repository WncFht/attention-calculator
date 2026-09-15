"""Parity harness for POST /decompose_inequality.

Two golden sets with different capture formats:

- ``bench/data/combo.jsonl`` — records carry ``raw`` (verbatim site response
  body) + ``http_status``; compared byte-for-byte like parity.py does for
  /calculate.
- ``bench/data/decompose.jsonl`` — records carry ``response`` (parsed JSON
  object) or ``{'success': False, 'transport_error': 'HTTP 400'}`` for site
  failures; compared as parsed JSON since the raw body was not captured.
  Transport-error records expect HTTP 400 with any ``{"error": ...}`` body.

Usage:
    python bench/parity_decompose.py [--out bench/out/parity_decompose.jsonl]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.server import app

DATA = Path(__file__).resolve().parent / "data"


def check_combo(client, rec: dict) -> dict:
    """Replay one combo record; byte-compare the response body."""
    out = {"problem": rec["problem"], "site_status": rec["http_status"]}
    resp = client.post("/decompose_inequality", data={"problem": rec["problem"]})
    out["ours_status"] = resp.status_code
    out["status_match"] = resp.status_code == rec["http_status"]
    out["body_match"] = resp.get_data(as_text=True) == rec["raw"]
    if not out["body_match"]:
        out["ours_body"] = resp.get_data(as_text=True)
        out["site_body"] = rec["raw"]
    return out


def check_decompose(client, rec: dict) -> dict:
    """Replay one decompose record; compare parsed JSON or expect a 400."""
    out = {"problem": rec["problem"], "note": rec.get("note")}
    site = rec.get("response") or {}
    resp = client.post("/decompose_inequality", data={"problem": rec["problem"]})
    out["ours_status"] = resp.status_code
    if "transport_error" in site:
        # 站端 400；采集脚本没存 body，只要求状态码 + error 键
        body = resp.get_json(silent=True) or {}
        out["status_match"] = resp.status_code == 400
        out["body_match"] = isinstance(body.get("error"), str)
        if not out["body_match"]:
            out["ours_body"] = resp.get_data(as_text=True)
        return out
    # records carry no http status; {"error":...} bodies are the site's 400s
    want = 400 if "error" in site else 200
    out["status_match"] = resp.status_code == want
    out["body_match"] = resp.get_json(silent=True) == site
    if not out["body_match"]:
        out["ours_body"] = resp.get_data(as_text=True)
        out["site_body"] = json.dumps(site, ensure_ascii=False)[:2000]
    return out


def main() -> None:
    """Run both sets and print parity metrics."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    client = app.test_client()
    results = []
    for name, checker in (("combo", check_combo), ("decompose", check_decompose)):
        with open(DATA / f"{name}.jsonl") as fh:
            for rec in (json.loads(line) for line in fh if line.strip()):
                results.append({"set": name, **checker(client, rec)})

    if args.out:
        with open(args.out, "w") as f:
            for r in results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")

    for name in ("combo", "decompose"):
        rs = [r for r in results if r["set"] == name]
        stat = sum(1 for r in rs if r["status_match"])
        body = sum(1 for r in rs if r["body_match"])
        print(f"{name}: total={len(rs)} status_match={stat} body_match={body}")
        for r in rs:
            if not (r["status_match"] and r["body_match"]):
                print(f"  DIFF {r['problem']}: site={r.get('site_status','?')} "
                      f"ours={r['ours_status']}")
                print(f"    ours: {r.get('ours_body','')[:300]}")
                print(f"    site: {r.get('site_body','')[:300]}")


if __name__ == "__main__":
    main()
