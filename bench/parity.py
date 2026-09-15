"""Parity harness: replay golden cases through the real server endpoints.

Each record is POSTed to /calculate exactly as the browser would send it, so
the comparison covers the site's full wire behavior — validation order, domain
messages, 500 catch-all — not just the solver. On success the proof equation is
fetched from /get_integral_image with the same fields the frontend sends.

Usage:
    python bench/parity.py [golden.jsonl] [--limit N] [--type pi]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator.server import IMAGE_KEYS, app


def _norm(tex: str) -> str:
    """Whitespace-insensitive LaTeX compare."""
    return "".join(tex.split())


def compare_record(client, rec: dict) -> dict:
    """Replay one golden record against the live endpoints; classify."""
    out = {
        "type": rec["type"],
        "power": rec["power"],
        "comparison": rec["comparison"],
        "rational": rec["rational"],
        "site_success": rec["success"],
        "site_error": rec.get("error"),
    }
    t0 = time.time()
    try:
        resp = client.post("/calculate", data={
            "type": rec["type"], "power": rec["power"],
            "comparison": rec["comparison"], "rational": rec["rational"]})
        body = resp.get_json(silent=True) or {}
    except Exception as exc:  # harness-level failure — record loudly
        out["ours_success"] = False
        out["ours_error"] = f"crash: {type(exc).__name__}: {exc}"
        return out
    out["elapsed_ms"] = round(1000 * (time.time() - t0), 1)
    out["status_match"] = resp.status_code == rec.get("http_status")

    if not body.get("success"):
        out["ours_success"] = False
        out["ours_error"] = body.get("error") or f"http {resp.status_code}"
    else:
        out["ours_success"] = True
        out["ours_parameters"] = body["parameters"]
        if rec["success"]:
            rc = rec.get("raw_calculate") or {}
            if isinstance(rc, str):
                rc = json.loads(rc)
            out["param_match"] = body["parameters"] == rc.get("parameters")
            out["solution_match"] = (
                body.get("equations", {}).get("solution")
                == rc.get("equations", {}).get("solution"))
            if rec.get("equation"):
                img = client.get("/get_integral_image", query_string={
                    **{k: str(body["parameters"][k]) for k in IMAGE_KEYS},
                    "type": rec["type"], "comparison": rec["comparison"],
                    "coef": rec["power"], "rational": rec["rational"]})
                eq = (img.get_json(silent=True) or {}).get("equation")
                if eq is None:
                    out["equation_match"] = None
                    out["equation_error"] = f"image http {img.status_code}"
                else:
                    out["equation_match"] = _norm(eq) == _norm(rec["equation"])
                    if not out["equation_match"]:
                        out["ours_equation"] = eq
    if not out["ours_success"] and not rec["success"]:
        out["error_match"] = out["ours_error"] == rec.get("error")
    return out


def main() -> None:
    """Replay every golden record and print parity metrics."""
    ap = argparse.ArgumentParser()
    ap.add_argument("golden", nargs="?", default="bench/data/golden.jsonl")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--type", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    records = [json.loads(line) for line in open(args.golden) if line.strip()]
    if args.type:
        records = [r for r in records if r["type"] == args.type]
    if args.limit:
        records = records[: args.limit]

    client = app.test_client()
    results = [compare_record(client, r) for r in records]
    if args.out:
        with open(args.out, "w") as f:
            for r in results:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")

    n = len(results)
    both_ok = sum(1 for r in results if r["site_success"] and r.get("ours_success"))
    site_ok_ours_fail = [
        r for r in results if r["site_success"] and not r.get("ours_success")
    ]
    site_fail_ours_ok = [
        r for r in results if not r["site_success"] and r.get("ours_success")
    ]
    exact = sum(1 for r in results if r.get("param_match"))
    sol = sum(1 for r in results if r.get("solution_match"))
    stat = sum(1 for r in results if r.get("status_match"))
    crashes = [r for r in results if str(r.get("ours_error", "")).startswith("crash")]

    eq_match = sum(1 for r in results if r.get("equation_match"))
    eq_total = sum(1 for r in results if r.get("equation_match") is not None)
    both_fail = [r for r in results if not r["site_success"] and not r.get("ours_success")]
    err_match = sum(1 for r in both_fail if r.get("error_match"))
    print(f"total={n} both_ok={both_ok} site_ok_ours_fail={len(site_ok_ours_fail)} "
          f"site_fail_ours_ok={len(site_fail_ours_ok)} param_exact={exact} "
          f"solution_match={sol} status_match={stat} "
          f"eq_match={eq_match}/{eq_total} both_fail={len(both_fail)} "
          f"err_match={err_match} crashes={len(crashes)}")
    err_diff = [r for r in both_fail if r.get("error_match") is False]
    if err_diff:
        print("both fail but error text differs:")
        for r in err_diff[:15]:
            print(f"  {r['type']} {r['power']} {r['comparison']} {r['rational']}: "
                  f"site={r.get('site_error','?')!r} ours={r.get('ours_error')!r}")
    if crashes:
        print("CRASHES (fix first):")
        for r in crashes[:20]:
            print(" ", r["type"], r["power"], r["comparison"], r["rational"], r["ours_error"])
    if site_ok_ours_fail:
        print("site proves, we don't:")
        for r in site_ok_ours_fail[:20]:
            print(" ", r["type"], r["power"], r["comparison"], r["rational"], r.get("ours_error"))


if __name__ == "__main__":
    main()
