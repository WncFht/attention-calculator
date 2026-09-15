"""Parity harness: run our solver over golden cases, compare with the site.

Metrics per bench/README.md: success agreement, parameter exact-match, and
(when integrand.py is available) validity of our produced identity.

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

from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.solve import prove


def compare_record(rec: dict) -> dict:
    """Run our solver on one golden record; classify the outcome."""
    out = {
        "type": rec["type"],
        "power": rec["power"],
        "comparison": rec["comparison"],
        "rational": rec["rational"],
        "site_success": rec["success"],
    }
    t0 = time.time()
    try:
        res = prove(rec["type"], rec["power"], rec["comparison"], rec["rational"])
        out["ours_success"] = True
        out["ours_parameters"] = res["parameters"]
        out["elapsed_ms"] = round(1000 * (time.time() - t0), 1)
        if rec["success"]:
            ours = res["parameters"]
            site = rec["parameters"]
            out["param_match"] = all(
                str(ours.get(k)) == str(site.get(k))
                for k in ("m", "n", "a_val", "b_val", "c_val", "u_val")
            )
    except WrongDirection:
        out["ours_success"] = False
        out["ours_error"] = "wrong_direction"
    except NoSolution:
        out["ours_success"] = False
        out["ours_error"] = "no_solution"
    except Exception as exc:  # solver bug — record loudly, don't hide
        out["ours_success"] = False
        out["ours_error"] = f"crash: {type(exc).__name__}: {exc}"
    return out


def main() -> None:
    """Compare our solver against every golden record and print metrics."""
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

    results = [compare_record(r) for r in records]
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
    crashes = [r for r in results if str(r.get("ours_error", "")).startswith("crash")]

    print(f"total={n} both_ok={both_ok} site_ok_ours_fail={len(site_ok_ours_fail)} "
          f"site_fail_ours_ok={len(site_fail_ours_ok)} param_exact={exact} crashes={len(crashes)}")
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
