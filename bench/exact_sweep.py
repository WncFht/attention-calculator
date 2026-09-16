#!/usr/bin/env python3
"""Sweep exact_check.verify over a golden-style jsonl corpus.

Every exact-mode kernel session rewrote this loop by hand (~20 instances
across transcripts): read records, keep successes, run the family checker,
tally identity_ok x nonneg.  This is the canned version.

    .venv/bin/python bench/exact_sweep.py bench/data/golden.jsonl
    .venv/bin/python bench/exact_sweep.py bench/data/golden.jsonl --types sin_q,cos_q
    .venv/bin/python bench/exact_sweep.py golden.jsonl --list-fail
    .venv/bin/python bench/exact_sweep.py golden.jsonl --exact-mode   # re-prove + verify smoke

Exit code 1 if any checked record fails (identity not ok or nonneg false),
so it can gate.
"""

import argparse
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from attention_calculator import solve
from attention_calculator.exact_check import verify, verify_response


def main() -> int:
    """CLI entry: load records, run the checker, print tallies."""
    ap = argparse.ArgumentParser(description="sweep exact_check.verify over a jsonl corpus")
    ap.add_argument("corpus", nargs="?", default="bench/data/golden.jsonl")
    ap.add_argument("--types", help="comma-separated type filter")
    ap.add_argument(
        "--include-fail",
        action="store_true",
        help="also check records with success=false (parameters may be absent)",
    )
    ap.add_argument(
        "--exact-mode",
        action="store_true",
        help="run solve.prove(exact=True) and verify the emitted params, not the stored ones",
    )
    ap.add_argument("--list-fail", action="store_true", help="print each failing record's inputs")
    args = ap.parse_args()

    want = set(args.types.split(",")) if args.types else None
    tally = Counter()
    by_type = defaultdict(Counter)
    fails = []
    skipped = Counter()

    for line in Path(args.corpus).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        kind = r["type"]
        if want and kind not in want:
            continue
        if not args.exact_mode and not r.get("success") and not args.include_fail:
            skipped["not-success"] += 1
            continue
        try:
            power, comp, bound = Fraction(r["power"]), r["comparison"], Fraction(r["rational"])
        except (KeyError, ValueError, ZeroDivisionError):
            skipped["unparseable-input"] += 1
            continue
        if args.exact_mode:
            try:
                resp = solve.prove(kind, r["power"], comp, r["rational"], exact=True)
            except Exception as e:
                key = f"prove:{type(e).__name__}"
                tally[key] += 1
                by_type[kind][key] += 1
                continue
            if not resp.get("parameters"):
                key = "prover:" + resp.get("prover", "?")
                tally[key] += 1
                by_type[kind][key] += 1
                continue
            res = verify_response(kind, power, comp, bound, resp)
        else:
            if not r.get("parameters"):
                skipped["no-parameters"] += 1
                continue
            try:
                res = verify(kind, power, comp, bound, r["parameters"])
            except Exception as e:
                key = f"check-raise:{type(e).__name__}"
                tally[key] += 1
                by_type[kind][key] += 1
                if args.list_fail:
                    fails.append((kind, r["power"], comp, r["rational"], key))
                continue
        ok = res["identity_ok"] and res["nonneg"]
        key = "ok" if ok else ("identity-false" if not res["identity_ok"] else "nondefinite")
        tally[key] += 1
        by_type[kind][key] += 1
        if not ok and args.list_fail:
            fails.append((kind, r["power"], comp, r["rational"], key))

    for k, n in tally.most_common():
        print(f"{n:>7}  {k}")
    if len(by_type) > 1:
        print("\nper-type (non-ok only):")
        for t in sorted(by_type):
            bad = {k: n for k, n in by_type[t].items() if k != "ok"}
            if bad:
                print(f"  {t}: {dict(bad)}")
    if skipped:
        print(f"\nskipped: {dict(skipped)}", file=sys.stderr)
    for f in fails:
        print("FAIL", *f)
    # Bad emitted proof or our-side crash = gate failure. WrongDirection /
    # NoSolution / ValueError / EqualClaim are honest outcomes (the corpus
    # deliberately contains false and unprovable claims).
    bad = ("identity-false", "nondefinite", "prove:InternalError")
    return 1 if any(k in bad or k.startswith("check-raise:") for k in tally) else 0


if __name__ == "__main__":
    sys.exit(main())
