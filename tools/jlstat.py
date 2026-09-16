#!/usr/bin/env python3
"""Ad-hoc stats over .jsonl records: count by field, filter, project.

Replaces the recurring ``[json.loads(l) for l in open(f)]`` + Counter
one-liner that was the main analysis idiom across sessions.

    .venv/bin/python tools/jlstat.py bench/out/judge.jsonl --key verdict
    .venv/bin/python tools/jlstat.py bench/data/golden.jsonl --key type --where success=false
    .venv/bin/python tools/jlstat.py out.jsonl --key type --key exact.outcome
    .venv/bin/python tools/jlstat.py out.jsonl --where verdict=BUG:crash --show type,power,rational
    cat x.jsonl | .venv/bin/python tools/jlstat.py - --key type
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

_MISSING = object()


def lookup(rec: dict, dotted: str):
    """Dotted-path lookup: 'exact.outcome' walks nested dicts; missing -> _MISSING."""
    cur = rec
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return _MISSING
        cur = cur[part]
    return cur


def parse_literal(text: str):
    """'true'/'false'/'null'/ints compare as their JSON scalars; else string."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def main() -> int:
    """CLI entry: load records, apply --where filters, count or project."""
    ap = argparse.ArgumentParser(description="counter/filter over jsonl records")
    ap.add_argument("file", help="jsonl path ('-' = stdin)")
    ap.add_argument(
        "--key", action="append", default=[], help="field to count by (dotted ok); repeatable"
    )
    ap.add_argument(
        "--where",
        action="append",
        default=[],
        metavar="K=V",
        help="keep records with field == V (JSON literal or string); repeatable",
    )
    ap.add_argument(
        "--show", metavar="F1,F2,...", help="project these fields per record instead of counting"
    )
    args = ap.parse_args()

    if args.file == "-":
        text = sys.stdin.read().splitlines()
    else:
        text = Path(args.file).read_text().splitlines()
    filters = []
    for w in args.where:
        k, _, v = w.partition("=")
        filters.append((k, parse_literal(v)))

    records = []
    for line in text:
        line = line.strip()
        if line:
            rec = json.loads(line)
            if all(lookup(rec, k) == v for k, v in filters):
                records.append(rec)

    if args.show:
        fields = args.show.split(",")
        for rec in records:
            row = {f: (None if lookup(rec, f) is _MISSING else lookup(rec, f)) for f in fields}
            print(json.dumps(row, ensure_ascii=False, default=str))
        print(f"# {len(records)} records", file=sys.stderr)
        return 0

    if not args.key:
        print(f"{len(records)} records")
        return 0

    ctr = Counter(tuple(lookup(r, k) for k in args.key) for r in records)
    for vals, n in sorted(ctr.items(), key=lambda kv: (-kv[1], str(kv[0]))):
        label = "\t".join("(missing)" if v is _MISSING else str(v) for v in vals)
        print(f"{n:>7}\t{label}")
    print(f"{len(records):>7}\ttotal", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
