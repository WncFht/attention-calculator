"""Aggregate parity + verify outputs into the final benchmark report.

Reads the JSONL emitted by ``parity.py --out`` and (when present) the JSONL
emitted by ``verify.py --out``, then writes a markdown summary per type:
success agreement, parameter exact-match rate, identity-validity rate.

Usage:
    python bench/report.py [--parity bench/out/parity.jsonl]
                           [--verify bench/out/verify.jsonl]
                           [--out bench/out/report.md]
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def load_jsonl(path: str | None) -> list[dict]:
    """Load a JSONL file, returning [] when the path is absent or missing."""
    if not path or not Path(path).exists():
        return []
    return [json.loads(line) for line in open(path) if line.strip()]


def parity_table(rows: list[dict]) -> str:
    """Per-type parity metrics from parity.py output rows."""
    by_type: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_type[r["type"]].append(r)

    lines = [
        "| type | n | site ok | both ok | site ok / ours fail"
        " | site fail / ours ok | param exact | crashes |",
        "|---|---|---|---|---|---|---|---|",
    ]
    tot = defaultdict(int)
    for t, rs in sorted(by_type.items()):
        n = len(rs)
        site_ok = sum(1 for r in rs if r["site_success"])
        both = sum(1 for r in rs if r["site_success"] and r.get("ours_success"))
        s_o_f = sum(1 for r in rs if r["site_success"] and not r.get("ours_success"))
        s_f_o = sum(1 for r in rs if not r["site_success"] and r.get("ours_success"))
        exact = sum(1 for r in rs if r.get("param_match"))
        crash = sum(1 for r in rs if str(r.get("ours_error", "")).startswith("crash"))
        lines.append(f"| {t} | {n} | {site_ok} | {both} | {s_o_f} | {s_f_o} | {exact} | {crash} |")
        for k, v in zip(
            ("n", "site_ok", "both", "s_o_f", "s_f_o", "exact", "crash"),
            (n, site_ok, both, s_o_f, s_f_o, exact, crash),
            strict=True,
        ):
            tot[k] += v
    lines.append(
        f"| **合计** | {tot['n']} | {tot['site_ok']} | {tot['both']} | "
        f"{tot['s_o_f']} | {tot['s_f_o']} | {tot['exact']} | {tot['crash']} |"
    )
    return "\n".join(lines)


def verify_table(rows: list[dict]) -> str:
    """Per-type identity-validity metrics from verify.py output rows."""
    if not rows:
        return "_verify.py 输出缺失，跳过。_"
    by_type: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_type[r["type"]].append(r)
    lines = [
        "| type | n | identity ok | sign ok | mismatch |",
        "|---|---|---|---|---|",
    ]
    for t, rs in sorted(by_type.items()):
        n = len(rs)
        ident = sum(1 for r in rs if r.get("identity_ok"))
        sign = sum(1 for r in rs if r.get("sign_ok"))
        lines.append(f"| {t} | {n} | {ident} | {sign} | {n - min(ident, sign)} |")
    return "\n".join(lines)


def main() -> None:
    """Write the aggregate markdown report."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--parity", default="bench/out/parity.jsonl")
    ap.add_argument("--verify", default="bench/out/verify.jsonl")
    ap.add_argument("--out", default="bench/out/report.md")
    args = ap.parse_args()

    parity_rows = load_jsonl(args.parity)
    verify_rows = load_jsonl(args.verify)

    parts = ["# benchmark 报告\n"]
    if parity_rows:
        parts += ["## 与线上一致性（parity）\n", parity_table(parity_rows), ""]
        bad = [r for r in parity_rows if r["site_success"] and not r.get("ours_success")]
        if bad:
            parts.append(f"## 本站可证而本方未证（{len(bad)} 条，前 30）\n")
            parts += [
                f"- `{r['type']}` {r['power']} {r['comparison']} {r['rational']}"
                f" → {r.get('ours_error')}"
                for r in bad[:30]
            ]
            parts.append("")
    parts += ["## 恒等式有效性（verify）\n", verify_table(verify_rows), ""]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(parts))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
