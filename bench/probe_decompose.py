"""Probe the site's /decompose_inequality endpoint and log every response.

Serial requests, >=0.7s apart, total count capped by --max (site etiquette:
the task budget is 120 requests, keep well under). Responses go to
bench/data/decompose.jsonl (append) and pretty JSON per request under
bench/data/decompose/.

Usage:
    python bench/probe_decompose.py                # run built-in probe list
    python bench/probe_decompose.py "pi+e<6" ...   # run ad-hoc problems
    python bench/probe_decompose.py --skip-existing
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://zhuyidao.net/decompose_inequality"
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "bench" / "data"
JSONL = DATA / "decompose.jsonl"
INTERVAL = 0.7  # seconds between requests (spec: >=0.5s, serial)

# Built-in probe series. Each tuple: (problem, what it tests).
PROBES = [
    # --- known samples (determinism + ground truth) ---
    ("pi^2+8*pi>35", "known sample"),
    ("e*pi+phi+sin(1)<11", "known sample"),
    ("e^2-e^(1/2)-cos(1)>5.2", "known sample"),
    # --- same LHS, sweeping RHS: how bounds move with slack ---
    ("pi^2+8*pi>30", "slack ~5.0"),
    ("pi^2+8*pi>34", "slack ~1.0"),
    ("pi^2+8*pi>34.5", "slack ~0.5"),
    ("pi^2+8*pi>34.9", "slack ~0.1"),
    ("pi^2+8*pi>35.001", "slack ~1e-3"),
    ("pi^2+8*pi>35.002", "slack ~3e-4"),
    ("pi^2+8*pi>36", "false statement -> expect error"),
    # --- two plain terms, slack sweep ---
    ("pi+e<7", "slack ~1.14"),
    ("pi+e<6", "slack ~0.14"),
    ("pi+e<5.87", "slack ~0.0101"),
    ("pi+e<5.86", "slack ~1e-4"),
    ("pi+e<5.8598", "slack ~-7e-5? false"),
    ("e+pi<6", "order swap of pi+e<6"),
    ("pi+phi<5", "slack ~0.24, phi denom watch"),
    ("pi+phi<4.77", "slack ~0.0104"),
    ("pi+phi<4.7597", "slack ~1e-4"),
    # --- single term: direct_basic? ---
    ("pi>3", "direct_basic"),
    ("pi<22/7", "direct_basic"),
    ("8*pi>25", "coefficient term: direct_basic or decompose?"),
    ("2*pi>6", "coefficient, slack 0.283"),
    ("pi^2>9", "power term alone"),
    ("sin(1)<0.85", "function term alone"),
    ("gamma<0.6", "named constant alone"),
    # --- products ---
    ("e*pi<9", "product, slack 0.46"),
    ("e*pi<8.6", "product, slack 0.06"),
    ("e*pi<8.55", "product, slack 0.0103"),
    ("e*pi<8.54", "product bound = 427/50 in sample; alone?"),
    ("e*pi<8.5398", "product, slack ~7e-5"),
    ("pi*e<9", "factor order swap"),
    ("2*e*pi<18", "coefficient on product"),
    ("pi*pi<10", "is pi*pi normalized to pi^2?"),
    ("e*e<8", "e*e vs e^2"),
    ("e*pi*pi<27", "three-factor product, e*pi^2=26.83"),
    ("phi+sin(1)+e*pi<11", "product as LAST term -> residual on product?"),
    ("e*pi+phi<10.2", "product + one term"),
    # --- negative terms / subtraction ---
    ("pi-e>0", "difference, slack 0.42"),
    ("pi-e>0.4", "difference, slack 0.023"),
    ("pi-e>0.42", "difference, slack 0.0033"),
    ("e-pi<0", "negative value <0"),
    ("pi^2-8*pi>-16", "neg result, slack 0.74"),
    ("-pi+8>4", "leading minus, 8-pi=4.86"),
    # --- coefficients ---
    ("8*pi>25.1", "slack 0.0327"),
    ("3*pi>9", "slack 0.42"),
    ("pi+2*e>8.5", "slack 0.078"),
    ("pi/2>1.5", "division in input? pi/2=1.5708"),
    ("pi^2+8*pi>35 ", "trailing space"),
    # --- 3+ plain terms ---
    ("pi+e+phi<8", "3 terms, slack 0.52"),
    ("pi+e+phi<7.5", "3 terms, slack 0.02"),
    ("pi+e+phi<7.49", "3 terms, slack 0.012"),
    ("pi+e+phi+gamma<9", "4 terms, slack 0.94"),
    ("pi+e+phi+gamma<8.1", "4 terms, slack 0.043"),
    # --- syntax coverage for parser ---
    ("ln(2)+pi<4", "ln syntax, slack 0.165"),
    ("ln(2)^2<0.5", "ln squared: type ln_q_square?"),
    ("zeta(3)<1.21", "zeta3 syntax"),
    ("catalan<0.92", "catalan syntax"),
    ("varpi<2.7", "varpi syntax"),
    ("gauss>0.8", "gauss syntax"),
    ("e^pi>23", "e_pi syntax"),
    ("sinh(1)<1.2", "sinh syntax"),
    ("arctan(1)<0.8", "arctan syntax"),
    ("sin(1)+cos(1)<2", "two functions"),
    ("cos(1)+sin(1)<1.39", "order swap, slack 0.007"),
    ("tan(1)<1.56", "tan alone"),
    ("sin(30°)<0.6", "degree syntax? maybe unsupported"),
    ("sin(pi/5)>0.58", "sin_pi_q syntax?"),
    # --- determinism / robustness ---
    ("pi^2+8*pi>35", "repeat: determinism"),
    ("e*pi+phi+sin(1)<11", "repeat: determinism"),
    ("pi+e<6", "repeat: determinism"),
]


def slug(problem: str, idx: int) -> str:
    """Filename-safe slug for a problem string."""
    s = re.sub(r"[^a-z0-9]+", "-", problem.lower()).strip("-")
    return f"{idx:03d}-{s[:60] or 'empty'}"


def post(problem: str) -> dict:
    """One POST to /decompose_inequality; returns parsed JSON or error dict."""
    body = urllib.parse.urlencode({"problem": problem}).encode()
    req = urllib.request.Request(
        BASE, data=body, headers={"Content-Type": "application/x-www-form-urlencoded"}
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode())
    except Exception as exc:
        return {"success": False, "transport_error": f"{type(exc).__name__}: {exc}"}


def main() -> None:
    """Run the probe list serially and record every response."""
    ap = argparse.ArgumentParser()
    ap.add_argument("problems", nargs="*", help="override built-in list")
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        help="skip problems already present in decompose.jsonl",
    )
    ap.add_argument("--max", type=int, default=120, help="request budget cap")
    args = ap.parse_args()

    items = [(p, "adhoc") for p in args.problems] if args.problems else PROBES

    seen: set[str] = set()
    if args.skip_existing and JSONL.exists():
        for line in JSONL.read_text().splitlines():
            if line.strip():
                seen.add(json.loads(line)["problem"])

    outdir = DATA / "decompose"
    outdir.mkdir(parents=True, exist_ok=True)
    sent = 0
    with JSONL.open("a") as log:
        for idx, (problem, note) in enumerate(items):
            if problem in seen:
                print(f"skip (already logged): {problem}")
                continue
            if sent >= args.max:
                print(f"budget {args.max} reached, stopping")
                break
            t0 = time.time()
            res = post(problem)
            elapsed = round(time.time() - t0, 2)
            rec = {"problem": problem, "note": note, "elapsed_s": elapsed, "response": res}
            log.write(json.dumps(rec, ensure_ascii=False) + "\n")
            log.flush()
            (outdir / f"{slug(problem, idx)}.json").write_text(
                json.dumps(rec, ensure_ascii=False, indent=2)
            )
            bounds = [s.get("bound") for s in res.get("steps", [])] if res.get("steps") else None
            print(
                f"[{sent + 1}] {problem} -> success={res.get('success')} "
                f"direct={res.get('direct_basic')} bounds={bounds} "
                f"err={res.get('error') or res.get('transport_error')}"
            )
            sent += 1
            time.sleep(INTERVAL)
    print(f"done, {sent} requests")


if __name__ == "__main__":
    main()
