"""Replayable probes of the live zhuyidao.net API.

Each batch POSTs /calculate cases serially (>= 0.5 s apart) and appends one
JSON line per HTTP request to bench/data/probes.jsonl. A case with
``image=True`` triggers a follow-up GET /get_integral_image built from the
returned ``parameters`` (needed to read the kernel exponent s, the basis
exponents, and reduced forms like artanh -> ln).

Re-running a batch re-fetches and appends fresh rows; ``tag`` + ``form``
identify the case. Usage:

    uv run python bench/probe_api.py --list
    uv run python bench/probe_api.py mn_order power_semantics
    uv run python bench/probe_api.py --all --delay 0.7
"""

import argparse
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

import requests

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "probes.jsonl"
DELAY = 0.6
LAST_CALL = 0.0

IMAGE_KEYS = ("m", "n", "a_val", "b_val", "c_val", "u_val", "au_val", "bu_val", "cu_val")


def calc(type_, power, comparison, rational, tag, image=False):
    """One /calculate case; ``tag`` marks the probing purpose."""
    return {
        "tag": tag,
        "endpoint": "/calculate",
        "form": {"type": type_, "power": power, "comparison": comparison, "rational": rational},
        "image": image,
    }


def send(case):
    """Send one case after the rate-limit gap; append request+response to OUT."""
    global LAST_CALL
    wait = DELAY - (time.monotonic() - LAST_CALL)
    if wait > 0:
        time.sleep(wait)
    t0 = time.monotonic()
    if case["endpoint"].startswith("/get"):
        r = requests.get(BASE + case["endpoint"], params=case["form"], timeout=60)
    else:
        r = requests.post(BASE + case["endpoint"], data=case["form"], timeout=60)
    LAST_CALL = time.monotonic()
    try:
        body = r.json()
    except ValueError:
        body = r.text
    rec = {
        "ts": datetime.now(UTC).isoformat(),
        "tag": case["tag"],
        "endpoint": case["endpoint"],
        "form": case["form"],
        "http": r.status_code,
        "elapsed_ms": round(1000 * (LAST_CALL - t0), 1),
        "response": body,
    }
    with OUT.open("a") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


def send_image(calc_rec, calc_form, tag):
    """GET /get_integral_image replaying the parameters a browser would send."""
    p = calc_rec["response"]["parameters"]
    form = {k: str(p[k]) for k in IMAGE_KEYS}
    form.update(
        type=calc_form["type"],
        coef=calc_form["power"],
        comparison=calc_form["comparison"],
        rational=calc_form["rational"],
    )
    return send({"tag": tag + ":img", "endpoint": "/get_integral_image", "form": form})


def run_batch(cases):
    """Send cases serially; follow image=True successes with /get_integral_image."""
    n = 0
    for case in cases:
        rec = send(case)
        n += 1
        ok = isinstance(rec["response"], dict) and rec["response"].get("success")
        if case.get("image") and ok:
            send_image(rec, case["form"], case["tag"])
            n += 1
        err = rec["response"].get("error") if isinstance(rec["response"], dict) else rec["http"]
        status = "ok" if ok else err
        print(f"[{case['tag']}] {case['form']} -> {status}", flush=True)
    print(f"{n} requests appended to {OUT}")


# --- batch definitions -----------------------------------------------------


def mn_order():
    """(m,n) growth vs bound tightness; grid-shape images for exotic types."""
    pi = [
        ("3/1", ">"),
        ("4/1", "<"),
        ("31/10", ">"),
        ("32/10", "<"),
        ("22/7", "<"),
        ("223/71", ">"),
        ("314159/100000", ">"),
        ("355/113", "<"),
    ]
    e = [
        ("2/1", ">"),
        ("3/1", "<"),
        ("5/2", ">"),
        ("11/4", "<"),
        ("27/10", ">"),
        ("271828/100000", ">"),
    ]
    ln2 = [
        ("3/5", ">"),
        ("7/10", "<"),
        ("694/1000", "<"),
        ("6931/10000", ">"),
    ]
    cases = [calc("pi", "1", c, r, "mn_order:pi") for r, c in pi]
    cases += [calc("e", "1", c, r, "mn_order:e") for r, c in e]
    cases += [calc("ln_q", "2", c, r, "mn_order:ln2") for r, c in ln2]
    cases += [
        calc("golden", "1", ">", "8/5", "mn_order:golden"),
        # grid-shape evidence: image needed to map returned m,n to exponents
        calc("pi_n", "2", ">", "197/20", "grid:pi_n", image=True),
        calc("zeta3", "1", ">", "6/5", "grid:zeta3", image=True),
        calc("varpi", "1", "<", "8/3", "grid:varpi", image=True),
        calc("gauss", "1", ">", "4/5", "grid:gauss", image=True),
        calc("sin_pi_q", "1/5", ">", "1/2", "grid:sin_pi_q", image=True),
        calc("e_pi", "1", ">", "23/1", "grid:e_pi", image=True),
        calc("gamma", "1", "<", "3/5", "grid:gamma", image=True),
    ]
    return cases


def exponent_caps():
    """Too-tight true-side bounds; error text names the exponent cap."""
    return [
        # boundary brackets: mid-tight bounds for the cap-30 candidates
        calc("pi", "1", "<", "312689/99532", "cap:pi:mid"),
        calc("e", "1", ">", "2718281828459/1000000000000", "cap:e:mid"),
        # near-maximal tightness under the 10^16 input cap -> cap message
        calc("pi", "1", "<", "3141592653589794/1000000000000000", "cap:pi:tight"),
        calc("e", "1", ">", "8154548655371357/3000000000000000", "cap:e:tight"),
        calc("pi_n", "2", ">", "9869604401089358/1000000000000000", "cap:pi_n:tight"),
        calc("e_q", "2", ">", "7389056098930650/1000000000000000", "cap:e_q:tight"),
        calc("e_pi", "1", ">", "2314069263277926/100000000000000", "cap:e_pi:tight"),
        calc("sin_q", "1", ">", "4207354924039482/5000000000000000", "cap:sin_q:tight"),
        calc("sin_pi_q", "1/5", ">", "587785252292473/1000000000000000", "cap:sin_pi_q:tight"),
        calc("arctan_q", "1", ">", "3926990816987241/5000000000000000", "cap:arctan_q:tight"),
    ]


def ln_denominator_power():
    """ln_q / ln_q_square: read s from the equation image; q<1 rejection."""
    qs = [("2", "7/10"), ("3", "11/10"), ("5", "161/100"), ("257", "555/100")]
    cases = [calc("ln_q", q, "<", r, f"ln_s:q={q}", image=True) for q, r in qs]
    cases += [
        calc("ln_q", "2", "<", "2773/4000", "ln_s:q=2:tight", image=True),
        # ln(q)<0 for q<1 but negative rationals are rejected; probe the q range check
        calc("ln_q", "1/2", ">", "0/1", "ln_s:q<1"),
        calc("ln_q_square", "2", "<", "49/100", "ln_s:sq:q=2", image=True),
        calc("ln_q_square", "3", "<", "121/100", "ln_s:sq:q=3", image=True),
    ]
    return cases


def trig_sign_edges():
    """tan/cot at sign-changing q; sin_q/cos_q q-range; degree types."""
    return [
        calc("tan_q", "2", "<", "0/1", "trig:tan2", image=True),
        calc("tan_q", "3/2", "<", "15/1", "trig:tan3/2", image=True),
        calc("cot_q", "2", "<", "0/1", "trig:cot2", image=True),
        calc("sin_q", "4", "<", "0/1", "trig:sin4"),
        calc("sin_q", "7", "<", "7/10", "trig:sin7"),
        calc("cos_q", "4", "<", "0/1", "trig:cos4"),
        calc("sin_q_degree", "45", ">", "7/10", "trig:sin_deg45", image=True),
        calc("cos_q_degree", "60", ">", "2/5", "trig:cos_deg60"),
    ]


def pi_multiple_range():
    """sin_pi_q/cos_pi_q outside (0,1/2): reduction, error, or wrong-direction."""
    return [
        calc("sin_pi_q", "2/3", ">", "4/5", "pimul:sin2/3"),
        calc("sin_pi_q", "1/2", ">", "9/10", "pimul:sin1/2"),
        calc("cos_pi_q", "2/3", "<", "0/1", "pimul:cos2/3"),
        calc("cos_pi_q", "1/2", "<", "1/10", "pimul:cos1/2"),
        calc("cos_pi_q", "3/5", ">", "0/1", "pimul:cos3/5"),
    ]


def power_semantics():
    """power as coefficient vs exponent; fractional power on pi_n/e_q."""
    return [
        calc("pi", "8", "<", "26/1", "power:pi8", image=True),
        calc("e_pi", "2", "<", "47/1", "power:e_pi2", image=True),
        calc("gamma", "2", ">", "1/2", "power:gamma2"),
        calc("pi_n", "5", ">", "300/1", "power:pi_n5", image=True),
        calc("pi_n", "5/2", ">", "100/1", "power:pi_n5/2"),
        calc("e_q", "3/2", "<", "9/2", "power:e_q3/2"),
    ]


def e_param_puzzle():
    """Coordinator's puzzle: e params may not satisfy the rendered identity."""
    return [
        calc("e", "1", ">", "8/3", "epuzzle:8/3", image=True),
        calc("e", "1", ">", "27/10", "epuzzle:27/10", image=True),
        calc("e", "1", ">", "3/1", "epuzzle:3", image=True),
        calc("e", "1", ">", "13/5", "epuzzle:13/5"),
        calc("e", "1", "<", "8/3", "epuzzle:lt8/3", image=True),
        calc("e", "2", ">", "16/3", "epuzzle:coef2"),
        calc("e_q", "2", "<", "15/2", "epuzzle:e_q2", image=True),
        calc("sinh_q", "1", ">", "23/20", "epuzzle:sinh", image=True),
    ]


def error_catalog():
    """Malformed input -> HTTP status + error text."""
    missing_type = {
        "tag": "err:no-type",
        "endpoint": "/calculate",
        "form": {"power": "1", "comparison": "<", "rational": "22/7"},
    }
    return [
        missing_type,
        calc("pi", "1", "<", "1/-2", "err:neg-denom"),
        calc("ln_q", "1", ">", "0/1", "err:ln1"),
        calc("arctan_q", "0", "<", "1/10", "err:atan0"),
    ]


def arc_reductions():
    """artanh/arcoth/arccot: parameters and reduced equation form."""
    return [
        calc("artanh_q", "1/2", "<", "11/20", "arc:artanh", image=True),
        calc("arcoth_q", "2", "<", "3/5", "arc:arcoth", image=True),
        calc("arccot_q", "2", ">", "23/50", "arc:arccot", image=True),
    ]


def type_coverage():
    """Remaining types for unified_form scan and timing stats."""
    return [
        calc("cosh_q", "1", ">", "3/2", "cover:cosh", image=True),
        calc("tanh_q", "1", "<", "4/5", "cover:tanh", image=True),
        calc("coth_q", "1", ">", "13/10", "cover:coth"),
    ]


BATCHES = {
    "mn_order": mn_order,
    "exponent_caps": exponent_caps,
    "ln_denominator_power": ln_denominator_power,
    "trig_sign_edges": trig_sign_edges,
    "pi_multiple_range": pi_multiple_range,
    "power_semantics": power_semantics,
    "e_param_puzzle": e_param_puzzle,
    "error_catalog": error_catalog,
    "arc_reductions": arc_reductions,
    "type_coverage": type_coverage,
}


def main():
    global DELAY
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("batches", nargs="*", choices=BATCHES)
    p.add_argument("--all", action="store_true", help="run every batch in order")
    p.add_argument("--list", action="store_true")
    p.add_argument("--delay", type=float, default=DELAY)
    p.add_argument(
        "--images-for",
        nargs="*",
        metavar="TAG",
        help="refetch /get_integral_image for successful calculates with these tags",
    )
    args = p.parse_args()
    if args.list:
        for name, f in BATCHES.items():
            print(f"{name}: {len(f())} cases")
        return
    DELAY = args.delay
    if args.images_for is not None:
        wanted = set(args.images_for)
        latest = {}
        for line in OUT.open():
            r = json.loads(line)
            if (
                r["endpoint"] == "/calculate"
                and r["tag"] in wanted
                and isinstance(r["response"], dict)
                and r["response"].get("success")
            ):
                latest[r["tag"]] = r
        for tag, rec in latest.items():
            send_image(rec, rec["form"], tag)
            print(f"[{tag}:img] refetched", flush=True)
        return
    names = list(BATCHES) if args.all else args.batches
    if not names:
        p.error("give batch names or --all")
    for name in names:
        run_batch(BATCHES[name]())


if __name__ == "__main__":
    sys.exit(main())
