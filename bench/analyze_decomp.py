"""Test allocation-rule hypotheses against every recorded decompose case."""

from __future__ import annotations

import json
import re
from fractions import Fraction
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
ROOT = Path(__file__).resolve().parent.parent

VALUES = {
    "pi": mp.pi,
    "e": mp.e,
    "golden": mp.phi,
    "gamma": mp.euler,
    "catalan": mp.catalan,
    "zeta3": mp.zeta(3),
    "varpi": mp.mpf("2.6220575542921198104648395898911194136827549514316"),
    "gauss": mp.mpf("0.83462684167407318630142973337359360395562655639648"),
}
FUNC = {
    "sin_q": mp.sin,
    "cos_q": mp.cos,
    "tan_q": mp.tan,
    "cot_q": mp.cot,
    "sinh_q": mp.sinh,
    "cosh_q": mp.cosh,
    "tanh_q": mp.tanh,
    "coth_q": mp.coth,
    "arctan_q": mp.atan,
    "arccot_q": lambda a: mp.pi / 2 - mp.atan(a),
    "artanh_q": mp.atanh,
    "arcoth_q": lambda a: mp.atanh(1 / a),
    "ln_q": mp.log,
}


def step_value(step: dict) -> mp.mpf:
    """True numeric value of one step's term: coefficient * constant/func."""
    kind = step["type"]
    if kind in VALUES:
        v = VALUES[kind]
    elif kind in FUNC:
        num = re.search(r"(\d+)(?:/(\d+))?", step["label"])
        arg = mp.mpf(int(num.group(1))) / mp.mpf(int(num.group(2) or 1))
        v = FUNC[kind](arg)
    elif kind == "pi_n":
        v = mp.pi ** int(step.get("power", 2))
    else:
        raise ValueError(f"unknown kind {kind}")
    c = Fraction(step["coefficient"])
    return mp.mpf(c.numerator) / c.denominator * v


def simplest(lo: mp.mpf, hi: mp.mpf, closed_hi: bool, want_upper=True):
    """Min-denominator fraction in (lo, hi) or (lo, hi]; tightest (nearest lo
    for upper bounds / nearest hi for lower bounds)."""
    for q in range(1, 300000):
        nlo = int(mp.floor(lo * q + mp.mpf("1e-30"))) + 1
        if closed_hi:
            nhi = int(mp.floor(hi * q - mp.mpf("1e-30")))
        else:
            nhi = int(mp.ceil(hi * q - mp.mpf("1e-30"))) - 1
        if nlo <= nhi:
            n = nlo if want_upper else nhi
            return Fraction(n, q)
    return None


def resid_bound(t: mp.mpf, resid: Fraction, rule: str):
    """Bound for the residual term under various rules."""
    r = mp.mpf(resid.numerator) / resid.denominator
    if r <= t:
        return None
    if rule == "exact":
        return resid
    if rule == "closed":
        return simplest(t, r, closed_hi=True)
    if rule == "open":
        return simplest(t, r, closed_hi=False)
    if rule == "closed_unless_int":
        if resid.denominator == 1:
            return simplest(t, r, closed_hi=False)
        return resid if simplest(t, r, True) == resid else simplest(t, r, True)
    raise ValueError(rule)


def run_order(ts, R, resid_pos, rule):
    """Sequential allocation: resid_pos = index (in processing order) of the
    residual term.  Non-residual terms get simplest-open bound with
    feasibility vs remaining true values; residual gets resid_bound().

    ts: list of mpf true values in PROCESSING order.  R: Fraction.
    Returns list of Fraction bounds in processing order, or None."""
    k = len(ts)
    bounds = [None] * k
    # process non-residual terms in order, resid last
    order = [i for i in range(k) if i != resid_pos] + [resid_pos]
    for j, i in enumerate(order):
        later = [order[x] for x in range(j + 1, k)]
        rest_true = sum(ts[x] for x in later)
        prev = sum(mp.mpf(b.numerator) / b.denominator for b in bounds if b is not None)
        if i == resid_pos:
            resid = R - sum(Fraction(b.numerator, b.denominator) for b in bounds if b is not None)
            b = resid_bound(ts[i], resid, rule)
        else:
            limit = mp.mpf(R.numerator) / R.denominator - prev - rest_true
            b = simplest(ts[i], limit, closed_hi=False)
        if b is None:
            return None
        bounds[i] = b
    return bounds


def main():
    cases = []
    for f in sorted((ROOT / "bench/data/decompose").glob("*.json")):
        rec = json.loads(f.read_text())
        r = rec["response"]
        if not r.get("success") or r.get("direct_basic"):
            continue
        steps = r.get("steps", [])
        if not steps:
            continue
        # reconstruct R from normalized latex tail ">K>0" etc — easier:
        # sum of step bounds compared with problem rhs is messy; instead
        # parse R out of the problem text directly.
        prob = rec["problem"]
        m = re.match(r"^(.*?)([<>])(-?[\d./]+)$", prob.replace(" ", ""))
        if not m:
            continue
        R = Fraction(m.group(3))
        cases.append(
            {
                "problem": prob,
                "comp": m.group(2),
                "R": R,
                "steps": steps,
                "ts": [step_value(s) for s in steps],
                "bounds": [Fraction(s["bound"]) for s in steps],
                "comps": [s["comparison"] for s in steps],
            }
        )
    print(f"{len(cases)} multi-step cases\n")
    for c in cases:
        lhs = sum(c["ts"])
        print(
            f"{c['problem']:28s} R={c['R']!s:>8} lhs={mp.nstr(lhs, 12)} "
            f"bounds={[str(b) for b in c['bounds']]} comps={c['comps']}"
        )


if __name__ == "__main__":
    main()
