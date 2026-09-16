"""Exhaustive (processing-order x resid-rule) search for the '<' sum cases.

Terms = sum terms (a product's factors collapse to one term whose bound is
the product of factor bounds -- factor split analysed separately).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
ROOT = Path(__file__).resolve().parent.parent

VALUES = {"pi": mp.pi, "e": mp.e, "golden": mp.phi, "gamma": mp.euler, "catalan": mp.catalan}
FUNC = {"sin_q": mp.sin, "cos_q": mp.cos, "ln_q": mp.log}


def term_specs(case):
    """Group steps into sum-terms: consecutive steps whose bounds multiply
    (product term) vs add.  We can't see grouping directly, so detect from
    the decomposition latex: '·' inside a term vs '+'/'-' between terms."""
    # split on + or - at top level (decomp latex has no nested +- inside terms)
    return None


def simplest_in(t: mp.mpf, limit: mp.mpf, closed: bool, upper=True):
    """Min-denominator fraction in (t,limit) / (t,limit]; tightest n."""
    for q in range(1, 500000):
        nlo = int(mp.floor(t * q + mp.mpf("1e-35"))) + 1
        if closed:
            nhi = int(mp.floor(limit * q - mp.mpf("1e-35")))
        else:
            nhi = int(mp.ceil(limit * q - mp.mpf("1e-35"))) - 1
        if nlo <= nhi:
            return Fraction(nlo if upper else nhi, q)
    return None


# ---- case table: (problem, display-order terms [(name, value)], site bounds
#      in display order, R)
# Terms listed in DISPLAY order (as in decomposition_latex LHS).
CASES = [
    # problem, [(name, t)], [site bounds display order], R
    ("pi+e<7", [("pi", mp.pi), ("e", mp.e)], ["7/2", "3"], Fraction(7)),
    ("pi+e<6", [("pi", mp.pi), ("e", mp.e)], ["13/4", "11/4"], Fraction(6)),
    ("pi+e<5.87", [("pi", mp.pi), ("e", mp.e)], ["355/113", "30/11"], Fraction(587, 100)),
    ("pi+e<5.86", [("pi", mp.pi), ("e", mp.e)], ["355/113", "193/71"], Fraction(293, 50)),
    ("pi+phi<5", [("phi", mp.phi), ("pi", mp.pi)], ["5/3", "10/3"], Fraction(5)),
    ("pi+phi<4.77", [("phi", mp.phi), ("pi", mp.pi)], ["13/8", "22/7"], Fraction(477, 100)),
    (
        "pi+phi<4.7597",
        [("phi", mp.phi), ("pi", mp.pi)],
        ["233/144", "355/113"],
        Fraction(47597, 10000),
    ),
    ("pi+gamma<4", [("gamma", mp.euler), ("pi", mp.pi)], ["2/3", "10/3"], Fraction(4)),
    ("pi+gamma<3.73", [("gamma", mp.euler), ("pi", mp.pi)], ["7/12", "22/7"], Fraction(373, 100)),
    ("phi+gamma<2.3", [("phi", mp.phi), ("gamma", mp.euler)], ["13/8", "2/3"], Fraction(23, 10)),
    ("phi+gamma<2.2", [("phi", mp.phi), ("gamma", mp.euler)], ["34/21", "11/19"], Fraction(11, 5)),
    ("e+phi<4.5", [("phi", mp.phi), ("e", mp.e)], ["5/3", "11/4"], Fraction(9, 2)),
    ("e+phi<4.34", [("phi", mp.phi), ("e", mp.e)], ["34/21", "68/25"], Fraction(217, 50)),
    ("2*pi+e<9.1", [("2pi", 2 * mp.pi), ("e", mp.e)], ["19/3", "11/4"], Fraction(91, 10)),
    ("pi+3*e<12", [("pi", mp.pi), ("3e", 3 * mp.e)], ["7/2", "17/2"], Fraction(12)),
    ("ln(2)+ln(3)<2", [("ln2", mp.log(2)), ("ln3", mp.log(3))], ["3/4", "5/4"], Fraction(2)),
    ("catalan+pi<4.1", [("C", mp.catalan), ("pi", mp.pi)], ["11/12", "19/6"], Fraction(41, 10)),
    ("pi/2+e<4.3", [("pi/2", mp.pi / 2), ("e", mp.e)], ["11/7", "30/11"], Fraction(43, 10)),
    ("sin(2)+sin(1)<2", [("sin1", mp.sin(1)), ("sin2", mp.sin(2))], ["1", "1"], Fraction(2)),
    ("sin(1)+cos(1)<2", [("sin1", mp.sin(1)), ("cos1", mp.cos(1))], ["1", "3/5"], Fraction(2)),
    (
        "cos(1)+sin(1)<1.39",
        [("sin1", mp.sin(1)), ("cos1", mp.cos(1))],
        ["16/19", "6/11"],
        Fraction(139, 100),
    ),
    (
        "pi+e+gamma<7",
        [("gamma", mp.euler), ("pi", mp.pi), ("e", mp.e)],
        ["2/3", "7/2", "11/4"],
        Fraction(7),
    ),
    (
        "pi+e+phi<8",
        [("phi", mp.phi), ("pi", mp.pi), ("e", mp.e)],
        ["5/3", "7/2", "11/4"],
        Fraction(8),
    ),
    (
        "pi+e+phi<7.5",
        [("phi", mp.phi), ("pi", mp.pi), ("e", mp.e)],
        ["13/8", "22/7", "30/11"],
        Fraction(15, 2),
    ),
    (
        "pi+e+phi<7.49",
        [("phi", mp.phi), ("pi", mp.pi), ("e", mp.e)],
        ["34/21", "22/7", "30/11"],
        Fraction(749, 100),
    ),
    (
        "pi+e+phi+gamma<9",
        [("phi", mp.phi), ("gamma", mp.euler), ("pi", mp.pi), ("e", mp.e)],
        ["5/3", "2/3", "7/2", "3"],
        Fraction(9),
    ),
    (
        "pi+e+phi+gamma<8.1",
        [("phi", mp.phi), ("gamma", mp.euler), ("pi", mp.pi), ("e", mp.e)],
        ["13/8", "3/5", "22/7", "30/11"],
        Fraction(81, 10),
    ),
    ("ln(2)+pi<4", [("ln2", mp.log(2)), ("pi", mp.pi)], ["3/4", "13/4"], Fraction(4)),
    (
        "e*pi+phi<10.2 (product)",
        [("epi", mp.e * mp.pi), ("phi", mp.phi)],
        ["60/7", "57/35"],
        Fraction(51, 5),
    ),
    (
        "e*pi+phi+sin(1)<11",
        [("epi", mp.e * mp.pi), ("phi", mp.phi), ("sin1", mp.sin(1))],
        ["427/50", "8173/5050", "85/101"],
        Fraction(11),
    ),
]


def resid_rules(name):
    def closed(t, resid):
        return simplest_in(t, mp.mpf(resid.numerator) / resid.denominator, True)

    def open_(t, resid):
        return simplest_in(t, mp.mpf(resid.numerator) / resid.denominator, False)

    def exact(t, resid):
        return resid

    def closed_notint(t, resid):
        if resid.denominator == 1:
            return open_(t, resid)
        b = closed(t, resid)
        return b

    def exact_if_int(t, resid):  # exact if integer else closed
        if resid.denominator == 1:
            return resid
        return closed(t, resid)

    return {
        "closed": closed,
        "open": open_,
        "exact": exact,
        "closed_notint": closed_notint,
        "exact_if_int": exact_if_int,
    }[name]


def allocate(ts, R, proc_order, resid_idx_in_order, rule):
    """proc_order: list of term indices processed in that order.
    resid term must be last processed."""
    k = len(ts)
    bounds = [None] * k
    order = list(proc_order)
    assert order[-1] == resid_idx_in_order
    for j, i in enumerate(order):
        later = order[j + 1 :]
        prev_sum = Fraction(0)
        for x in order[:j]:
            prev_sum += bounds[x]
        rest_true = sum(ts[x] for x in later)
        Rf = mp.mpf(R.numerator) / R.denominator
        if i == resid_idx_in_order:
            resid = R - prev_sum
            if mp.mpf(resid.numerator) / resid.denominator <= ts[i]:
                return None
            b = resid_rules(rule)(ts[i], resid)
        else:
            limit = Rf - mp.mpf(prev_sum.numerator) / prev_sum.denominator - rest_true
            b = simplest_in(ts[i], limit, closed=False)
        if b is None:
            return None
        bounds[i] = b
    return bounds


def main():
    rules = ["closed", "open", "exact", "closed_notint", "exact_if_int"]
    # for every case, try every permutation; resid = last processed
    report = []
    for name, terms, bounds_s, R in CASES:
        ts = [t for _, t in terms]
        want = [Fraction(b) for b in bounds_s]
        k = len(ts)
        good = []
        for perm in permutations(range(k)):
            for rule in rules:
                got = allocate(ts, R, perm, perm[-1], rule)
                if got == want:
                    good.append((perm, rule))
        report.append((name, good, terms, want))
        disp = [terms[i][0] for i in range(k)]
        nice = {(tuple(terms[i][0] for i in p), r) for p, r in good}
        print(f"{name:32s} display={disp} ok={sorted(nice)}")

    # intersect rules over cases: which (order-desc, rule) explains all?
    print("\n=== per-rule success count ===")
    for rule in rules:
        n = sum(1 for _, good, _, _ in report if any(r == rule for _, r in good))
        print(rule, n, "/", len(report))


if __name__ == "__main__":
    main()
