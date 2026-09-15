"""Trace our 7 divergent edge cases: per-plan exact vs float sign verdicts.

For each request, rebuild the kernel's (plans, target) exactly as prove() does
and print every plan index where exact and float-of-exact verdicts differ, plus
the first sign-definite plan under each rule.
"""

import sys
from fractions import Fraction

sys.path.insert(0, "src")
sys.path.insert(0, "bench")
from sim_mechanisms import basis_and_target

from attention_calculator.engine import mn_order, poly_nonneg, solve_moment

CASES = [
    ("equal:pi-float-lt", "pi", "1", "<", "884279719003555/281474976710656"),
    ("equal:e-float-lt", "e", "1", "<", "6121026514868073/2251799813685248"),
    ("cap:pi-gt-tight", "pi", "1", ">", "6134899525417045/1952799169684491"),
    ("cap:e-gt-tight", "e", "1", ">", "2124008553358849/781379079653017"),
    ("cap:e_q-tight", "e_q", "1", ">", "2124008553358849/781379079653017"),
    ("cap:sin_q-tight", "sin_q", "1", ">", "8199121568317184/9743795943467975"),
    ("cap:ln_q-tight", "ln_q", "2", ">", "1554903831458736/2243252046704767"),
    # controls
    ("pi-gt-cf", "pi", "1", ">", "884279719003555/281474976710656"),
    ("e-gt-cf", "e", "1", ">", "6121026514868073/2251799813685248"),
]


def f_nonneg(c):
    return poly_nonneg([float(v) for v in c])


def f_nonpos(c):
    return poly_nonneg([-float(v) for v in c])


def run(name, kind, pw, comp, rs, maxplans=3000):
    bound = Fraction(rs)
    limit, mk, target, _B = basis_and_target(kind, pw, comp, bound)
    print(f"== {name}: {kind} {pw} {comp} {rs}  limit={limit}")
    first_e = first_f = None
    diffs = []
    for i, (m, n) in enumerate(mn_order(limit)):
        basis = mk(m, n)
        try:
            u = solve_moment(basis, target)
        except ValueError:
            continue
        nn, np_ = poly_nonneg(u), poly_nonneg([-v for v in u])
        fnn, fnp = f_nonneg(u), f_nonpos(u)
        e = "NN" if nn else ("NP" if np_ else "..")
        f = "NN" if fnn else ("NP" if fnp else "..")
        if first_e is None and e != "..":
            first_e = (i, m, n, e)
        if first_f is None and f != "..":
            first_f = (i, m, n, f)
        if e != f:
            diffs.append((i, m, n, e, f, [str(v) for v in u]))
        if first_e and first_f and len(diffs) > 8:
            break
        if i > maxplans:
            break
    print(f"  first exact : {first_e}")
    print(f"  first float : {first_f}")
    for d in diffs[:12]:
        print(f"    plan#{d[0]} ({d[1]},{d[2]}) exact={d[3]} float={d[4]} u={d[5]}")


for c in CASES:
    run(*c)
