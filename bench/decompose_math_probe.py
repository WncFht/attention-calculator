"""Measure the (m,n)-search provable frontier: bound margin vs needed depth.

For each (kind, power, comp) this scans a geometric ladder of rational
bounds approaching the true constant and records the minimal depth
m+n at which ``solve.prove(..., exact=True)`` first emits a proof.
Because engine.mn_order iterates depth-major and a true claim can never
hit a mid-scan WrongDirection, the returned m+n is exactly the minimal
feasible depth for that bound.

Run: ``python bench/decompose_math_probe.py`` — prints a markdown table
consumed by docs/2026-09-16-decompose-math.md.
"""

from fractions import Fraction

import mpmath as mp

from attention_calculator import solve
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.integrand import constant_mpf

DPS = 120

# representative (kind, power) units reachable from the decompose grammar
CASES = [
    ("pi", "1"),
    ("pi", "8"),
    ("e", "1"),
    ("e_q", "2"),
    ("e_q", "1/2"),
    ("e_pi", "1"),
    ("pi_n", "2"),
    ("pi_n", "3"),
    ("pi_n", "5/2"),
    ("ln_q", "2"),
    ("ln_q", "3"),
    ("ln_q", "3/2"),
    ("sin_q", "1"),
    ("cos_q", "1"),
    ("tan_q", "1"),
    ("arctan_q", "1"),
    ("arctan_q", "2"),
    ("sinh_q", "1"),
    ("tanh_q", "1"),
    ("gamma", "1"),
    ("golden", "1"),
    ("catalan", "1"),
    ("zeta3", "1"),
    ("varpi", "1"),
    ("gauss", "1"),
]

# margin ladder: delta = 10^-k
LADDER = [10**-k for k in range(0, 16)]


def rat_below(c, delta):
    """Rational < c-delta, exact decimal snap at 80 significant digits."""
    return Fraction(mp.nstr(c - mp.mpf(delta), 80))


def rat_above(c, delta):
    return Fraction(mp.nstr(c + mp.mpf(delta), 80))


def depth_of(kind, power, comp, r):
    """(depth, status) of the first-hit proof, or the exception name."""
    try:
        resp = solve.prove(kind, power, comp, str(r), exact=True)
    except (NoSolution, WrongDirection) as exc:
        return None, type(exc).__name__
    except Exception as exc:
        return None, type(exc).__name__
    p = resp["parameters"]
    return int(p["m"]) + int(p["n"]), "ok"


def probe(kind, power, comp):
    """One row: minimal depth per margin rung for this constant."""
    q = Fraction(power)
    with mp.workdps(DPS):
        c = constant_mpf(kind, q)
    row = []
    for delta in LADDER:
        with mp.workdps(DPS):
            r = rat_below(c, delta) if comp == ">" else rat_above(c, delta)
        d, status = depth_of(kind, power, comp, r)
        row.append(d if status == "ok" else status)
    return row


def main():
    header = "| kind | power | comp | " + " | ".join(f"1e-{k}" for k in range(16)) + " |"
    sep = "|" + "---|" * (len(LADDER) + 3)
    print(header)
    print(sep)
    for kind, power in CASES:
        for comp in (">", "<"):
            row = probe(kind, power, comp)
            cells = ["-" if d is None else (d if isinstance(d, str) else str(d)) for d in row]
            print(f"| {kind} | {power} | {comp} | " + " | ".join(cells) + " |", flush=True)


if __name__ == "__main__":
    main()
