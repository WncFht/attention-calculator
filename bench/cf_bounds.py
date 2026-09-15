"""Generate tight rational bounds via continued fractions.

Site validation requires num,den < 1e16, so sub-ulp bounds must be found by
Diophantine approximation: enumerate convergents/semiconvergents of each
constant and keep those within a window of float64(C).
"""
import math
import sys
from fractions import Fraction

sys.path.insert(0, "src")
from mpmath import mp
mp.dps = 120
from attention_calculator.integrand import constant_mpf


def cf_terms(x, n):
    """First n partial quotients of mpf x."""
    out = []
    for _ in range(n):
        a = int(mp.floor(x))
        out.append(a)
        x = x - a
        if x == 0:
            break
        x = 1 / x
    return out


def convergents(cf):
    """All convergents p/q."""
    p0, p1 = mp.mpf(0), mp.mpf(1)
    q0, q1 = mp.mpf(1), mp.mpf(0)
    for a in cf:
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
        yield int(p1), int(q1)


def semiconvergents(cf, maxden):
    """Convergents plus intermediate fractions (p_{k}+t*p_{k-1})/(q_k+t*q_{k-1})."""
    terms = cf
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    out = []
    for a in terms:
        for t in range(1, a + 1):
            p = t * p1 + p0
            q = t * q1 + q0
            if q > maxden:
                break
            out.append(Fraction(p, q))
        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0
    return out


def bounds_near(kind, power, maxden=9999999999999999):
    """Rationals with den<=maxden within half-ulp of float64(C), sorted by gap."""
    C = constant_mpf(kind, Fraction(power))
    Cf = float(C)
    half = math.ulp(Cf) / 2
    cf = cf_terms(C, 400)
    out = []
    for r in semiconvergents(cf, maxden):
        if r.numerator > maxden:
            continue
        rf = float(r)
        if abs(rf - Cf) < half and rf == Cf:
            gap = C - mp.mpf(r.numerator) / r.denominator
            out.append((gap, r))
    out.sort(key=lambda t: abs(t[0]))
    return out


if __name__ == "__main__":
    for kind, power in [("pi", "1"), ("e", "1"), ("sin_q", "1"),
                        ("ln_q", "2"), ("ln_q_square", "2"),
                        ("cos_q", "1"), ("e_q", "3"), ("zeta3", "1"),
                        ("arctan_q", "3"), ("pi_n", "5/2")]:
        C = constant_mpf(kind, Fraction(power))
        Cf = float(C)
        near = bounds_near(kind, power)
        print(f"== {kind} {power}  Cf={Cf!r} (C-Cf={mp.nstr(C-mp.mpf(Fraction(Cf).numerator)/Fraction(Cf).denominator,4)})")
        for gap, r in near[:8]:
            print(f"   gap={mp.nstr(gap,4):>10s}  r={r}  float_eq={float(r)==Cf}")
