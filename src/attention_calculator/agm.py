"""AGM-iteration prover for ``gauss`` and ``varpi`` (plan doc W7).

Second prover for the two lemniscate constants, whose (m, n)-search floor
is the coarsest of all families (varpi '>' cannot even prove varpi > 2.6,
gauss '>' stops near 0.825).  The math is the Gauss AGM: a_0 = sqrt(2),
b_0 = 1,

    a_{k+1} = (a_k + b_k)/2,   b_{k+1} = sqrt(a_k b_k),

and the true iterates satisfy b_n < M < a_n for every finite n, where
M = AGM(1, sqrt(2)) = 1.1981402347...  The site constants are
G = 1/M (gauss) and varpi = pi/M = pi*G, so a rigorous interval on M --
plus one pi bound via the quadlog kernel as oracle for varpi -- proves
both.  Quadratic convergence reaches ~1e-80 in six steps; the decimal
scale is the only error floor, escalated through LADDER until the
comparison decides.

The b-iterates are irrational: each is carried as a rational interval at a
fixed decimal scale D = 10^digits via integer-floor sqrt enclosures
(isqrt on scaled numerators), so every bound is an exact Fraction and the
interval ends always contain the true iterate.

Certificates:

    gauss    {prover:"agm", kind, power, comparison, bound,
              agm_iter, agm_digits, lo, hi}
             certifies power*G comp bound from lo < M < hi (strict: the
             AGM never attains M in finitely many steps)
    varpi    gamma_special-shaped composite DAG, rule "pi_div_agm":
             children pi comp A (quadlog kernel cert) and G comp B (an
             "agm" cert); the transfer needs A*B on R's side with A,B > 0

verify_cert replays the interval AGM over QQ and rechecks every claimed
value; no mpmath anywhere on the proof path.
"""

from fractions import Fraction
from math import isqrt

# Escalating decimal scales for the isqrt enclosures.  64 digits already
# resolves any gap above ~1e-60 (the AGM converges quadratically and is
# never the bottleneck); the rest of the ladder exists for absurdly tight
# bounds.  MAX_DIGITS caps what verify_cert will replay.
LADDER = (64, 512, 4096)
MAX_ITER = 64
MAX_DIGITS = LADDER[-1]


def sqrt_floor(x: Fraction, D: int) -> Fraction:
    """Largest k/D with (k/D)^2 <= x, for x >= 0."""
    return Fraction(isqrt(D * D * x.numerator // x.denominator), D)


def sqrt_ceil(x: Fraction, D: int) -> Fraction:
    """Smallest k/D with (k/D)^2 >= x, for x >= 0."""
    m = D * D * x.numerator
    k = isqrt(m // x.denominator)
    if k * k * x.denominator != m:
        k += 1
    return Fraction(k, D)


def agm_iter(D: int):
    """Yield (lo, hi) with lo < M < hi for n = 0, 1, ... at scale D.

    The a-runs are exact means of rationals; the b-runs carry isqrt
    floor/ceil approximants of sqrt(a*b), so [a_lo, a_hi] and [b_lo, b_hi]
    always contain the true iterates and the enclosure [b_lo, a_hi]
    inherits the AGM sandwich b_n < M < a_n.
    """
    a_lo, a_hi = sqrt_floor(Fraction(2), D), sqrt_ceil(Fraction(2), D)
    b_lo = b_hi = Fraction(1)
    while True:
        yield b_lo, a_hi
        a_lo, a_hi, b_lo, b_hi = (
            (a_lo + b_lo) / 2,
            (a_hi + b_hi) / 2,
            sqrt_floor(a_lo * b_lo, D),
            sqrt_ceil(a_hi * b_hi, D),
        )


def agm_enclosure(n: int, D: int) -> tuple[Fraction, Fraction]:
    """The n-th (lo, hi) pair of agm_iter; verify_cert replays this."""
    it = agm_iter(D)
    for _ in range(n + 1):
        pair = next(it)
    return pair


def flip(comp: str) -> str:
    """Opposite comparison."""
    return "<" if comp == ">" else ">"


def decide_gauss(comp: str, R: Fraction) -> tuple[int, int, Fraction, Fraction] | None:
    """AGM enclosure deciding ``G comp R``; (n, digits, lo, hi) or None.

    G = 1/M sits strictly inside (1/hi, 1/lo): '>' holds iff hi*R <= 1
    (covers R <= 0 too), '<' iff lo*R >= 1.  Once the near end crosses 1/R
    the claim is refuted outright -- lo*R >= 1 for '>' resp. hi*R <= 1 for
    '<' -- so false claims exit at the first scale.
    """
    for digits in LADDER:
        for n, (lo, hi) in enumerate(agm_iter(10**digits)):
            if n > MAX_ITER:
                break
            if comp == ">":
                if hi * R <= 1:
                    return n, digits, lo, hi
                if lo * R >= 1:
                    return None
            else:
                if lo * R >= 1:
                    return n, digits, lo, hi
                if hi * R <= 1:
                    return None
    return None


def converged(digits: int) -> tuple[int, Fraction, Fraction]:
    """Iterate the scale-10^digits AGM to its noise floor; (n, lo, hi).

    The isqrt drift is O(n)/D per end, so once the enclosure width drops
    below 64/D the remaining width is scale noise -- more iterations only
    add drift.  Rigour never depends on this choice; it just picks the
    tightest pair the scale supports.
    """
    D = 10**digits
    for n, (lo, hi) in enumerate(agm_iter(D)):
        if hi - lo <= Fraction(64, D) or n >= MAX_ITER:
            break
    return n, lo, hi


def gauss_cert(
    power: Fraction,
    comp: str,
    bound: Fraction,
    n: int,
    digits: int,
    lo: Fraction,
    hi: Fraction,
) -> dict:
    """Assemble the agm certificate for ``power*G comp bound``."""
    return {
        "prover": "agm",
        "kind": "gauss",
        "power": power,
        "comparison": comp,
        "bound": bound,
        "agm_iter": n,
        "agm_digits": digits,
        "lo": lo,
        "hi": hi,
    }


def pi_cert(comp: str, A: Fraction) -> dict | None:
    """The quadlog kernel's certificate for ``pi comp A``; None if unproved."""
    from . import solve

    try:
        resp = solve.prove("pi", "1", comp, str(A), exact=True)
    except Exception:
        return None
    return resp["certificate"]


def prove_varpi(c: str, R: Fraction, power: Fraction, comp: str, bound: Fraction) -> dict | None:
    """Composite cert for the normalized ``varpi c R``; None if unproved.

    varpi = pi*G: the cert pairs a pi kernel cert for ``pi c A`` with an
    agm cert for ``G c B``, choosing B as the enclosure end and A = R/B so
    the product sits exactly on R -- the strictness the transfer needs is
    inherited from the children.  For '>' with R <= 0 any positive A
    discharges the product, so A = 1.
    """
    if c == "<" and R <= 0:
        return None
    for digits in LADDER:
        n, lo, hi = converged(digits)
        B = Fraction(1) / hi if c == ">" else Fraction(1) / lo
        A = R / B
        if A <= 0:
            A = Fraction(1)
        child_pi = pi_cert(c, A)
        if child_pi is None:
            continue
        child_g = cert_jsonable(gauss_cert(Fraction(1), c, B, n, digits, lo, hi))
        return {
            "prover": "composite",
            "kind": "varpi",
            "comp": comp,
            "q": str(power),
            "p": str(bound),
            "rule": "pi_div_agm",
            "witness": {"A": str(A), "B": str(B)},
            "expect": [
                {"kind": "pi", "power": "1", "comp": c, "bound": str(A)},
                {"kind": "gauss", "power": "1", "comp": c, "bound": str(B)},
            ],
            "children": [child_pi, child_g],
        }
    return None


def prove(kind: str, power: Fraction, comp: str, bound: Fraction) -> dict | None:
    """AGM second prover: ``power*C comp bound``, C in {G, varpi}; cert|None.

    power is the coefficient: q*C comp r normalizes to C comp r/q with the
    comparison flipped when q < 0; q = 0 is a rational-only claim and
    rejected.  Returns the certificate dict (agm schema for gauss,
    composite schema for varpi), or None when the ladder cannot decide --
    a false or sub-floor claim.  Direction certification is the caller's
    job (solve.certified_cmp in mode=exact).
    """
    if comp not in (">", "<"):
        raise ValueError(f"bad comparison {comp!r}")
    if power == 0:
        raise ValueError("常数系数不能为0：命题退化为有理数比较")
    c = comp if power > 0 else flip(comp)
    R = bound / power
    if kind == "gauss":
        hit = decide_gauss(c, R)
        if hit is None:
            return None
        return gauss_cert(power, comp, bound, *hit)
    if kind == "varpi":
        return prove_varpi(c, R, power, comp, bound)
    raise ValueError(f"agm prover does not cover kind {kind!r}")


def verify_cert(cert: dict) -> bool:
    """Re-verify an agm certificate over QQ (parsed form; cert_parse first).

    Replays the claimed iteration: the stored (lo, hi) must equal the n-th
    enclosure at scale 10^digits verbatim, so a tampered iteration count,
    scale, or interval end fails.  Then lo < M < hi strictly (the AGM
    never attains M), hence G in (1/hi, 1/lo) strictly, and the claim
    power*G comp bound reduces to one Fraction comparison after
    normalizing to G comp bound/power (flipped for power < 0).
    """
    if cert.get("kind") != "gauss" or cert.get("prover") != "agm":
        return False
    comp = cert["comparison"]
    q, r = cert["power"], cert["bound"]
    n, digits = cert["agm_iter"], cert["agm_digits"]
    lo, hi = cert["lo"], cert["hi"]
    if comp not in (">", "<") or q == 0:
        return False
    if not isinstance(n, int) or not isinstance(digits, int):
        return False
    if not 0 <= n <= MAX_ITER or not 1 <= digits <= MAX_DIGITS:
        return False
    if not 0 < lo < hi or agm_enclosure(n, 10**digits) != (lo, hi):
        return False
    c = comp if q > 0 else flip(comp)
    R = r / q
    return hi * R <= 1 if c == ">" else lo * R >= 1


_FRAC_KEYS = ("power", "bound", "lo", "hi")


def cert_jsonable(cert: dict) -> dict:
    """JSON-safe copy of an agm certificate: every Fraction as "n/d" text."""
    return {**cert, **{k: str(cert[k]) for k in _FRAC_KEYS}}


def cert_parse(cert: dict) -> dict:
    """Inverse of cert_jsonable; also accepts already-Fraction fields."""
    return {**cert, **{k: Fraction(cert[k]) for k in _FRAC_KEYS}}
