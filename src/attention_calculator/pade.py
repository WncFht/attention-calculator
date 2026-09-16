"""Padé-interpolation prover for ``ln q`` and ``arctan q`` (plan doc W4).

An independent second proof engine for the ln/arctan bounds the site's
(m, n)-search budget cannot reach.  Topsoe's theorem: for x > 0 the Padé
approximants of ln(1+x) interlace around it -- ``l_n = [n/n]`` is a strict
lower bound, ``u_n = [n+1/n]`` a strict upper bound -- and the error
functions are manifestly non-negative rationals

    s_n(x) = (ln(1+x) - l_n(x))' = c x^{2n}   / ((1+x)  D_n(x)^2)
    t_n(x) = (u_n(x) - ln(1+x))' = c x^{2n+1} / ((1+x)  D~_n(x)^2)

with c > 0 and D(0) = 1 (in this normalization c is observed to be
C(2n,n)^{-2} resp. C(2n+1,n)^{-2}).  For arctan the parity-reduced
approximants ``l_n = [2n/2n]`` (lower) and ``u_n = [2n+1/2n+1]`` (upper)
give ``s_n = c x^{4n}/((1+x^2) D_n^2)``, ``t_n = c x^{4n+2}/((1+x^2) D~_n^2)``.

To prove ``ln(1+q) > p`` pick approximants bracketing p -- ``l_n(q) <= p``
below and ``l_m(q) > p`` above -- and interpolate: ``a + b = 1`` with
``a*l_n(q) + b*l_m(q) = p`` turns

    ln(1+q) - p = a (ln - l_n)(q) + b (ln - l_m)(q)
                = integral_0^q (a s_n + b s_m) dx  >  0

into an integral of a rational function that is strictly positive on
(0, q].  Following the author's tutorial examples (article-audit.md par.5)
the bracket is always the weakest approximant on the needed side against
the first one crossing p: for '>' that is ``l_1`` against the first
``l_m(q) > p`` (reproduces ``ln(8/5) > 47/100`` with weight 4499/4500 on
l_3), for '<' ``u_0`` against the first ``u_m(q) < p`` (reproduces
``arctan 2 < 6/5`` with weight 1749/1918 on u_2).  When no approximant on
the needed side exists (p below l_1 resp. above u_0) a single-term
certificate with a positive residual is emitted instead:

    ln(1+q) - p = integral_0^q s_m dx + resid,   resid = l_m(q) - p >= 0.

The certificate ``{kind, comp, q, p, n, m, a, b, resid, serr}`` is all
Fractions; ``verify_cert`` re-derives the approximants and error functions
and rechecks the whole argument over QQ.
"""

from fractions import Fraction
from functools import cache

from .engine import gauss_solve

Poly = list[Fraction]  # ascending coefficients

# Index budget for the crossing scan; n <= 50 already reaches ~1e-30 for
# moderate q (ln 2: n=20, arctan 2: n=37) and each step costs one Fraction
# linear solve of size ~2n, so scans are a few seconds at worst.
MAX_N = 50


# ------------------------------------------------------------------- polys


def pmul(a: Poly, b: Poly) -> Poly:
    """Multiply two ascending-coefficient polynomials over QQ."""
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def psub(a: Poly, b: Poly) -> Poly:
    """a - b for ascending-coefficient polynomials."""
    n = max(len(a), len(b))
    a = a + [Fraction(0)] * (n - len(a))
    b = b + [Fraction(0)] * (n - len(b))
    return [x - y for x, y in zip(a, b, strict=True)]


def pder(a: Poly) -> Poly:
    """Derivative of an ascending-coefficient polynomial."""
    return [a[k] * k for k in range(1, len(a))]


def peval(a: Poly, x: Fraction) -> Fraction:
    """Horner evaluation at a rational point."""
    out = Fraction(0)
    for c in reversed(a):
        out = out * x + c
    return out


# ---------------------------------------------------------------- Padé solve


def ln_series(kmax: int) -> list[Fraction]:
    """Series coefficients of ln(1+x) = sum_{k>=1} (-1)^{k+1} x^k/k up to x^kmax."""
    return [Fraction(0)] + [Fraction((-1) ** (k + 1), k) for k in range(1, kmax + 1)]


def atan_series(kmax: int) -> list[Fraction]:
    """Series coefficients of arctan(x) = sum (-1)^j x^{2j+1}/(2j+1) up to x^kmax."""
    return [Fraction(0) if k % 2 == 0 else Fraction((-1) ** (k // 2), k) for k in range(kmax + 1)]


def pade_approx(c: list[Fraction], num_deg: int, den_deg: int) -> tuple[Poly, Poly]:
    """[num_deg/den_deg] Padé approximant of f(x) = sum_k c_k x^k, over QQ.

    Solves sum_{j=1..M} q_j c_{k-j} = -c_k for k = L+1..L+M and truncates
    Q*f to degree L.  Returns (P, Q) ascending coefficients with Q(0) = 1;
    P(0) = c_0 = 0 for both ln and arctan.
    """
    L, M = num_deg, den_deg
    rows = [
        [(c[k - j] if k - j >= 0 else Fraction(0)) for j in range(1, M + 1)]
        for k in range(L + 1, L + M + 1)
    ]
    rhs = [-c[k] for k in range(L + 1, L + M + 1)]
    q = [Fraction(1), *gauss_solve(rows, rhs)]
    p = [sum(q[j] * c[k - j] for j in range(min(k, M) + 1)) for k in range(L + 1)]
    return p, q


@cache
def approx(func: str, num_deg: int, den_deg: int) -> tuple[Poly, Poly]:
    """Cached Padé approximant of the named series ('ln' | 'atan')."""
    c = (ln_series if func == "ln" else atan_series)(num_deg + den_deg)
    return pade_approx(c, num_deg, den_deg)


# ------------------------------------------------------------- family table

# Per (site kind, comparison) the approximant degrees, the first index that
# counts as a real bound (l_0 = 0 is trivial and skipped for '>'), and the
# exponent of the error-function monomial.  'upper' selects the t-side
# error identity (u - f)' = -E/(G Q^2).
FAMILY = {
    ("ln_q", ">"): {
        "func": "ln",
        "deg": lambda n: (n, n),
        "start": 1,
        "eexp": lambda n: 2 * n,
        "upper": False,
    },
    ("ln_q", "<"): {
        "func": "ln",
        "deg": lambda n: (n + 1, n),
        "start": 0,
        "eexp": lambda n: 2 * n + 1,
        "upper": True,
    },
    ("arctan_q", ">"): {
        "func": "atan",
        "deg": lambda n: (2 * n, 2 * n),
        "start": 1,
        "eexp": lambda n: 4 * n,
        "upper": False,
    },
    ("arctan_q", "<"): {
        "func": "atan",
        "deg": lambda n: (2 * n + 1, 2 * n + 1),
        "start": 0,
        "eexp": lambda n: 4 * n + 2,
        "upper": True,
    },
}

# f' = 1/G: G = 1+x for ln(1+x), 1+x^2 for arctan(x).
BASE = {"ln": [Fraction(1), Fraction(1)], "atan": [Fraction(1), Fraction(0), Fraction(1)]}


def err_form(P: Poly, Q: Poly, func: str, upper: bool) -> tuple[Fraction, int] | None:
    """Monomial data of the Padé error function, or None if malformed.

    (f - P/Q)' = E/(G Q^2) with E = Q^2 - G (P'Q - P Q').  The lower-bound
    error s_n is +E/(G Q^2) and the upper-bound error t_n is -E/(G Q^2), so
    for the upper family -E is checked instead.  Returns (c, e) when the
    relevant sign of E is a single monomial c x^e with c > 0.
    """
    wronsk = psub(pmul(pder(P), Q), pmul(P, pder(Q)))
    e_poly = psub(pmul(Q, Q), pmul(BASE[func], wronsk))
    if upper:
        e_poly = [-v for v in e_poly]
    nz = [(k, v) for k, v in enumerate(e_poly) if v]
    if len(nz) != 1 or nz[0][1] <= 0:
        return None
    return nz[0][1], nz[0][0]


def approx_value(func: str, spec: dict, n: int, q: Fraction) -> Fraction:
    """Value at q of the n-th approximant of the family selected by spec."""
    P, Q = approx(func, *spec["deg"](n))
    return peval(P, q) / peval(Q, q)


def term_desc(func: str, spec: dict, n: int, w: Fraction) -> dict:
    """Self-auditing description of one certificate term's error function:
    err(x) = c x^e / (G(x) den(x)^2) with the stored approximant weight w."""
    P, Q = approx(func, *spec["deg"](n))
    c, e = err_form(P, Q, func, spec["upper"])
    return {"w": w, "n": n, "c": c, "e": e, "den": Q}


# ------------------------------------------------------------------- prover


def prove_point(kind: str, q: Fraction, comp: str, p: Fraction, max_n: int) -> dict | None:
    """Prove the claim at the Padé evaluation point q; cert dict or None.

    Scans the one-sided approximant family (lower bounds for '>', upper for
    '<'): the certificate pairs the first index still on the needed side of
    p (n, weight a) with the first index crossing p (m, weight b), both
    nonneg with a + b = 1 and a A_n(q) + b A_m(q) = p.  If the very first
    index already crosses p a single-term certificate is emitted whose
    value-side slack is recorded in resid (still >= 0).  Returns None when
    no approximant crosses p within max_n -- the claim is then either false
    or too tight for the index budget.
    """
    spec = FAMILY[(kind, comp)]
    func, start = spec["func"], spec["start"]

    # cheap refutation: a small opposite-side approximant already settling
    # the direction means the claim is false and a deep scan is pointless
    opp = FAMILY[(kind, "<" if comp == ">" else ">")]
    veto = approx_value(func, opp, max(opp["start"], 8), q)
    if (veto <= p) if comp == ">" else (veto >= p):
        return None

    weak = cross = None
    vweak = vcross = Fraction(0)
    for k in range(start, max_n + 1):
        v = approx_value(func, spec, k, q)
        if comp == ">":
            if v <= p and weak is None:
                weak, vweak = k, v
            if v > p:
                cross, vcross = k, v
                break
        else:
            if v >= p and weak is None:
                weak, vweak = k, v
            if v < p:
                cross, vcross = k, v
                break
    if cross is None:
        return None

    if weak is None:
        n = m = cross
        a, b = Fraction(0), Fraction(1)
    else:
        n, m = weak, cross
        a = (vcross - p) / (vcross - vweak)
        b = (p - vweak) / (vcross - vweak)
    resid = a * vweak + b * vcross - p if comp == ">" else p - a * vweak - b * vcross

    serr = [term_desc(func, spec, i, w) for w, i in ((a, n), (b, m)) if w > 0]
    return {
        "kind": kind,
        "comp": comp,
        "q": q,
        "p": p,
        "n": n,
        "m": m,
        "a": a,
        "b": b,
        "resid": resid,
        "serr": serr,
    }


def prove(
    kind: str, power: Fraction, comp: str, bound: Fraction, max_n: int = MAX_N
) -> dict | None:
    """Dispatch a site claim ``const(power) ⋚ bound`` to the Padé prover.

    Covers ``ln_q`` (constant ln(power), evaluation point q = power - 1)
    and ``arctan_q`` (constant arctan(power), q = power).  Returns the
    certificate dict, or None for out-of-domain input and for claims no
    approximant crosses within max_n (false or too tight).  Direction
    certification is the caller's job (solve.certified_cmp in mode=exact).
    """
    if comp not in (">", "<"):
        raise ValueError(f"bad comparison {comp!r}")
    if kind == "ln_q":
        if power <= 1:
            return None
        return prove_point("ln_q", power - 1, comp, bound, max_n)
    if kind == "arctan_q":
        if power <= 0:
            return None
        return prove_point("arctan_q", power, comp, bound, max_n)
    raise ValueError(f"pade prover does not cover kind {kind!r}")


# ------------------------------------------------------------------ checker


def verify_cert(cert: dict) -> bool:
    """Re-verify a Padé certificate over QQ; every check is exact Fraction.

    Rechecks: weights a, b >= 0 with a + b == 1; resid >= 0; q > 0; every
    nonzero-weight approximant is re-derived and (i) evaluates to the
    bracketing value used on the resid identity a A_n(q) + b A_m(q) - p =
    resid ('>') resp. p - a A_n - b A_m = resid ('<'), (ii) has error
    function a single positive monomial c x^e/(G D^2) with the family's
    expected exponent e, (iii) has an all-nonneg-coefficient denominator,
    hence pole-free on x > 0 so the integrand is continuous on [0, q]; and
    the stored serr description matches the recomputed error data.  Then
    claim - bound = integral_0^q (a s_n + b s_m) + resid with a strictly
    positive integrand on (0, q] and resid >= 0, so the inequality holds.
    """
    spec = FAMILY.get((cert["kind"], cert["comp"]))
    if spec is None:
        return False
    func, comp = spec["func"], cert["comp"]
    q, p = cert["q"], cert["p"]
    a, b, n, m, resid = cert["a"], cert["b"], cert["n"], cert["m"], cert["resid"]
    if q <= 0 or a < 0 or b < 0 or a + b != 1 or resid < 0:
        return False

    sval = Fraction(0)
    serr = []
    for w, i in ((a, n), (b, m)):
        if w == 0:
            continue
        if i < 0:
            return False
        P, Q = approx(func, *spec["deg"](i))
        if P[0] != 0 or Q[0] <= 0 or any(v < 0 for v in Q):
            return False
        form = err_form(P, Q, func, spec["upper"])
        if form is None or form[1] != spec["eexp"](i):
            return False
        sval += w * peval(P, q) / peval(Q, q)
        serr.append({"w": w, "n": i, "c": form[0], "e": form[1], "den": Q})
    expect = sval - p if comp == ">" else p - sval
    return resid == expect and serr == cert["serr"]


# ------------------------------------------------------------------ wire JSON

FRAC_KEYS = ("q", "p", "a", "b", "resid")


def cert_jsonable(cert: dict) -> dict:
    """JSON-safe copy of a Padé certificate: every Fraction as "n/d" text."""
    out = dict(cert)
    for k in FRAC_KEYS:
        out[k] = str(cert[k])
    out["serr"] = [
        {**t, "w": str(t["w"]), "c": str(t["c"]), "den": [str(v) for v in t["den"]]}
        for t in cert["serr"]
    ]
    return out


def cert_parse(cert: dict) -> dict:
    """Inverse of cert_jsonable; also accepts already-Fraction fields."""
    out = dict(cert)
    for k in FRAC_KEYS:
        out[k] = Fraction(cert[k])
    out["serr"] = [
        {**t, "w": Fraction(t["w"]), "c": Fraction(t["c"]), "den": [Fraction(v) for v in t["den"]]}
        for t in cert["serr"]
    ]
    return out
