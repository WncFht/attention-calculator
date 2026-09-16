"""Euler--Maclaurin second prover for the ``gamma`` type (plan doc W7).

The site gamma kernel spends its whole (m, n) budget on moderately tight
bounds because the ln(N+1) leftover is discharged by another budgeted
site search.  This engine instead builds a certified rational enclosure
of Euler's constant directly: Euler--Maclaurin on f(x) = 1/x over [1, N]
gives

    gamma = A_J(N) - ln N + E_J(N),
    A_J(N) = H_N - 1/(2N) + sum_{j=1..J} B_{2j} / (2j N^{2j})

where the remainder has a definite sign and an elementary bound (proof in
docs/2026-09-16-euler-gamma-impl-notes.md): writing s(x) = B_{2J+2}(x) -
B_{2J+2}, the periodic extension s({x}) keeps the sign (-1)^{J+1} on
(0, 1), so

    E_J(N) = -int_N^inf s({x}) / x^{2J+3} dx
    sign(E_J) = (-1)^J,   |E_J| < T := (sum_k |c_k|) / ((2J+2) N^{2J+2})

with c_k the ascending coefficients of s.  Taking N = 2^t leaves ln N =
t ln 2 as the only irrational piece, discharged by two child certificates
bounding ln 2 from below (l) and above (u):

    J even:  gamma in (RA - t u,      RA - t l + T)   -- proves gamma > R
    J odd:   gamma in (RA - t u - T,  RA - t l)       -- proves gamma < R

with RA the Fraction A_J(N) + ln N = H_N - 1/(2N) + sum.  The comparison
holds whenever the binding endpoint clears R, so the tail T never eats
into the proof side -- J only has to make the *feasibility* condition
|E_J| < |gamma - R| hold.

The certificate carries the chosen (n, j, tail, lo, hi), the two child
claims in ``expect``, and the child certificates verbatim in ``children``
(same DAG contract as kernels/gamma_special).  verify_cert replays all
arithmetic over QQ: it recomputes H_N, the Bernoulli numbers (von
Staudt-free exact recurrence), the tail constant and the enclosure,
re-verifies both children recursively, and -- so the sign lemma is not
taken on trust -- certifies the one-signedness of s on [0, 1] itself via
Bernstein coefficients (a sufficient, fully elementary check performed
per certificate).
"""

from fractions import Fraction
from math import comb

import mpmath as mp

from . import pade
from .kernels.gamma_special import child_claim

T_CAP = 10  # N = 2^t runs t = 1..T_CAP
J_CAP = 768  # Bernoulli-term budget per t
LN_MAX_N = 140  # Pade index budget for the ln-2 children

# ----------------------------------------------------------------- QQ tables

B = [Fraction(1)]  # incremental Bernoulli table; B[odd > 1] = 0


def bernoulli(n: int) -> Fraction:
    """B_n over QQ via sum_{k<=m} C(m+1, k) B_k = 0 (B_0 = 1, B_1 = -1/2)."""
    while len(B) <= n:
        m = len(B)
        B.append(-Fraction(1, m + 1) * sum(Fraction(comb(m + 1, k)) * B[k] for k in range(m)))
    return B[n]


def bpoly_minus(n: int) -> list[Fraction]:
    """Ascending coefficients of s(x) = B_n(x) - B_n (constant term drops)."""
    bernoulli(n)
    c = [Fraction(comb(n, k)) * B[n - k] for k in range(n + 1)]
    c[0] -= B[n]
    return c


def bernstein(c: list[Fraction]) -> list[Fraction]:
    """Bernstein coefficients of an ascending-coefficient polynomial on [0,1]."""
    d = len(c) - 1
    return [
        sum(Fraction(comb(i, k), comb(d, k)) * c[k] for k in range(i + 1)) for i in range(d + 1)
    ]


def tail_bound(n_harm: int, j: int) -> Fraction:
    """T = (sum_k C(2J+2,k) |B_k|) / ((2J+2) N^{2J+2}); |E_J(N)| < T."""
    n = 2 * j + 2
    bernoulli(n)
    s = sum(Fraction(comb(n, k)) * abs(B[k]) for k in range(n))
    return s / (n * n_harm**n)


def rat_part(n_harm: int, j: int) -> Fraction:
    """RA = H_N - 1/(2N) + sum_{j<=J} B_{2j}/(2j N^{2j}) over QQ."""
    bernoulli(2 * j + 2)
    h = sum(Fraction(1, i) for i in range(1, n_harm + 1))
    s = sum(B[2 * i] / (2 * i * n_harm ** (2 * i)) for i in range(1, j + 1))
    return h - Fraction(1, 2 * n_harm) + s


def sign_certified(j: int) -> bool:
    """Bernstein check: s = B_{2J+2}(x) - B_{2J+2} keeps sign (-1)^{J+1}.

    All Bernstein coefficients of s on [0, 1] must lie in {0, (-1)^{J+1}};
    then s(x) = sum b_i B_{i,d}(x) has that strict sign on (0, 1) since
    some interior b_i is nonzero.  Sufficient and fully elementary.
    """
    want_pos = j % 2 == 1  # (-1)^{J+1}: J odd -> positive
    return all(b == 0 or (b > 0) == want_pos for b in bernstein(bpoly_minus(2 * j + 2)))


# -------------------------------------------------------------------- prover


def flip(comp: str) -> str:
    """Opposite comparison."""
    return "<" if comp == ">" else ">"


def snap(x, digits: int) -> Fraction:
    """Snap an mpf to a Fraction at the given significant digits."""
    return Fraction(mp.nstr(x, digits))


def ln_cert(comp: str, bound: Fraction, max_n: int) -> dict | None:
    """Wire-form Pade certificate for ``ln 2 comp bound``; None if unproved."""
    cert = pade.prove("ln_q", Fraction(2), comp, bound, max_n=max_n)
    return None if cert is None else pade.cert_jsonable(cert)


def pick_j(n_harm: int, comp: str, want) -> int | None:
    """Smallest right-parity J with tail < want (mpf heuristic), or None.

    J must be even for '>' (E_J > 0) and odd for '<' (E_J < 0).  The scan
    stops when the asymptotic tail starts growing again -- the minimum at
    this N is then already passed and J_CAP is irrelevant.
    """
    best = None
    for j in range(0 if comp == ">" else 1, J_CAP + 1, 2):
        t = tail_bound(n_harm, j)
        if best is not None and t >= best:
            return None  # tail bottomed out above want at this N
        best = t
        if mp.mpf(t.numerator) / t.denominator < want:
            return j if sign_certified(j) else None
    return None


def prove_const(comp: str, bound: Fraction, max_n: int) -> dict | None:
    """Certificate for the normalized claim ``gamma comp bound``, or None."""
    dps = max(80, len(str(bound.numerator)) + len(str(bound.denominator)) + 40)
    with mp.workdps(dps):
        gam = mp.euler
        bn = mp.mpf(bound.numerator) / bound.denominator
        gap = gam - bn if comp == ">" else bn - gam
        ln2 = mp.log(2)
        if gap <= 0:
            return None  # numerically false -- or below heuristic resolution
        for t in range(1, T_CAP + 1):
            n_harm = 1 << t
            j = pick_j(n_harm, comp, gap / 4)
            if j is None:
                continue
            ra = rat_part(n_harm, j)
            lim = (ra - bound) / t  # '>' needs u <= lim; '<' needs l >= lim
            limn = mp.mpf(lim.numerator) / lim.denominator
            room = limn - ln2 if comp == ">" else ln2 - limn
            if room <= 0:
                continue
            digits = max(40, int(-mp.log10(room)) + 15)
            cands = []
            # a non-positive lim makes the '<' lower bound vacuous: any
            # positive l < ln 2 discharges it -- lead with a pretty one
            if comp == "<" and lim <= 0:
                cands.append(Fraction(69, 100))
            for f in (mp.mpf(1) / 2, mp.mpf(3) / 4, mp.mpf(1) / 4, mp.mpf(15) / 16):
                pt = ln2 + f * room if comp == ">" else ln2 - f * room
                cands += [snap(pt, digits), snap(pt, digits).limit_denominator(10**6)]
            for b in cands:
                if comp == ">":
                    if ra - t * b < bound:
                        continue  # candidate past lim (snap noise)
                    u = b
                    low = snap(ln2 - (mp.mpf(u.numerator) / u.denominator - ln2), digits)
                    if low <= 0:
                        low = Fraction(69, 100)  # loose filler, still < ln 2
                else:
                    # low must stay in (0, ln 2); a non-positive lim makes the
                    # tight side vacuous -- any positive low < ln 2 discharges it
                    if b <= 0 or ra - t * b > bound:
                        continue
                    low = b
                    u = snap(ln2 + (ln2 - mp.mpf(low.numerator) / low.denominator), digits)
                hi_child = ln_cert("<", u, max_n)
                if hi_child is None:
                    continue
                lo_child = ln_cert(">", low, max_n)
                if lo_child is None:
                    continue
                tail = tail_bound(n_harm, j)
                e_lo = Fraction(0) if j % 2 == 0 else -tail
                e_hi = tail if j % 2 == 0 else Fraction(0)
                return {
                    "prover": "euler_gamma",
                    "kind": "gamma",
                    "comp": comp,
                    "q": "1",
                    "p": str(bound),
                    "n": n_harm,
                    "j": j,
                    "tail": str(tail),
                    "lo": str(ra - t * u + e_lo),
                    "hi": str(ra - t * low + e_hi),
                    "expect": [
                        {"kind": "ln_q", "power": "2", "comp": ">", "bound": str(low)},
                        {"kind": "ln_q", "power": "2", "comp": "<", "bound": str(u)},
                    ],
                    "children": [lo_child, hi_child],
                }
    return None


def prove(
    kind: str, power: Fraction, comp: str, bound: Fraction, max_n: int = LN_MAX_N
) -> dict | None:
    """prove("gamma", q, comp, r) -> wire-form certificate dict, or None.

    q*gamma ⋚ r normalizes to gamma ⋚ r/q with the comparison flipped when
    q < 0.  q = 0 is a rational-only claim and rejected as degenerate.  None
    means the (t, j, ln-index) budget was honestly exhausted -- the claim is
    false or tighter than the ladder reaches (~1e-700).
    """
    if kind != "gamma":
        raise ValueError(f"euler_gamma prover does not cover kind {kind!r}")
    if comp not in (">", "<"):
        raise ValueError(f"bad comparison {comp!r}")
    if power == 0:
        raise ValueError("γ型系数不能为0：命题退化为有理数比较")
    c = comp if power > 0 else flip(comp)
    cert = prove_const(c, bound / power, max_n)
    if cert is None:
        return None
    # the cert's claim fields record the literal request; the normalized
    # bound is recovered by the verifier as R = p/q with the q<0 flip
    cert["q"], cert["p"], cert["comp"] = str(power), str(bound), comp
    return cert


# -------------------------------------------------------------------- checker


def verify_cert(cert: dict) -> bool:
    """Re-verify an euler_gamma certificate over QQ.

    Replays: schema (kind gamma, comp, literal q/p, n = 2^t, j, tail,
    enclosure lo/hi); both child certificates recursively with claims
    matching expect verbatim; the tail constant and the rational part
    recomputed from exact Bernoulli numbers; the one-signedness of
    s(x) = B_{2j+2}(x) - B_{2j+2} on (0,1) certified by its Bernstein
    coefficients; the enclosure endpoints recomputed; and finally the
    normalized comparison lo >= R resp. hi <= R with R = p/q.
    """
    try:
        kind, comp = cert["kind"], cert["comp"]
        q, p = Fraction(cert["q"]), Fraction(cert["p"])
        n, j = int(cert["n"]), int(cert["j"])
        tail, lo, hi = (Fraction(cert[k]) for k in ("tail", "lo", "hi"))
        expect = [
            (e["kind"], Fraction(e["power"]), e["comp"], Fraction(e["bound"]))
            for e in cert["expect"]
        ]
        children = cert["children"]
    except (KeyError, TypeError, ValueError):
        return False
    if (
        cert.get("prover") != "euler_gamma"
        or kind != "gamma"
        or comp not in (">", "<")
        or q == 0
        or j < 0
        or len(expect) != 2
        or len(children) != 2
        or expect[0][:3] != ("ln_q", Fraction(2), ">")
        or expect[1][:3] != ("ln_q", Fraction(2), "<")
    ):
        return False
    low, u = expect[0][3], expect[1][3]
    t = n.bit_length() - 1
    if n != 1 << t or not 0 < low < u:
        return False

    from .certificate import verify_cert as verify_any

    for e, ch in zip(expect, children, strict=True):
        try:
            if child_claim(ch) != e or not verify_any(ch):
                return False
        except Exception:
            return False

    if tail != tail_bound(n, j) or not sign_certified(j):
        return False
    ra = rat_part(n, j)
    e_lo = Fraction(0) if j % 2 == 0 else -tail
    e_hi = tail if j % 2 == 0 else Fraction(0)
    if (lo, hi) != (ra - t * u + e_lo, ra - t * low + e_hi):
        return False
    c = comp if q > 0 else flip(comp)
    return lo >= p / q if c == ">" else hi <= p / q


# ------------------------------------------------------------------ wire JSON

FRAC_KEYS = ("q", "p", "tail", "lo", "hi")


def cert_jsonable(cert: dict) -> dict:
    """JSON-safe copy; fields are already 'n/d' strings (idempotent)."""
    out = dict(cert)
    for k in FRAC_KEYS:
        out[k] = str(cert[k])
    return out


def cert_parse(cert: dict) -> dict:
    """Fractions restored for q/p/tail/lo/hi; children left in wire form."""
    out = dict(cert)
    for k in FRAC_KEYS:
        out[k] = Fraction(cert[k])
    return out
