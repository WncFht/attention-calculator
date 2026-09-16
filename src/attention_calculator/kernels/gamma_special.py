"""Composite Γ-value proofs: gamma14, gamma34, gamma12 (exact-mode only).

No beta-moment kernel reaches Γ(1/4) itself — the quarter-lattice parity
lemma locks g := Γ(1/4) into even powers (docs/2026-09-16-w3-research-
gamma-quarter.md).  What exists is the author's own mechanism: g² = S·T
with S := 2ϖ (the varpi kernel's constant doubled) and T := √(2π) (a pi
bound on π vs t²/2), so

    q·g ⋚ r   ⟺   g ⋚ R,   R := r/q   (comparison flips when q < 0)

reduces to one rational bound per factor plus an exact square-root
transfer.  h := Γ(3/4) = V/g with V := π√2 (the pi_sqrt2 kernel) needs one
more level: h ⋚ R ⟺ V ⋚ R·g, discharged as V ⋚ A and g ⋛ B with
A ≷ R·B.  s := Γ(1/2) = √π is a single pi bound through a square root.

Each certificate is a small proof DAG: ``rule`` names the transfer lemma,
``witness`` holds the chosen rational sub-bounds, ``expect`` states every
child claim in /calculate units, and ``children`` carries the child
certificates verbatim.  verify_cert re-verifies each child recursively
(certificate.verify_cert dispatch) and rechecks the transfer arithmetic
over QQ — the only non-QQ content is the constant identities
g² = S·T, h = V/g, s² = π, recorded once here:

    sqrt_mul      g ⋚ R:   2ϖ ⋚ s1 (varpi),  π ⋚ t2²/2 (pi)
                  '>': s1, t2 > 0 and s1·t2 ≥ R²;  '<': same with ≤
    sqrt_div      h ⋚ R:   g ⋛ B (gamma14 cert),  V ⋚ A (pi_sqrt2)
                  '>': A ≥ R·B;  '<': A ≤ R·B;  B > 0 both ways
    sqrt          s ⋚ R:   π ⋚ R² (pi), R > 0
    pos_transfer  C > R with R ≤ 0: child cert proves C > w, w > 0
    pi_div_agm    ϖ ⋚ R:   π ⋚ A (pi kernel), G ⋚ B (agm cert)
                  '>': A·B ≥ R;  '<': A·B ≤ R;  A, B > 0 both ways
"""

from fractions import Fraction

import mpmath as mp

from ..engine import NoSolution

DPS_LADDER = (80, 300, 1200)

_RULES = ("sqrt_mul", "sqrt_div", "sqrt", "pos_transfer", "pi_div_agm")


def _flip(comp: str) -> str:
    """Opposite comparison."""
    return "<" if comp == ">" else ">"


def _snap(x) -> Fraction:
    """Snap an mpf to a Fraction at ~40 significant digits."""
    return Fraction(mp.nstr(x, 40))


def _snap_candidates(x) -> list:
    """Interior-snap candidates, pretty rationals first: the child oracle and
    the exact product checks are the real gatekeepers, so a cheap
    small-denominator snap is tried before the full-precision one."""
    pretty = Fraction(mp.nstr(x, 30)).limit_denominator(10**6)
    fine = _snap(x)
    return [pretty, fine] if pretty != fine else [fine]


def _const_mpf(name: str):
    """Unit constants at ambient precision: g, h, s, S, T, V."""
    if name == "g":
        return mp.gamma(mp.mpf(1) / 4)
    if name == "h":
        return mp.gamma(mp.mpf(3) / 4)
    if name == "s":
        return mp.sqrt(mp.pi)
    if name == "S":
        return mp.gamma(mp.mpf(1) / 4) ** 2 / mp.sqrt(2 * mp.pi)
    if name == "T":
        return mp.sqrt(2 * mp.pi)
    return mp.pi * mp.sqrt(2)  # V


def _child_cert(kind: str, power: Fraction, comp: str, bound: Fraction):
    """solve.prove a child claim; (certificate, ok)."""
    from .. import solve

    try:
        resp = solve.prove(kind, str(power), comp, str(bound), exact=True)
    except Exception:
        return None, False
    return resp["certificate"], True


# ------------------------------------------------------------------ builders


def _cert(
    kind: str, comp: str, R: Fraction, rule: str, witness: dict, expect: list, children: list
) -> dict:
    """Assemble a composite cert for the normalized claim ``C(kind) comp R``."""
    return {
        "prover": "composite",
        "kind": kind,
        "comp": comp,
        "q": "1",
        "p": str(R),
        "rule": rule,
        "witness": {k: str(v) for k, v in witness.items()},
        "expect": [
            {"kind": k, "power": str(pw), "comp": c, "bound": str(b)} for k, pw, c, b in expect
        ],
        "children": children,
    }


def _sqrt_mul(comp: str, R: Fraction) -> dict | None:
    """Witness + child certs for ``g comp R`` via g² = S·T; None if unproved."""
    for dps in DPS_LADDER:
        with mp.workdps(dps):
            Sv, Tv = _const_mpf("S"), _const_mpf("T")
            Rm = mp.mpf(R.numerator) / R.denominator
            R2 = Rm * Rm
            if comp == ">":
                lo, hi = R2 / Sv, Tv  # t2 must exceed R²/S, stay below T
            else:
                lo, hi = Tv, R2 / Sv
            if not lo < hi:
                return None  # numerically infeasible at this precision
            for frac_t in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)):
                t2_pt = lo + (hi - lo) * mp.mpf(frac_t.numerator) / frac_t.denominator
                for t2 in _snap_candidates(t2_pt):
                    if t2 <= 0:
                        continue
                    # t2 inside (lo,hi) is only snap-accurate; the child oracle
                    # and the exact product check below are the real gatekeepers
                    half = t2 * t2 / 2
                    pi_cert, ok = _child_cert("pi", Fraction(1), comp, half)
                    if not ok:
                        continue
                    s1_lim = Fraction(R) * R / t2  # s1 range endpoint from the product
                    if comp == ">":
                        s_lo, s_hi = s1_lim, _snap(Sv)  # s1 ∈ [R²/t2, S)
                    else:
                        s_lo, s_hi = _snap(Sv), s1_lim  # s1 ∈ (S, R²/t2]
                    if not s_lo < s_hi:
                        continue
                    for frac_s in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)):
                        for s1 in _snap_candidates(s_lo + (s_hi - s_lo) * frac_s):
                            if s1 <= 0:
                                continue
                            if comp == ">" and not s1 * t2 >= R * R:
                                continue
                            if comp == "<" and not s1 * t2 <= R * R:
                                continue
                            vp_cert, ok = _child_cert("varpi", Fraction(2), comp, s1)
                            if not ok:
                                continue
                            return _cert(
                                "gamma14",
                                comp,
                                R,
                                "sqrt_mul",
                                {"s1": s1, "t2": t2},
                                [
                                    ("varpi", Fraction(2), comp, s1),
                                    ("pi", Fraction(1), comp, half),
                                ],
                                [vp_cert, pi_cert],
                            )
    return None


def _sqrt_div(comp: str, R: Fraction) -> dict | None:
    """Witness + child certs for ``h comp R`` via h = V/g."""
    for dps in DPS_LADDER:
        with mp.workdps(dps):
            gv, Vv = _const_mpf("g"), _const_mpf("V")
            Rm = mp.mpf(R.numerator) / R.denominator
            if comp == ">":
                lo, hi = gv, Vv / Rm  # B in (g, V/R)
            else:
                lo, hi = Vv / Rm, gv  # B in (V/R, g)
            if not lo < hi:
                return None
            for frac_b in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 4)):
                b_pt = lo + (hi - lo) * mp.mpf(frac_b.numerator) / frac_b.denominator
                for B in _snap_candidates(b_pt):
                    if B <= 0:
                        continue
                    g_cert, ok = _child_cert("gamma14", Fraction(1), _flip(comp), B)
                    if not ok:
                        continue
                    # A must sit on V's side of R·B: midpoint of (R·B, V) for
                    # '>' resp. (V, R·B) for '<'; the exact side check gates it
                    with mp.workdps(dps):
                        Am = (
                            mp.mpf(R.numerator)
                            / R.denominator
                            * (mp.mpf(B.numerator) / B.denominator)
                        )
                        a_pt = (Am + Vv) / 2
                    for A in _snap_candidates(a_pt):
                        if comp == ">" and not A >= R * B:
                            continue
                        if comp == "<" and not A <= R * B:
                            continue
                        v_cert, ok = _child_cert("pi_sqrt2", Fraction(1), comp, A)
                        if not ok:
                            continue
                        return _cert(
                            "gamma34",
                            comp,
                            R,
                            "sqrt_div",
                            {"B": B, "A": A},
                            [
                                ("gamma14", Fraction(1), _flip(comp), B),
                                ("pi_sqrt2", Fraction(1), comp, A),
                            ],
                            [g_cert, v_cert],
                        )
    return None


def _sqrt(comp: str, R: Fraction) -> dict | None:
    """Witness + child cert for ``s comp R`` via s² = π; R > 0."""
    half = R * R
    pi_cert, ok = _child_cert("pi", Fraction(1), comp, half)
    if not ok:
        return None
    return _cert(
        "gamma12",
        comp,
        R,
        "sqrt",
        {},
        [("pi", Fraction(1), comp, half)],
        [pi_cert],
    )


def _pos_transfer(kind: str, R: Fraction) -> dict | None:
    """``C > R`` for R ≤ 0 via a child proving C > w for some w > 0."""
    for dps in DPS_LADDER:
        with mp.workdps(dps):
            w = _snap(_const_mpf({"gamma14": "g", "gamma34": "h", "gamma12": "s"}[kind]) / 2)
        child, ok = _child_cert(kind, Fraction(1), ">", w)
        if not ok:
            continue
        return _cert(
            kind,
            ">",
            R,
            "pos_transfer",
            {"w": w},
            [(kind, Fraction(1), ">", w)],
            [child],
        )
    return None


def _prove_const(kind: str, comp: str, R: Fraction) -> dict | None:
    """Composite cert for the normalized unit claim, or None."""
    if comp == ">" and R <= 0:
        return _pos_transfer(kind, R)
    if comp == "<" and R <= 0:
        return None  # false claim — certified_cmp should have caught it
    if kind == "gamma14":
        return _sqrt_mul(comp, R)
    if kind == "gamma34":
        return _sqrt_div(comp, R)
    return _sqrt(comp, R)


def prove(kind: str, power: Fraction, comp: str, bound: Fraction, exact: bool = True) -> dict:
    """prove("gamma14", q, comp, bound) -> composite /calculate-shaped dict.

    power is the rational coefficient: q·C ⋚ r normalizes to C ⋚ r/q with
    the comparison flipped when q < 0.  q = 0 is a rational-only claim and
    rejected as degenerate.
    """
    if power == 0:
        raise ValueError("Γ型系数不能为0：命题退化为有理数比较")
    R = bound / power
    c = comp if power > 0 else _flip(comp)
    cert = _prove_const(kind, c, R)
    if cert is None:
        raise NoSolution
    # the cert's own claim fields record the literal request; the normalized
    # unit claim is recoverable by re-deriving R = p/q inside the verifier
    cert["q"], cert["p"], cert["comp"] = str(power), str(bound), comp
    return {"type": kind, "prover": "composite", "certificate": cert}


# ------------------------------------------------------------------ verifier


def child_claim(cert: dict) -> tuple[str, Fraction, str, Fraction]:
    """The (kind, power, comp, bound) a certificate certifies, in /calculate units."""
    if cert.get("prover") in ("composite", "euler_gamma"):
        return cert["kind"], Fraction(cert["q"]), cert["comp"], Fraction(cert["p"])
    if "serr" in cert:
        q = Fraction(cert["q"])
        power = q + 1 if cert["kind"] == "ln_q" else q
        return cert["kind"], power, cert["comp"], Fraction(cert["p"])
    return cert["kind"], Fraction(cert["power"]), cert["comparison"], Fraction(cert["bound"])


def verify_cert(cert: dict) -> bool:
    """Re-verify a composite certificate: children recursive + transfer over QQ."""
    try:
        kind, comp = cert["kind"], cert["comp"]
        q, p = Fraction(cert["q"]), Fraction(cert["p"])
        rule = cert["rule"]
        witness = {k: Fraction(v) for k, v in cert["witness"].items()}
        expect = [
            (e["kind"], Fraction(e["power"]), e["comp"], Fraction(e["bound"]))
            for e in cert["expect"]
        ]
        children = cert["children"]
    except (KeyError, TypeError, ValueError):
        return False
    if kind not in ("gamma14", "gamma34", "gamma12", "varpi") or comp not in (">", "<"):
        return False
    if rule not in _RULES or q == 0 or len(expect) != len(children):
        return False

    from ..certificate import verify_cert as verify_any

    for e, ch in zip(expect, children, strict=True):
        try:
            if child_claim(ch) != e or not verify_any(ch):
                return False
        except Exception:
            return False

    R = p / q
    c = comp if q > 0 else _flip(comp)
    if rule == "pos_transfer":
        w = witness["w"]
        return c == ">" and R <= 0 and w > 0 and expect == [(kind, Fraction(1), ">", w)]
    if rule == "sqrt":
        return kind == "gamma12" and R > 0 and expect == [("pi", Fraction(1), c, R * R)]
    if rule == "sqrt_mul":
        s1, t2 = witness["s1"], witness["t2"]
        if kind != "gamma14" or s1 <= 0 or t2 <= 0:
            return False
        prod = s1 * t2
        side_ok = prod >= R * R if c == ">" else (R > 0 and prod <= R * R)
        return side_ok and expect == [
            ("varpi", Fraction(2), c, s1),
            ("pi", Fraction(1), c, t2 * t2 / 2),
        ]
    if rule == "pi_div_agm":
        # varpi = pi*G: children pi comp A (quadlog cert) and G comp B
        # (agm cert); both strict, so the product is strict already at
        # A*B == R -- the boundary form is what the prover emits
        A, B = witness["A"], witness["B"]
        if kind != "varpi" or A <= 0 or B <= 0:
            return False
        side_ok = A * B >= R if c == ">" else A * B <= R
        return side_ok and expect == [
            ("pi", Fraction(1), c, A),
            ("gauss", Fraction(1), c, B),
        ]
    # sqrt_div
    B, A = witness["B"], witness["A"]
    if kind != "gamma34" or B <= 0 or R <= 0:
        return False
    side_ok = A >= R * B if c == ">" else A <= R * B
    return side_ok and expect == [
        ("gamma14", Fraction(1), _flip(c), B),
        ("pi_sqrt2", Fraction(1), c, A),
    ]
