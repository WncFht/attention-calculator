"""Exact-mode tests for the gauss-erf family (W3 types, exact-only).

gaussint_q / dawson_q / erfiint_q share one moment machine on [0,1]
(kernels.gauss_erf): span{C_q, e^{±q²}, 1} with deg-2 P. Every emitted proof
is re-verified by exact_check.gauss_erf.check — dict-equal moments plus the
sign certificate. Constants are evaluated by quadrature of the defining
integrals: this box's mp.erf/mp.erfi leak ~1e-17 float64 noise, while
mp.quad agrees with the incomplete-gamma closed form at full precision.

Tight-bound coverage within m,n <= 10: 1e-6 floors/ceilings prove for
gaussint_q up to q = 3/2 (q = 2 needs more budget), dawson_q up to q = 3,
erfiint_q up to q = 3/2 (H(2) ~= 16.45 still proves at 1e-3, q = 3 is
honestly NoSolution — the window closes as the kernel steepens).
"""

from fractions import Fraction

import mpmath as mp
import pytest

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify
from attention_calculator.exact_check.gauss_erf import check
from attention_calculator.kernels import EXACT_TYPES
from attention_calculator.kernels.gauss_erf import prove, render_equation

mp.dps = 60
_ZERO = mp.mpf(0)

KINDS = ("gaussint_q", "dawson_q", "erfiint_q")


def const_mpf(kind: str, q: Fraction):
    """C(q) via mp.quad of the defining integral (mpf endpoints)."""
    mqf = mp.mpf(q.numerator) / q.denominator
    if kind == "gaussint_q":
        return mp.quad(lambda t: mp.exp(-t * t), [_ZERO, mqf])
    h = mp.quad(lambda t: mp.exp(t * t), [_ZERO, mqf])
    return mp.exp(-mqf * mqf) * h if kind == "dawson_q" else h


def tight_bound(kind: str, q: Fraction, comp: str, digs: int) -> Fraction:
    """``digs``-place floor of C(q) for '>', ceil for '<' (floor is exact for
    negative C too: claims on q < 0 take negative bounds)."""
    c = const_mpf(kind, q)
    v = int(mp.floor(c * 10**digs))
    out = Fraction(v, 10**digs)
    return out + Fraction(1, 10**digs) if comp == "<" else out


def assert_verified(kind: str, q: str, comp: str, bound: Fraction) -> dict:
    """prove() must emit and check() must confirm the exact identity."""
    resp = prove(kind, Fraction(q), comp, bound)
    res = check(kind, Fraction(q), comp, bound, resp["parameters"])
    assert res["identity_ok"], (kind, q, comp, bound, res)
    assert res["nonneg"], (kind, q, comp, bound, res)
    return resp


@pytest.mark.parametrize(
    ("kind", "q"),
    [
        ("gaussint_q", "1/4"),
        ("gaussint_q", "1/2"),
        ("gaussint_q", "1"),
        ("gaussint_q", "3/2"),
        ("gaussint_q", "2"),
        ("dawson_q", "1/4"),
        ("dawson_q", "1/2"),
        ("dawson_q", "1"),
        ("dawson_q", "3/2"),
        ("dawson_q", "2"),
        ("dawson_q", "3"),
        ("erfiint_q", "1/4"),
        ("erfiint_q", "1/2"),
        ("erfiint_q", "1"),
        ("erfiint_q", "3/2"),
        ("erfiint_q", "2"),
    ],
)
@pytest.mark.parametrize("comp", [">", "<"])
def test_tight_bounds_both_directions(kind, q, comp):
    """1e-3 floor/ceil of C(q) proves in both directions."""
    bound = tight_bound(kind, Fraction(q), comp, 3)
    assert_verified(kind, q, comp, bound)


@pytest.mark.parametrize(
    ("kind", "q", "digs"),
    [
        ("gaussint_q", "1/2", 6),
        ("gaussint_q", "1", 6),
        ("gaussint_q", "3/2", 6),
        ("dawson_q", "1", 6),
        ("dawson_q", "3/2", 6),
        ("dawson_q", "3", 6),
        ("erfiint_q", "1/2", 6),
        ("erfiint_q", "1", 6),
        ("erfiint_q", "3/2", 6),
    ],
)
@pytest.mark.parametrize("comp", [">", "<"])
def test_microtight_bounds(kind, q, digs, comp):
    """1e-6 convergent bounds still prove within the exponent budget."""
    bound = tight_bound(kind, Fraction(q), comp, digs)
    assert_verified(kind, q, comp, bound)


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        # the research note's verified identities: G(1) - 7/10 lands at (2,0),
        # F(1) - 1/2 at (0,1) with integrand x(1-x)^2 e^{x^2-1}
        ("gaussint_q", "1", ">", "7/10"),
        ("gaussint_q", "1", "<", "3/4"),
        ("gaussint_q", "1", ">", "74682413/100000000"),  # 1e-8 floor of G(1)
        ("dawson_q", "1", ">", "1/2"),
        ("dawson_q", "1", "<", "53807952/100000000"),  # 1e-8 ceil of F(1)
        ("dawson_q", "2", ">", "3/10"),
        ("erfiint_q", "1", ">", "7/5"),
        ("erfiint_q", "1", "<", "3/2"),
        ("erfiint_q", "1/2", ">", "1/2"),
    ],
)
def test_named_bounds_verify(kind, q, comp, bound):
    assert_verified(kind, q, comp, Fraction(bound))


def test_doc_identity_lands_at_2_0():
    """G(1) - 7/10 = int x^2(4-7x+4x^2)e^{-x^2}/5 dx, verified in the note."""
    resp = assert_verified("gaussint_q", "1", ">", Fraction(7, 10))
    assert (resp["parameters"]["m"], resp["parameters"]["n"]) == (2, 0)
    assert resp["solution"] == "a = 4/5, b = -7/5, c= 4/5"


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("gaussint_q", "1", ">", "4/5"),  # G(1) ~= 0.7468 < 4/5
        ("gaussint_q", "1", "<", "7/10"),
        ("gaussint_q", "-1", ">", "-7/10"),  # G(-1) = -0.7468 > -0.7 false
        ("gaussint_q", "1", "<", "0"),
        ("dawson_q", "1", ">", "7/10"),  # F(1) ~= 0.538
        ("dawson_q", "1", "<", "1/2"),
        ("dawson_q", "-1", ">", "-1/2"),
        ("erfiint_q", "1", "<", "7/5"),  # H(1) ~= 1.4627
        ("erfiint_q", "1", ">", "3/2"),
        ("erfiint_q", "-1", "<", "-3/2"),
    ],
)
def test_false_claim_rejected(kind, q, comp, bound):
    """A false claim dies as WrongDirection: the solved P is uniformly
    non-positive, certifying the opposite inequality."""
    with pytest.raises(WrongDirection):
        prove(kind, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("gaussint_q", "2", ">", "88208139/100000000"),  # 1e-8 floor, budget ends ~1e-6
        ("gaussint_q", "3", "<", "886208/1000000"),
        ("erfiint_q", "2", "<", "16452628/1000000"),  # 1e-6 ceil, beyond budget
        ("erfiint_q", "3", ">", "1444"),  # H(3) ~= 1444.5, loose but unprovable
    ],
)
def test_beyond_budget_is_honest_no_solution(kind, q, comp, bound):
    with pytest.raises(NoSolution):
        prove(kind, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("gaussint_q", "-1", "<", "-7/10"),  # -0.7468 < -0.7
        ("gaussint_q", "-1", ">", "-4/5"),
        ("gaussint_q", "-1", "<", "0"),  # bound 0 exercises the zero-filter
        ("dawson_q", "-1", "<", "-1/2"),  # -0.538 < -0.5
        ("dawson_q", "-3/2", ">", "-1/2"),  # -0.4282 > -0.5
        ("erfiint_q", "-1", "<", "-7/5"),  # -1.4627 < -1.4
        ("erfiint_q", "-2", ">", "-17"),
    ],
)
def test_negative_q_reduction(kind, q, comp, bound):
    """C is odd in q: claims on negative arguments prove via |q| moments."""
    assert_verified(kind, q, comp, Fraction(bound))


@pytest.mark.parametrize("kind", KINDS)
def test_zero_q_rejected(kind):
    """q = 0 degenerates every kernel of the family (1/q in I_0)."""
    with pytest.raises(ValueError):
        prove(kind, Fraction(0), ">", Fraction(0))


@pytest.mark.parametrize("kind", KINDS)
def test_zero_integrand_guard(kind):
    """All-zero params: poly_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {"m": 0, "n": 0, "au_val": "0", "bu_val": "0", "cu_val": "0", "u_val": "1"}
    res = check(kind, Fraction(1), ">", Fraction(0), params)
    assert res["integrand"] == {}
    assert not res["nonneg"]
    assert not res["identity_ok"]


def test_corrupted_params_fail_identity():
    resp = prove("gaussint_q", Fraction(1), ">", Fraction(7, 10))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check("gaussint_q", Fraction(1), ">", Fraction(7, 10), params)
    assert not res["identity_ok"]


def test_params_against_wrong_comp_fail():
    resp = prove("dawson_q", Fraction(1), ">", Fraction(1, 2))
    res = check("dawson_q", Fraction(1), "<", Fraction(1, 2), resp["parameters"])
    assert not res["identity_ok"]


def test_emitted_params_shape():
    resp = prove("gaussint_q", Fraction(1), ">", Fraction(7, 10))
    p = resp["parameters"]
    assert set(p) == {
        "m",
        "n",
        "a_val",
        "b_val",
        "c_val",
        "au_val",
        "bu_val",
        "cu_val",
        "u_val",
        "unified_form",
    }
    assert resp["type"] == "gaussint_q"
    assert resp["solution"].startswith("a = ")
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u
    assert Fraction(p["c_val"]) == Fraction(p["cu_val"]) / u


def test_render_equation_doc_identity():
    """G(1) > 7/10 renders the note's verified integrand x^2(4-7x+4x^2)e^{-x^2}/5."""
    resp = prove("gaussint_q", Fraction(1), ">", Fraction(7, 10))
    tex = render_equation(resp["parameters"], "gaussint_q", Fraction(1), ">", Fraction(7, 10))
    assert tex == (
        "\\mathrm{GaussInt}\\left(1\\right) - \\dfrac{7}{10} = \\int_0^1 "
        "\\frac{x^{2} \\left(4 x^{2} - 7 x + 4\\right) e^{- x^{2}}}{5} \\mathrm{d} x > 0"
    )


def test_render_equation_dawson_prefactor():
    """dawson's normalized kernel q*e^{q^2(x^2-1)} prints its |q| prefactor."""
    resp = prove("dawson_q", Fraction(-3, 2), ">", Fraction(-1, 2))
    tex = render_equation(resp["parameters"], "dawson_q", Fraction(-3, 2), ">", Fraction(-1, 2))
    assert "\\mathrm{Dawson}\\left(\\dfrac{-3}{2}\\right)" in tex
    assert "e^{\\frac{9 x^{2}}{4} - \\frac{9}{4}}" in tex
    assert tex.endswith("\\mathrm{d} x > 0")


def test_render_equation_erfiint():
    resp = prove("erfiint_q", Fraction(1, 2), ">", Fraction(1, 2))
    tex = render_equation(resp["parameters"], "erfiint_q", Fraction(1, 2), ">", Fraction(1, 2))
    assert tex.startswith(
        "\\mathrm{ErfiInt}\\left(\\dfrac{1}{2}\\right) - \\dfrac{1}{2} = \\int_0^1"
    )
    assert "e^{\\frac{x^{2}}{4}}" in tex


@pytest.mark.parametrize("kind", KINDS)
def test_verify_dispatch(monkeypatch, kind):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, kind, "gauss_erf")
    bound = {"gaussint_q": "7/10", "dawson_q": "1/2", "erfiint_q": "7/5"}[kind]
    resp = prove(kind, Fraction(1), ">", Fraction(bound))
    res = verify(kind, Fraction(1), ">", Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize("kind", KINDS)
def test_full_pipeline_when_registered(kind):
    """Once kernels.EXACT_TYPES / solve.FAMILY / integrand.constant_mpf wire
    the family, solver.prove(exact=True) must emit self-checking proofs."""
    if kind not in solver.FAMILY or kind not in EXACT_TYPES:
        pytest.skip(f"{kind} not yet registered (merge wiring is the leader's)")
    bound = {"gaussint_q": "7/10", "dawson_q": "1/2", "erfiint_q": "7/5"}[kind]
    resp = solver.prove(kind, "1", ">", bound, exact=True)
    res = check(kind, Fraction(1), ">", Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
