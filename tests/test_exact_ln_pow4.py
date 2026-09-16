"""Exact-mode tests for ln_q_quad (kernels.ln_pow, quartic branch).

The type is exact-mode-only and not yet registered in kernels.EXACT_TYPES /
solve.FAMILY (merge wiring), so tests drive kernels.ln_pow.prove and
exact_check.ln_pow.check directly; the dispatch and full-pipeline tests are
skip-guarded until then. The power slot carries q itself (constant is
(ln q)^4 — mirror ln_q_cube), and the quartic P is sign-tested by the
Sturm odd-part rule ln_pow.quartic_nonneg.
"""

from fractions import Fraction

import pytest
from mpmath import floor, log, mp, mpf

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import ln_pow as check_mod
from attention_calculator.exact_check import verify
from attention_calculator.kernels import EXACT_TYPES, ln_pow

KIND = "ln_q_quad"


def run(q, comp, bound):
    """prove + check: every emitted proof must verify exactly."""
    resp = ln_pow.prove(KIND, Fraction(q), comp, Fraction(bound))
    res = check_mod.check(KIND, Fraction(q), comp, Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
    return resp


def grid_bound(q: Fraction, comp: str, digits: int) -> Fraction:
    """Tight decimal bound: ``digits``-place floor for '>', ceil for '<'."""
    with mp.workdps(80):
        c = log(mpf(q.numerator) / q.denominator) ** 4
        v = int(floor(c * 10**digits))
    out = Fraction(v, 10**digits)
    return out + Fraction(1, 10**digits) if comp == "<" else out


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("2", "<", "1/4"),  # ln^4 2 = 0.2308350...
        ("2", ">", "23/100"),
        ("3", "<", "3/2"),  # ln^4 3 = 1.456725...
        ("3", ">", "29/20"),
        ("4", "<", "37/10"),  # ln^4 4 = 16 ln^4 2 = 3.6933...
        ("4", ">", "7/2"),
        ("5", "<", "671/100"),  # ln^4 5 = 6.7096...
        ("5", ">", "67/10"),
        ("6", ">", "10"),  # ln^4 6 = 10.306...
        ("3/2", ">", "27/1000"),  # ln^4(3/2) = 0.0270280...
        ("3/2", "<", "7/250"),
        ("4/3", ">", "3/500"),  # ln^4(4/3) = 0.0068494...
        ("4/3", "<", "7/1000"),
        ("5/4", ">", "1/500"),  # ln^4(5/4) = 0.0024793...
        ("5/4", "<", "1/400"),
        ("8/7", ">", "3/10000"),  # ln^4(8/7) = 0.00031817...
        ("8/7", "<", "1/2500"),
        ("9/8", ">", "1/10000"),  # ln^4(9/8) = 0.0001924..., (0,0) proof
        ("9/8", "<", "1/5000"),
        ("16/9", ">", "1/10"),  # ln^4(16/9) = 0.10959...
        ("16/9", "<", "11/100"),
    ],
)
def test_proofs_verify(q, comp, bound):
    run(q, comp, bound)


@pytest.mark.parametrize(
    ("q", "comp", "digits"),
    [
        ("2", "<", 12),
        ("2", ">", 13),  # deepest in-budget hit: (10,10)
        ("3", "<", 6),
        ("3", ">", 7),
        ("3/2", "<", 12),
        ("5/4", ">", 12),
        ("5", "<", 4),
        ("5", ">", 3),
    ],
)
def test_tight_bounds_verify(q, comp, digits):
    run(q, comp, grid_bound(Fraction(q), comp, digits))


def test_beyond_budget_reports_no_solution():
    # 14-digit lower bound on ln^4 2: true by ~1e-14 but nothing in the
    # (10,10) budget proves it — honest NoSolution, not a fake proof.
    with pytest.raises(NoSolution):
        ln_pow.prove(KIND, Fraction(2), ">", grid_bound(Fraction(2), ">", 14))
    # q = 10: ln^4 10 = 28.1101... — 5623/200 = 28.115 is a true claim but
    # every (m,n) <= (10,10) solves to a sign-indefinite P: NoSolution.
    with pytest.raises(NoSolution):
        ln_pow.prove(KIND, Fraction(10), "<", Fraction(5623, 200))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("2", ">", "1/4"),  # ln^4 2 < 1/4: opposite holds
        ("2", "<", "2/10"),
        ("3", "<", "7/5"),
        ("3", ">", "3/2"),
        ("3/2", ">", "3/100"),
        ("4/3", "<", "1/200"),
    ],
)
def test_false_claim_rejected(q, comp, bound):
    """False claims die as WrongDirection: the solved uniformly non-positive
    P certifies the opposite inequality."""
    with pytest.raises(WrongDirection):
        ln_pow.prove(KIND, Fraction(q), comp, Fraction(bound))


def test_corrupted_params_fail_identity():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 4))
    params = dict(resp["parameters"], eu_val=str(int(resp["parameters"]["eu_val"]) + 1))
    res = check_mod.check(KIND, Fraction(2), "<", Fraction(1, 4), params)
    assert not res["identity_ok"]


def test_params_against_wrong_direction_fail():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 4))
    res = check_mod.check(KIND, Fraction(2), ">", Fraction(1, 4), resp["parameters"])
    assert not res["identity_ok"]


def test_sign_indefinite_params_rejected():
    """A sign-indefinite P is refused even when its moments were honestly
    reported: P = 1 - 17x^2(1-x)^2 has P(0) = P(1) = 1 but dips below zero
    inside (0,1) — the Sturm odd-part count must catch it."""
    params = {
        "m": 0,
        "n": 0,
        "au_val": "1",
        "bu_val": "0",
        "cu_val": "-17",
        "du_val": "34",
        "eu_val": "-17",
        "u_val": "1",
    }
    res = check_mod.check(KIND, Fraction(2), "<", Fraction(1, 4), params)
    assert not res["nonneg"]


def test_input_validation():
    with pytest.raises(ValueError, match="右侧有理数格式无效"):
        ln_pow.prove(KIND, Fraction(2), "<", Fraction(-1))
    with pytest.raises(ValueError, match="请在ln后输入一个大于1的数"):
        ln_pow.prove(KIND, Fraction(1), ">", Fraction(0))
    with pytest.raises(ValueError, match="请在ln后输入一个大于1的数"):
        ln_pow.prove(KIND, Fraction(1, 2), ">", Fraction(0))


def test_render_equation_smoke():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 4))
    tex = ln_pow.render_equation(resp["parameters"], KIND, Fraction(2), "<", Fraction(1, 4))
    assert tex.startswith("\\dfrac{1}{4} - \\ln^{4}2 = \\int_0^1")
    assert "\\ln" in tex and tex.endswith("> 0")


def test_verify_dispatch(monkeypatch):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, KIND, "ln_pow")
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 4))
    res = verify(KIND, Fraction(2), "<", Fraction(1, 4), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


def test_full_pipeline_when_registered():
    """Once kernels.EXACT_TYPES / solve.FAMILY / integrand.constant_mpf wire
    the type, solver.prove(exact=True) must emit self-checking proofs."""
    if KIND not in solver.FAMILY or KIND not in EXACT_TYPES:
        pytest.skip(f"{KIND} not yet registered (merge wiring is the leader's)")
    resp = solver.prove(KIND, "2", "<", "1/4", exact=True)
    res = check_mod.check(KIND, Fraction(2), "<", Fraction(1, 4), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("coeffs", "expected"),
    [
        ((1, 0, 0, 0, 0), True),  # constant
        ((0, 0, 0, 0, -1), False),  # -x^4 <= 0
        ((0, 0, 0, 0, 0), True),  # P ≡ 0
        # 1 - 17x^2(1-x)^2: P(0) = P(1) = 1, dips to -1/16 at x = 1/2 —
        # two interior simple roots; not caught by endpoint values
        ((1, 0, -17, 34, -17), False),
        # (2x-1)^2(x^2+2)/4: interior double root at 1/2 touches, no dip
        ((Fraction(1, 2), -2, Fraction(9, 4), -1, 1), True),
        # x^3(1-x): odd-multiplicity roots at both endpoints only
        ((0, 0, 0, 1, -1), True),
        # x(1-x)(x^2+1): endpoint roots, positive inside
        ((0, 1, -1, 1, -1), True),
        # x(1-x)(2x-1)^2: endpoint roots + interior double root
        ((0, 1, -5, 8, -4), True),
        # (x^2-1)^2 = (x-1)^2(x+1)^2: even roots at both endpoints
        ((1, 0, -2, 0, 1), True),
        # -x^3-x^2/2-x/2+1/2 = (x-1/2)(-(x^2+x+1)): simple root at 1/2
        ((Fraction(1, 2), Fraction(-1, 2), Fraction(-1, 2), -1, 0), False),
        # (x-1/4)(x-1/2): two interior simple roots — dips then returns
        ((Fraction(1, 8), Fraction(-3, 4), 1, 0, 0), False),
        # x^3-x+1/4 + 0x^4 stays a cubic with interior dip (mixed degree)
        ((Fraction(1, 4), -1, 0, 1, 0), False),
        # (x-1/2)^4: quadruple interior root — still nonneg
        ((Fraction(1, 16), Fraction(-1, 2), Fraction(3, 2), -2, 1), True),
        # (1/2-x)^3(x+1): triple interior root — sign changes at 1/2
        ((Fraction(1, 8), Fraction(-5, 8), Fraction(3, 4), Fraction(1, 2), -1), False),
        # degenerate to lower degree: x^2-2x+1 = (x-1)^2
        ((1, -2, 1, 0, 0), True),
    ],
)
def test_quartic_nonneg(coeffs, expected):
    assert ln_pow.quartic_nonneg([Fraction(v) for v in coeffs]) is expected
