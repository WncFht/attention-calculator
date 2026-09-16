"""Exact-mode tests for ln_q_cube (kernels.ln_pow).

The type is exact-mode-only and not yet registered in kernels.TYPES /
solve.FAMILY (merge wiring), so tests drive kernels.ln_pow.prove and
exact_check.ln_pow.check directly. Once wired, direction decisions run
through solve.certified_cmp (constant_mpf entry ``mp.log(q)**3``); at
kernel level a false claim already dies as WrongDirection — the solved
uniformly non-positive P certifies the opposite inequality — or as
NoSolution, never a fake proof.
"""

from fractions import Fraction

import pytest
from mpmath import floor, log, mp, mpf

from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import ln_pow as check_mod
from attention_calculator.kernels import ln_pow

KIND = "ln_q_cube"


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
        c = log(mpf(q.numerator) / q.denominator) ** 3
        v = int(floor(c * 10**digits))
    out = Fraction(v, 10**digits)
    return out + Fraction(1, 10**digits) if comp == "<" else out


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("2", "<", "1/3"),  # ln^3 2 = 0.3330246...
        ("2", ">", "33/100"),
        ("3", "<", "133/100"),
        ("3", ">", "13/10"),
        ("3/2", ">", "3/50"),  # ln^3(3/2) = 0.0684658...
        ("3/2", "<", "7/100"),
        ("5", ">", "4"),  # ln^3 5 = 4.1699...
        ("5", "<", "21/5"),
        ("10", "<", "611/50"),  # ln^3 10 = 12.219...
        ("10", ">", "61/5"),
        ("4/3", ">", "1/50"),  # ln^3(4/3) = 0.02380...
        ("4/3", "<", "3/125"),
        ("5/4", "<", "139/12500"),  # ln^3(5/4) = 0.011122...
    ],
)
def test_proofs_verify(q, comp, bound):
    run(q, comp, bound)


@pytest.mark.parametrize(
    ("q", "comp", "digits"),
    [("2", ">", 6), ("2", "<", 6), ("2", ">", 12), ("2", "<", 12), ("3", "<", 8), ("3/2", ">", 8)],
)
def test_tight_bounds_verify(q, comp, digits):
    run(q, comp, grid_bound(Fraction(q), comp, digits))


def test_beyond_budget_reports_no_solution():
    # 15-digit bound: true by ~1e-15 but no (m,n) <= (10,10) proves it —
    # the honest answer is NoSolution, not a fake proof.
    with pytest.raises(NoSolution):
        ln_pow.prove(KIND, Fraction(2), ">", Fraction(333024651988929, 10**15))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("2", ">", "1/3"),  # ln^3 2 < 1/3: opposite holds
        ("2", "<", "3/10"),
        ("3", "<", "13/10"),
        ("3/2", ">", "7/100"),
        ("5", "<", "4"),
    ],
)
def test_false_claim_rejected(q, comp, bound):
    with pytest.raises(WrongDirection):
        ln_pow.prove(KIND, Fraction(q), comp, Fraction(bound))


def test_corrupted_params_fail_identity():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 3))
    params = dict(resp["parameters"], du_val=str(int(resp["parameters"]["du_val"]) + 1))
    res = check_mod.check(KIND, Fraction(2), "<", Fraction(1, 3), params)
    assert not res["identity_ok"]


def test_params_against_wrong_direction_fail():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 3))
    res = check_mod.check(KIND, Fraction(2), ">", Fraction(1, 3), resp["parameters"])
    assert not res["identity_ok"]


def test_input_validation():
    with pytest.raises(ValueError, match="右侧有理数格式无效"):
        ln_pow.prove(KIND, Fraction(2), "<", Fraction(-1))
    with pytest.raises(ValueError, match="请在ln后输入一个大于1的数"):
        ln_pow.prove(KIND, Fraction(1), ">", Fraction(0))
    with pytest.raises(ValueError, match="请在ln后输入一个大于1的数"):
        ln_pow.prove(KIND, Fraction(1, 2), ">", Fraction(0))


def test_render_equation_smoke():
    resp = ln_pow.prove(KIND, Fraction(2), "<", Fraction(1, 3))
    tex = ln_pow.render_equation(resp["parameters"], KIND, Fraction(2), "<", Fraction(1, 3))
    assert tex.startswith("\\dfrac{1}{3} - \\ln^{3}2 = \\int_0^1")
    assert "\\ln" in tex and tex.endswith("> 0")


@pytest.mark.parametrize(
    ("coeffs", "expected"),
    [
        ((1, 0, 0, 0), True),  # constant
        ((0, -1, 0, 1), False),  # x^3-x dips below 0 inside (0,1)
        ((1, -3, 0, 1), False),  # P(1) = -1: endpoint catches it
        # (2x-1)^2(x+1)/4: interior zero touch, square D
        ((Fraction(1, 4), Fraction(-3, 4), 0, 1), True),
        ((Fraction(1, 4), -1, 0, 1), False),  # x^3-x+1/4: interior min < 0, nonsquare D
        ((0, 0, 0, -1), False),  # -x^3 <= 0
        ((0, 1, 0, 1), True),  # x^3+x: P' has no real root
        ((0, 1, -2, 1), True),  # x(x-1)^2: interior critical point is a max
    ],
)
def test_cubic_nonneg(coeffs, expected):
    assert ln_pow.cubic_nonneg([Fraction(v) for v in coeffs]) is expected
