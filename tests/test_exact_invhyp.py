"""Exact-mode tests for the invhyp family: arsinh_q (new W3 type).

Drives kernels.invhyp.prove + exact_check.invhyp.check directly — the type
is wired into solve.FAMILY/EXACT_TYPES and has no site counterpart.
Bounds are CF convergents of arsinh(q) chosen with mpmath:

    arsinh(1/2) ≈ 0.4812118250596   arsinh(1) ≈ 0.8813735870195
    arsinh(2) ≈ 1.4436354751788     arsinh(3/5) ≈ 0.5688248987322

Every emitted proof must pass check() exactly. False claims raise
WrongDirection mid-scan: the identity int f = ±(C−r) holds exactly for every
solved candidate, so a false claim's first sign-definite P is uniformly
non-positive. A true claim beyond the (m,n) <= 10 budget exhausts to
NoSolution. Negative q is admitted (arsinh odd, kernel sees only q^2); q = 0
crashes on M_0's 1/q like sinh_q's own 1/q.
"""

from fractions import Fraction

import pytest

from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check.invhyp import check
from attention_calculator.kernels import invhyp


def run(q, comp, bound):
    resp = invhyp.prove("arsinh_q", Fraction(q), comp, Fraction(bound))
    return check("arsinh_q", Fraction(q), comp, Fraction(bound), resp["parameters"])


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        # arsinh(1/2) '<' upper bounds, loose -> 7.6e-12 tight
        ("1/2", "<", "13/27"),  # err +2.7e-4 -> (2,3)
        ("1/2", "<", "397/825"),  # err +3.0e-7 -> (4,5)
        ("1/2", "<", "2049/4258"),  # err +1.1e-8 -> (4,6)
        ("1/2", "<", "9848/20465"),  # err +7.6e-12 -> (8,6)
        # arsinh(1) '<' and '>'
        ("1", "<", "9/10"),  # err +1.9e-2 -> (1,2)
        ("1", "<", "15/17"),  # err +9.8e-4 -> (2,4)
        ("1", "<", "795/902"),  # err +1.1e-6 -> (7,7)
        ("1", ">", "4/5"),  # err -8.1e-2 -> (1,2)
        ("1", ">", "7/8"),  # err -6.4e-3 -> (2,1)
        # q spread: 3/5, 2, 99/100, 7/3, 1/10
        ("3/5", "<", "1186/2085"),  # err +4.1e-8 -> (6,6)
        ("3/5", ">", "1/2"),
        ("2", "<", "3/2"),  # -> (1,4)
        ("2", "<", "13/9"),  # err +8.1e-4 -> (3,10)
        ("2", ">", "7/5"),  # -> (2,2)
        ("99/100", "<", "9/10"),  # -> (1,2)
        ("99/100", "<", "7/8"),  # err +7.2e-4 -> (2,4)
        ("99/100", "<", "153/175"),  # err +9.0e-7 -> (7,7)
        ("7/3", "<", "8/5"),  # -> (2,8)
        ("1/10", "<", "1/10"),  # err +1.7e-4 -> (1,1)
        ("1/10", "<", "301/3015"),  # err +8.4e-8 -> (2,2)
        # negative q: arsinh odd -> mirror of the positive claims
        ("-1", "<", "-4/5"),  # arsinh(-1) < -0.8 -> (1,2)
        ("-1", ">", "-9/10"),  # -> (1,2)
        ("-1", ">", "-795/902"),  # err -1.1e-6 -> (7,7)
        # bound 0: arsinh(q) > 0 for q > 0, target carries no "1" key
        ("1", ">", "0"),
        ("2", ">", "0"),
    ],
)
def test_emitted_proofs_verify(q, comp, bound):
    res = run(q, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("1/2", "<", "64/133"),  # 64/133 > arsinh(1/2) by 8.8e-6
        ("1/2", ">", "397/825"),  # convergent is below the constant
        ("1", "<", "52/59"),  # err -1.8e-5
        ("1", ">", "795/902"),  # 795/902 < arsinh(1)
        ("3/5", ">", "33/58"),  # 33/58 > arsinh(3/5) by 1.4e-4
        ("99/100", "<", "146/167"),  # err -3.3e-5
        ("1/10", ">", "301/3015"),  # 301/3015 < arsinh(1/10)
        ("2", ">", "13/9"),  # 13/9 > arsinh(2)
        ("-1", "<", "-795/902"),  # mirror of the q=1 case above
    ],
)
def test_false_claim_raises_wrong_direction(q, comp, bound):
    # prove() does no direction screening; the mid-scan WrongDirection is the
    # exact-moment certificate that the flipped inequality holds
    with pytest.raises(WrongDirection):
        invhyp.prove("arsinh_q", Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("2", "<", "397/275"),  # true (err +8.9e-7) but needs (m,n) > 10
        ("1", "<", "50768/57601"),  # true (err +2.4e-10), beyond budget
        ("10", "<", "3"),  # true (err +1.8e-3), q large -> slow convergence
    ],
)
def test_true_claim_beyond_budget_reports_no_solution(q, comp, bound):
    with pytest.raises(NoSolution):
        invhyp.prove("arsinh_q", Fraction(q), comp, Fraction(bound))


def test_zero_argument_hits_the_moment_wall():
    # arsinh(0) = 0 is rational; M_0 = arsinh(q)/q cannot start at q = 0 —
    # the same natural crash sinh_q surfaces for its 1/q
    with pytest.raises(ZeroDivisionError):
        invhyp.prove("arsinh_q", Fraction(0), "<", Fraction(1, 2))


def test_corrupted_params_fail_identity():
    resp = invhyp.prove("arsinh_q", Fraction(1), "<", Fraction(9, 10))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check("arsinh_q", Fraction(1), "<", Fraction(9, 10), params)
    assert not res["identity_ok"]


def test_proof_does_not_verify_flipped_claim():
    resp = invhyp.prove("arsinh_q", Fraction(1), "<", Fraction(9, 10))
    res = check("arsinh_q", Fraction(1), ">", Fraction(9, 10), resp["parameters"])
    assert not res["identity_ok"]


def test_render_equation_smoke():
    resp = invhyp.prove("arsinh_q", Fraction(1), "<", Fraction(9, 10))
    tex = invhyp.render_equation(resp["parameters"], "arsinh_q", Fraction(1), "<", Fraction(9, 10))
    assert tex.startswith("\\dfrac{9}{10} - \\operatorname{arsinh}\\left(1\\right) = \\int_0^1 ")
    assert tex.endswith("\\mathrm{d} x > 0")
