"""Exact-mode tests for the arcsin_q kernel (W3 type extension, exact-only).

Every emitted proof is re-verified by exact_check.arcsin.check -- dict-equal
moments plus the sign certificate. q with 1 - q^2 a rational square (3/5,
4/5, 5/13, 7/25, ...) takes the folded 2-dim path with linear P; other q the
symbolic-w 3-dim path with quadratic P. Coverage degrades for q -> 1: the
provable window only approaches arcsin q through strongly asymmetric (m, n),
so near-1 q at tight bounds honestly report NoSolution.
"""

from fractions import Fraction

import pytest
from mpmath import mp

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check.arcsin import check
from attention_calculator.kernels.arcsin import prove


@pytest.fixture(autouse=True, scope="module")
def module_dps():
    """本模块数值校验的 mpmath 精度。"""
    with mp.workdps(60):
        yield


def tight_bounds(q: Fraction, digs: int) -> tuple[Fraction, Fraction]:
    """floor/ceil of arcsin(q) at 10^-digs: provable '>' / '<' bounds."""
    c = mp.asin(mp.mpf(q.numerator) / q.denominator)
    d = 10**digs
    return Fraction(int(mp.floor(c * d)), d), Fraction(int(mp.ceil(c * d)), d)


def assert_verified(q: str, comp: str, bound: str) -> dict:
    """prove() must emit and check() must confirm the exact identity."""
    resp = prove("arcsin_q", Fraction(q), comp, Fraction(bound))
    res = check("arcsin_q", Fraction(q), comp, Fraction(bound), resp["parameters"])
    assert res["identity_ok"], (q, comp, bound, res)
    assert res["nonneg"], (q, comp, bound, res)
    return resp


@pytest.mark.parametrize("q", ["1/10", "1/5", "1/4", "1/3", "2/7", "1/2", "9/16"])
@pytest.mark.parametrize("digs", [3, 6])
def test_tight_bounds_both_directions(q, digs):
    """Irrational-w q: convergent-tight bounds in both directions."""
    lo, hi = tight_bounds(Fraction(q), digs)
    assert_verified(q, ">", str(lo))
    assert_verified(q, "<", str(hi))


@pytest.mark.parametrize("q", ["3/5", "4/5", "5/13", "7/25", "12/13", "8/17"])
@pytest.mark.parametrize("digs", [4, 6])
def test_perfect_square_q_folds_to_linear_p(q, digs):
    """w rational: 2-dim moment space -> linear P (cu_val stays '0')."""
    lo, hi = tight_bounds(Fraction(q), digs)
    for comp, bound in ((">", lo), ("<", hi)):
        resp = assert_verified(q, comp, str(bound))
        assert resp["parameters"]["cu_val"] == "0"


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("9/10", ">", "111/100"),
        ("9/10", "<", "9/8"),
        ("17/20", "<", "26/25"),
        ("99/100", ">", "6/5"),
        ("99/100", ">", "19/16"),
    ],
)
def test_near_one_q_loose_bounds(q, comp, bound):
    """q -> 1 still proves loose bounds via asymmetric (m, n)."""
    assert_verified(q, comp, bound)


def test_near_one_q_tight_bound_is_honest_no_solution():
    """asin(0.99) < 36/25 is true but unprovable within budget -- the
    construction's window closes toward q = 1; NoSolution is honest."""
    with pytest.raises(NoSolution):
        prove("arcsin_q", Fraction(99, 100), "<", Fraction(36, 25))


@pytest.mark.parametrize(
    ("comp", "bound"),
    [
        (">", "3/5"),  # asin(1/2) ~= 0.5236 < 3/5
        ("<", "1/2"),  # asin(1/2) ~= 0.5236 > 1/2
        (">", "1"),
        ("<", "0"),
    ],
)
def test_false_claim_rejected(comp, bound):
    """False claims never emit: nonpos solves give WrongDirection, fully
    indefinite scans give NoSolution -- both are honest rejections."""
    with pytest.raises((WrongDirection, NoSolution)):
        prove("arcsin_q", Fraction(1, 2), comp, Fraction(bound))


@pytest.mark.parametrize("q", ["1/2", "3/5"])
def test_zero_bound_claim(q):
    """bound = 0 drops the "1" key from target (Moments omit zero coeffs);
    arcsin(q) > 0 must still emit and verify on both search paths."""
    assert_verified(q, ">", "0")


@pytest.mark.parametrize("q", ["0", "1", "3/2", "-1/2"])
def test_out_of_domain_rejected(q):
    with pytest.raises(ValueError):
        prove("arcsin_q", Fraction(q), ">", Fraction(0))


def test_emitted_solution_string_shape():
    resp = prove("arcsin_q", Fraction(1, 3), "<", Fraction(3, 5))
    assert resp["solution"].startswith("a = ")
    assert "b = " in resp["solution"]


def test_render_equation_matches_author_example():
    """The author's published (1,1) proof of asin(1/3) < 3/5 renders his
    single-fraction form: numerator (au,bu,cu)/u, denominator u/d * sqrt."""
    from attention_calculator.kernels.arcsin import render_equation

    params = {
        "m": 1,
        "n": 1,
        "au_val": "1710",
        "bu_val": "37",
        "cu_val": "-236",
        "u_val": "1080",
        "a_val": "19/12",
        "b_val": "37/1080",
        "c_val": "-59/270",
    }
    out = render_equation(params, "arcsin_q", Fraction(1, 3), "<", Fraction(3, 5))
    assert out == (
        "\\dfrac{3}{5} - \\arcsin\\left(\\dfrac{1}{3}\\right)"
        " = \\int_0^1 \\dfrac{x \\left(1 - x\\right) \\left(- 236 x^{2} + 37 x + 1710\\right)}"
        "{360 \\sqrt{9 - x^{2}}} \\mathrm{d} x > 0"
    )


# 11/21 > asin(1/2)=0.523598… would be a false claim; the '>' bound must sit
# just below the true value (261799/500000 = 0.523598 < 0.5235987…)
@pytest.mark.parametrize(("comp", "bound"), [(">", "261799/500000"), ("<", "11/20")])
def test_full_pipeline(comp, bound):
    """Through the registered type, the certified path must emit only
    self-checking proofs."""
    resp = solver.prove("arcsin_q", "1/2", comp, bound, exact=True)
    res = check("arcsin_q", Fraction(1, 2), comp, Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
