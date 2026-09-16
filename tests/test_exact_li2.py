"""Exact-mode tests for li2_q (W3 type, exact-only).

The kernel proves ``Li_2(q) ~ bound`` for rational q < 1, q != 0 via
int_0^1 (1-x)^n P(x) sigma(phi - L) dx >= 0, phi = -ln(1-qx)/x.  Moments live
in span{li2_q, ln_1mq, 1}; the Li_2 symbol sits only in A_0, so m == 0 and
a == +1 in every emitted proof.  Every emitted proof is re-verified by
exact_check.li2.check — dict-equal moments plus the family sign lemma
(kernels.li2.certified).

Coverage: 1e-6 floors/ceilings prove for |q| <= 9/10 in both directions
(d axis reaches ~57 at q = 9/10).  q < -1 has only the linear chord/tan0
witnesses (the series diverges at x > 1/|q|), so tight bounds there are
honestly NoSolution — the alternating families cover -1 <= q < 0.
"""

from fractions import Fraction

import pytest
from mpmath import mp

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify
from attention_calculator.exact_check.li2 import check
from attention_calculator.kernels.li2 import prove, render_equation

KIND = "li2_q"


def const_mpf(q: Fraction):
    with mp.workdps(60):
        return mp.polylog(2, mp.mpf(q.numerator) / q.denominator)


def tight_bound(q: Fraction, comp: str, digs: int) -> Fraction:
    """``digs``-place floor of Li_2(q) for '>', ceil for '<'."""
    c = const_mpf(q)
    v = int(mp.floor(c * 10**digs))
    out = Fraction(v, 10**digs)
    return out + Fraction(1, 10**digs) if comp == "<" else out


def assert_verified(q: str, comp: str, bound: Fraction) -> dict:
    """prove() must emit and check() must confirm the exact identity."""
    resp = prove(KIND, Fraction(q), comp, bound)
    res = check(KIND, Fraction(q), comp, bound, resp["parameters"])
    assert res["identity_ok"], (q, comp, bound, res)
    assert res["nonneg"], (q, comp, bound, res)
    return resp


@pytest.mark.parametrize(
    ("q", "comp", "bound", "fam", "d", "n"),
    [
        # the research note's verified identities
        ("1/2", ">", "4/9", "taylor", 0, 0),
        ("1/2", ">", "29/50", "taylor", 2, 0),
        ("1/2", "<", "59/100", "chord", 0, 0),
        ("-1/2", ">", "-9/20", "alt", 3, 0),
        ("-1/2", "<", "-2/5", "alt", 2, 0),
    ],
)
def test_doc_identities(q, comp, bound, fam, d, n):
    resp = assert_verified(q, comp, Fraction(bound))
    p = resp["parameters"]
    assert (p["fam"], p["d"], int(p["m"]), int(p["n"])) == (fam, d, 0, n)
    assert p["a_val"] == "1"  # the Li_2 column pins P(0) = +1


def test_doc_identity_coefficients():
    """Li_2(1/2) - 4/9 = int (1 - 8x/3 + 16x^2/9) phi dx from the note."""
    resp = assert_verified("1/2", ">", Fraction(4, 9))
    assert resp["solution"] == "a = 1, b = -8/3, c= 16/9"


@pytest.mark.parametrize("q", ["1/3", "1/2", "2/3", "9/10", "-1/2", "-9/10", "-1"])
@pytest.mark.parametrize("comp", [">", "<"])
def test_tight_bounds_both_directions(q, comp):
    """1e-3 floor/ceil of Li_2(q) proves in both directions."""
    assert_verified(q, comp, tight_bound(Fraction(q), comp, 3))


@pytest.mark.parametrize("q", ["1/2", "2/3", "9/10", "-1/2", "-9/10", "-1"])
@pytest.mark.parametrize("comp", [">", "<"])
def test_microtight_bounds(q, comp):
    """1e-6 convergent bounds still prove within the (n, d) budget."""
    assert_verified(q, comp, tight_bound(Fraction(q), comp, 6))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("-2", "<", "-1/2"),  # q < -1: tan0 tangent witness (series diverges)
        ("-2", ">", "-3/2"),  # q < -1: chord sub-function
        ("-3", "<", "-9/5"),  # Li_2(-3) ~= -1.939
        ("-3", ">", "-2"),
    ],
)
def test_q_below_minus_one(q, comp, bound):
    """Beyond the series' radius only the concavity witnesses apply."""
    assert_verified(q, comp, Fraction(bound))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("1/2", ">", "3/5"),  # Li_2(1/2) ~= 0.5822 < 3/5
        ("1/2", "<", "7/16"),
        ("-1/2", ">", "-2/5"),  # Li_2(-1/2) ~= -0.448 < -0.4
        ("-1/2", "<", "-9/20"),
        ("9/10", ">", "13/10"),  # Li_2(0.9) ~= 1.2997
    ],
)
def test_false_claim_rejected(q, comp, bound):
    """A false claim cannot produce a nonneg P on a certified kernel —
    the scan finds only sign-indefinite solutions (NoSolution) or a
    non-positive one certifying the opposite (WrongDirection)."""
    with pytest.raises((WrongDirection, NoSolution)):
        prove(KIND, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("99/100", ">", "1588/1000"),  # 1e-3 floor: d budget exhausted (q->1^-
        ("99/100", "<", "159/100"),  # weak spot; loose '<' also dies: C' ~ 1/(1-q))
    ],
)
def test_near_one_beyond_budget(q, comp, bound):
    """True claims past the d budget report NoSolution, not a wrong proof."""
    with pytest.raises(NoSolution):
        prove(KIND, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize("q", ["0", "1", "3/2", "2"])
def test_domain_rejected(q):
    """q = 0 degenerates the kernel; q >= 1 puts the branch point in-domain."""
    with pytest.raises(ValueError):
        prove(KIND, Fraction(q), ">", Fraction(0))


def test_zero_integrand_guard():
    """All-zero params: poly_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {
        "m": 0,
        "n": 0,
        "au_val": "0",
        "bu_val": "0",
        "cu_val": "0",
        "u_val": "1",
        "fam": "taylor",
        "d": 0,
    }
    res = check(KIND, Fraction(1, 2), ">", Fraction(0), params)
    assert res["integrand"] == {}
    assert not res["nonneg"]
    assert not res["identity_ok"]


def test_corrupted_params_fail_identity():
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check(KIND, Fraction(1, 2), ">", Fraction(4, 9), params)
    assert not res["identity_ok"]


def test_corrupted_family_tag_fails_nonneg():
    """An illegal (fam, d, comp, q) pairing breaks the kernel-sign lemma:
    alt with even d for '>' claims p_d <= phi, which is the wrong side."""
    resp = prove(KIND, Fraction(-1, 2), ">", Fraction(-9, 20))
    bad = dict(resp["parameters"], d=4)  # emitted d is odd; even d is a super-fn
    res = check(KIND, Fraction(-1, 2), ">", Fraction(-9, 20), bad)
    assert not res["nonneg"]


def test_corrupted_m_fails_identity():
    """m > 0 kills the Li_2 component — the identity cannot hold."""
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
    bad = dict(resp["parameters"], m=1)
    res = check(KIND, Fraction(1, 2), ">", Fraction(4, 9), bad)
    assert not res["identity_ok"]


def test_params_against_wrong_comp_fail():
    """Under '<' the same params give a true identity (sigma flips the whole
    equation) but a vacuous one: taylor is not a certified '<' family, so the
    sign lemma rejects it."""
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
    res = check(KIND, Fraction(1, 2), "<", Fraction(4, 9), resp["parameters"])
    assert res["identity_ok"]
    assert not res["nonneg"]


def test_emitted_params_shape():
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
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
        "fam",
        "d",
    }
    assert resp["type"] == "li2_q"
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u == 1


def test_render_equation_doc_identity():
    """The bare-kernel identity renders P(x) * (-log(1-x/2)/x)."""
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
    tex = render_equation(resp["parameters"], KIND, Fraction(1, 2), ">", Fraction(4, 9))
    assert tex.startswith("\\mathrm{Li}_2\\left(\\dfrac{1}{2}\\right) - \\dfrac{4}{9} = \\int_0^1")
    assert tex.endswith("\\mathrm{d} x > 0")
    assert "16 x^{2} - 24 x + 9" in tex


def test_verify_dispatch(monkeypatch):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, KIND, "li2")
    resp = prove(KIND, Fraction(1, 2), ">", Fraction(4, 9))
    res = verify(KIND, Fraction(1, 2), ">", Fraction(4, 9), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


def test_full_pipeline():
    """Through the registered type, solver.prove(exact=True) must emit
    self-checking proofs."""
    resp = solver.prove(KIND, "1/2", ">", "4/9", exact=True)
    res = check(KIND, Fraction(1, 2), ">", Fraction(4, 9), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
    assert "certificate" in resp
