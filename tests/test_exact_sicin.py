"""Exact-mode tests for the si_q / cin_q family (W3 types, exact-only).

kernels/sicin.py proves Si(q) ⋚ r and Cin(q) ⋚ r over the 4-symbol span
{C_q, sin_q, cos_q, 1} with cubic P, the (1+x)^m(1-x)^n basis, and the
Taylor-remainder ladder kernels indexed by t — derivations and coverage
measurements in docs/2026-09-16-w3-research-sicin.md. Constants are
evaluated by mpmath: mp.si directly; Cin(q) = γ + ln|q| − Ci(|q|).

Structural note: the F(0) pinning forces P(0) = +1 on every solution
(the C_q row of the moment system is ±δ_{j0}), so WrongDirection is
unreachable — false claims exhaust honestly as NoSolution.
"""

from fractions import Fraction

import pytest
from mpmath import mp

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution
from attention_calculator.exact_check import verify
from attention_calculator.exact_check.sicin import check
from attention_calculator.kernels.sicin import prove, render_equation


@pytest.fixture(autouse=True, scope="module")
def module_dps():
    """本模块数值校验的 mpmath 精度。"""
    with mp.workdps(60):
        yield


KINDS = ("si_q", "cin_q")


def const_mpf(kind: str, q: Fraction):
    """C(q): mp.si for si_q; Cin(q) = γ + ln|q| − Ci(|q|) for cin_q."""
    mqf = mp.mpf(q.numerator) / q.denominator
    if kind == "si_q":
        return mp.si(mqf)
    a = abs(mqf)
    return mp.euler + mp.log(a) - mp.ci(a)


def tight_bound(kind: str, q: Fraction, comp: str, digs: int) -> Fraction:
    """``digs``-place floor of C(q) for '>', ceil for '<'."""
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
        ("si_q", "1/4"),
        ("si_q", "1/2"),
        ("si_q", "1"),
        ("si_q", "2"),
        ("cin_q", "1/4"),
        ("cin_q", "1/2"),
        ("cin_q", "1"),
        ("cin_q", "2"),
    ],
)
@pytest.mark.parametrize("comp", [">", "<"])
def test_tight_bounds_both_directions(kind, q, comp):
    """1e-3 floor/ceil of C(q) proves in both directions."""
    bound = tight_bound(kind, Fraction(q), comp, 3)
    assert_verified(kind, q, comp, bound)


@pytest.mark.parametrize(
    ("kind", "q"),
    [("si_q", "1"), ("si_q", "1/2"), ("cin_q", "1"), ("cin_q", "1/2")],
)
@pytest.mark.parametrize("comp", [">", "<"])
def test_microtight_bounds(kind, q, comp):
    """1e-6 convergent bounds prove via ladder levels t >= 1."""
    bound = tight_bound(kind, Fraction(q), comp, 6)
    assert_verified(kind, q, comp, bound)


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        # research note's verified identities; the first two land at (0,0,t=0)
        ("si_q", "1", ">", "9/10"),
        ("si_q", "1", "<", "1"),
        ("si_q", "7/2", "<", "19/10"),  # q > π '<' instance
        ("si_q", "7/2", ">", "183/100"),  # q > π '>' still provable
        ("cin_q", "1", "<", "1/4"),
        ("cin_q", "1", ">", "23/100"),
        # deep-oscillation regime: '>' fails honestly, '<' still reaches
        ("si_q", "20", "<", "9"),
        ("cin_q", "10", ">", "28/10"),
        ("cin_q", "10", "<", "3"),
    ],
)
def test_named_bounds_verify(kind, q, comp, bound):
    assert_verified(kind, q, comp, Fraction(bound))


def test_doc_identity_si_at_origin():
    """Si(1) - 9/10 = int (1+48x-32x^2+16x^3)(sin x - x + x^3/6)/x dx."""
    resp = assert_verified("si_q", "1", ">", Fraction(9, 10))
    p = resp["parameters"]
    assert (p["m"], p["n"], p["t_val"]) == (0, 0, "0")
    assert resp["solution"] == "a = 1, b = 48, c= -32, d= 16"


def test_doc_identity_cin_at_origin():
    """1/4 - Cin(1) = int (x^2/2 - (1-cos x))/x dx: the unit cubic P = 1."""
    resp = assert_verified("cin_q", "1", "<", Fraction(1, 4))
    p = resp["parameters"]
    assert (p["m"], p["n"], p["t_val"]) == (0, 0, "0")
    assert resp["solution"] == "a = 1, b = 0, c= 0, d= 0"


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("si_q", "1", ">", "1"),  # Si(1) ~= 0.9461 < 1
        ("si_q", "1", "<", "9/10"),
        ("si_q", "-1", ">", "-9/10"),  # Si(-1) ~= -0.9461
        ("cin_q", "1", "<", "23/100"),  # Cin(1) ~= 0.2398
        ("cin_q", "1", ">", "1/4"),
        ("cin_q", "2", ">", "8/5"),  # Cin(2) ~= 0.8474
    ],
)
def test_false_claim_is_no_solution(kind, q, comp, bound):
    """P(0) = +1 is forced by the C_q row, so WrongDirection cannot fire —
    a false claim exhausts the search as honest NoSolution."""
    with pytest.raises(NoSolution):
        prove(kind, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("si_q", "20", ">", "1"),  # Si(20) ~= 1.5482: true but '>' ladder needs t ~ q/2π
    ],
)
def test_beyond_budget_is_honest_no_solution(kind, q, comp, bound):
    with pytest.raises(NoSolution):
        prove(kind, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("si_q", "-1", ">", "-1"),  # Si odd: -0.9461 > -1
        ("si_q", "-1", "<", "-9/10"),
        ("si_q", "-7/2", "<", "-9/5"),  # -1.8331 < -1.8
        ("cin_q", "-1", ">", "23/100"),  # Cin even: same claim as q = 1
        ("cin_q", "-1", "<", "1/4"),
    ],
)
def test_negative_q_reduction(kind, q, comp, bound):
    """Si is odd (folds into kernel direction), Cin is even (same |q| claim)."""
    assert_verified(kind, q, comp, Fraction(bound))


@pytest.mark.parametrize("kind", KINDS)
def test_zero_q_rejected(kind):
    """q = 0 degenerates the family (1/q in the T/U recurrences)."""
    with pytest.raises(ValueError):
        prove(kind, Fraction(0), ">", Fraction(0))


@pytest.mark.parametrize("kind", KINDS)
def test_negative_t_rejected(kind):
    """t_val < 0 would drop the polynomial tail onto the bare kernel —
    sin(qx)/x is not sign-definite for q > π — so the checker rejects it."""
    params = {
        "m": 0,
        "n": 0,
        "t_val": "-1",
        "au_val": "1",
        "bu_val": "0",
        "cu_val": "0",
        "du_val": "0",
        "u_val": "1",
    }
    with pytest.raises(ValueError):
        check(kind, Fraction(1), ">", Fraction(9, 10), params)


@pytest.mark.parametrize("kind", KINDS)
def test_zero_integrand_guard(kind):
    """All-zero params: cubic_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {
        "m": 0,
        "n": 0,
        "t_val": "0",
        "au_val": "0",
        "bu_val": "0",
        "cu_val": "0",
        "du_val": "0",
        "u_val": "1",
    }
    res = check(kind, Fraction(1), ">", Fraction(0), params)
    assert res["integrand"] == {}
    assert not res["nonneg"]
    assert not res["identity_ok"]


@pytest.mark.parametrize("kind", KINDS)
def test_corrupted_params_fail_identity(kind):
    resp = prove(kind, Fraction(1), ">", Fraction(9, 10) if kind == "si_q" else Fraction(23, 100))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    bound = Fraction(9, 10) if kind == "si_q" else Fraction(23, 100)
    res = check(kind, Fraction(1), ">", bound, params)
    assert not res["identity_ok"]


@pytest.mark.parametrize("kind", KINDS)
def test_params_against_wrong_comp_fail(kind):
    bound = Fraction(9, 10) if kind == "si_q" else Fraction(23, 100)
    resp = prove(kind, Fraction(1), ">", bound)
    res = check(kind, Fraction(1), "<", bound, resp["parameters"])
    assert not res["identity_ok"]


def test_emitted_params_shape():
    resp = prove("si_q", Fraction(1), ">", Fraction(9, 10))
    p = resp["parameters"]
    assert set(p) == {
        "m",
        "n",
        "a_val",
        "b_val",
        "c_val",
        "d_val",
        "au_val",
        "bu_val",
        "cu_val",
        "du_val",
        "u_val",
        "t_val",
        "unified_form",
    }
    assert resp["type"] == "si_q"
    assert resp["solution"].startswith("a = ")
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u
    assert Fraction(p["d_val"]) == Fraction(p["du_val"]) / u


def test_render_equation_doc_identity_si():
    """Si(1) > 9/10 renders the t=0 kernel (sin x - x + x^3/6)/x times the cubic."""
    resp = prove("si_q", Fraction(1), ">", Fraction(9, 10))
    tex = render_equation(resp["parameters"], "si_q", Fraction(1), ">", Fraction(9, 10))
    assert tex == (
        "\\mathrm{Si}\\left(1\\right) - \\dfrac{9}{10} = \\int_0^1 "
        "\\frac{\\left(\\frac{x^{3}}{6} - x + \\sin{\\left(x \\right)}\\right) "
        "\\left(16 x^{3} - 32 x^{2} + 48 x + 1\\right)}{x} \\mathrm{d} x > 0"
    )


def test_render_equation_doc_identity_cin():
    """1/4 - Cin(1) renders the complementary t=0 kernel (x^2/2 - (1-cos x))/x."""
    resp = prove("cin_q", Fraction(1), "<", Fraction(1, 4))
    tex = render_equation(resp["parameters"], "cin_q", Fraction(1), "<", Fraction(1, 4))
    assert tex == (
        "\\dfrac{1}{4} - \\mathrm{Cin}\\left(1\\right) = \\int_0^1 "
        "\\frac{\\frac{x^{2}}{2} + \\cos{\\left(x \\right)} - 1}{x} \\mathrm{d} x > 0"
    )


def test_render_equation_negative_q_uses_abs():
    """The kernel prints on |q| — the oddness fold lives in the solved target."""
    resp = prove("si_q", Fraction(-1), ">", Fraction(-1))
    tex = render_equation(resp["parameters"], "si_q", Fraction(-1), ">", Fraction(-1))
    assert "\\mathrm{Si}\\left(-1\\right)" in tex
    assert "\\sin{\\left(x \\right)}" in tex  # |q| = 1, not sin(-x)


@pytest.mark.parametrize("kind", KINDS)
def test_verify_dispatch(monkeypatch, kind):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, kind, "sicin")
    bound = "9/10" if kind == "si_q" else "23/100"
    resp = prove(kind, Fraction(1), ">", Fraction(bound))
    res = verify(kind, Fraction(1), ">", Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize("kind", KINDS)
def test_full_pipeline(kind):
    """Through the registered family, solver.prove(exact=True) must emit
    self-checking proofs."""
    bound = "9/10" if kind == "si_q" else "23/100"
    resp = solver.prove(kind, "1", ">", bound, exact=True)
    res = check(kind, Fraction(1), ">", Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
