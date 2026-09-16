"""Exact-mode tests for psi1_q — the trigamma telescope kernel pair (W3).

kernels.trigamma proves psi'(q) vs a rational bound for q in QQ_{>0} on
[0,1]: '>' runs x^{q-1}(x-1-ln x)/(1-x), '<' runs x^{q+N-2}(1-x+x ln x)/(1-x)
with N the least integer giving q+N > 1. The moment space is
span{psi1_q, 1} and the (1-x)^n freedom is structurally dead, so emitted
proofs are always the one-dimensional (m, 0) hit under the site rule
a >= 0 and a+b >= 0 — which covers every true bound with no constant tail.

Every emitted proof is re-verified by exact_check.trigamma.check — dict-equal
moments plus the sign certificate. Constants are mp.psi(1, q) directly.
"""

from fractions import Fraction

import mpmath as mp
import pytest

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify
from attention_calculator.exact_check.trigamma import check
from attention_calculator.kernels import EXACT_TYPES
from attention_calculator.kernels.trigamma import prove, render_equation

mp.mp.dps = 60
KIND = "psi1_q"


def const_mpf(q: Fraction):
    """psi'(q) — the trigamma value mp.psi(1, q)."""
    return mp.psi(1, mp.mpf(q.numerator) / q.denominator)


def tight_bound(q: Fraction, comp: str, digs: int) -> Fraction:
    """``digs``-place floor of psi'(q) for '>', ceil for '<'."""
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
    "q",
    ["1/10", "1/3", "1/2", "2/3", "1", "3/2", "2", "5/2", "7"],
)
@pytest.mark.parametrize("comp", [">", "<"])
def test_tight_bounds_both_directions(q, comp):
    """1e-3 floor/ceil of psi'(q) proves in both directions."""
    bound = tight_bound(Fraction(q), comp, 3)
    assert_verified(q, comp, bound)


@pytest.mark.parametrize(
    ("q", "comp"),
    [("1/3", ">"), ("1/2", "<"), ("2/3", ">"), ("1", ">"), ("1", "<"), ("2", "<")],
)
def test_microtight_bounds(q, comp):
    """1e-4 convergent bounds still prove within the m <= 256 budget
    (the skipped mirror cases sit at delta ~ 2e-6, needing m ~ 450+)."""
    bound = tight_bound(Fraction(q), comp, 4)
    assert_verified(q, comp, bound)


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        # research-doc verified identities: psi'(1) > 8/5 at m=2,
        # psi'(1) < 7/4 at m=0 (P=x), psi'(1/3) < 21/2 at m=0 (a=1/9)
        ("1", ">", "8/5"),
        ("1", ">", "164/100"),
        ("1", "<", "7/4"),
        ("1", "<", "33/20"),
        ("1/3", ">", "10"),
        ("1/3", "<", "21/2"),
        ("1/3", "<", "11"),
        ("1/2", ">", "49/10"),
        ("1/2", "<", "5"),
        # q > 1 exercises the '<' kernel's N=0 branch (exponent q-2)
        ("2", ">", "3/5"),
        ("2", "<", "13/20"),
        ("5/2", "<", "1/2"),
        # small q: '>' loose side starts at 1/q, '<' exponent is q-1
        ("1/10", ">", "100"),
        ("1/10", "<", "102"),
    ],
)
def test_named_bounds_verify(q, comp, bound):
    assert_verified(q, comp, Fraction(bound))


def test_doc_identity_lands_at_2_0():
    """psi'(1) - 8/5 = int x^2(2+3x)/5 * (x-1-ln x)/(1-x) dx (m=2)."""
    resp = assert_verified("1", ">", Fraction(8, 5))
    assert (resp["parameters"]["m"], resp["parameters"]["n"]) == (2, 0)
    assert resp["solution"] == "a = 2/5, b = 3/5"


def test_doc_identity_lt_is_p_eq_x():
    """7/4 - psi'(1) = int x(1-x+x ln x)/(1-x) dx — the doc's P=x hit."""
    resp = assert_verified("1", "<", Fraction(7, 4))
    assert (resp["parameters"]["m"], resp["parameters"]["n"]) == (0, 0)
    assert resp["solution"] == "a = 0, b = 1"


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("1", ">", "1/2"),  # r < 1/q: loose side proves at m=0 via b < 0
        ("1", ">", "0"),  # bound 0 exercises the zero-filter
        ("1", ">", "-5"),
        ("1", "<", "2"),  # r = U_0 exactly: a=1, b=0
        ("1", "<", "100"),  # r > U_0: loose side proves at m=0 via b < 0
        ("1/10", ">", "1"),  # deep loose side, 1/q = 10
    ],
)
def test_loose_side_needs_no_tail(q, comp, bound):
    """Under the site rule (a>=0, a+b>=0; a+b=1 always) the linear P already
    covers the doc's constant-tail ranges — no tail encoding exists."""
    resp = assert_verified(q, comp, Fraction(bound))
    assert resp["parameters"]["m"] == 0


@pytest.mark.parametrize(
    ("q", "comp", "bound"),
    [
        ("1", ">", "7/4"),  # psi'(1) ~= 1.6449
        ("1", "<", "8/5"),
        ("1/3", ">", "11"),  # psi'(1/3) ~= 10.0955
        ("1/3", "<", "10"),
        ("1/2", "<", "49/10"),  # psi'(1/2) = pi^2/2 ~= 4.9348
        ("2", "<", "3/5"),  # psi'(2) = pi^2/6 - 1 ~= 0.6449
    ],
)
def test_false_claim_is_honest_no_solution(q, comp, bound):
    """False claims solve P with a+b=1 and a<0: never nonneg, never nonpos —
    the m scan exhausts to NoSolution (solver.prove reports WrongDirection
    earlier, via certified_cmp)."""
    with pytest.raises(NoSolution):
        prove(KIND, Fraction(q), comp, Fraction(bound))


@pytest.mark.parametrize("q", ["0", "-1", "-1/2", "-7/3"])
def test_nonpositive_q_rejected(q):
    """psi' poles at 0,-1,-2,...; negative non-integers are out of scope."""
    with pytest.raises(ValueError):
        prove(KIND, Fraction(q), ">", Fraction(0))
    with pytest.raises(ValueError):
        prove(KIND, Fraction(q), "<", Fraction(1))


def test_beyond_budget_is_honest_no_solution():
    """psi'(1) vs its 1e-6 floor needs m ~ 700 — past the 256 budget."""
    bound = tight_bound(Fraction(1), ">", 6)
    with pytest.raises(NoSolution):
        prove(KIND, Fraction(1), ">", bound)


def test_zero_integrand_guard():
    """All-zero params: poly_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {"m": 0, "n": 0, "au_val": "0", "bu_val": "0", "u_val": "1"}
    res = check(KIND, Fraction(1), ">", Fraction(0), params)
    assert res["integrand"] == {}
    assert not res["nonneg"]
    assert not res["identity_ok"]


def test_forged_n_fails_identity():
    """A hand-forged n>=1 kills the psi' row (the (1-x)^n moments are purely
    rational) — dict equality fails honestly rather than matching."""
    resp = prove(KIND, Fraction(1), ">", Fraction(8, 5))
    params = dict(resp["parameters"], n=1)
    res = check(KIND, Fraction(1), ">", Fraction(8, 5), params)
    assert not res["identity_ok"]


def test_corrupted_params_fail_identity():
    resp = prove(KIND, Fraction(1), ">", Fraction(8, 5))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check(KIND, Fraction(1), ">", Fraction(8, 5), params)
    assert not res["identity_ok"]


def test_params_against_wrong_comp_fail():
    """'< ' params checked under '>' rebuild the other kernel's moments."""
    resp = prove(KIND, Fraction(1), "<", Fraction(7, 4))
    res = check(KIND, Fraction(1), ">", Fraction(7, 4), resp["parameters"])
    assert not res["identity_ok"]


def test_emitted_params_shape():
    resp = prove(KIND, Fraction(1), ">", Fraction(8, 5))
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
    assert resp["type"] == KIND
    assert resp["solution"].startswith("a = ")
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u
    assert Fraction(p["b_val"]) == Fraction(p["bu_val"]) / u


def test_render_equation_gt():
    """psi'(1) > 8/5 prints the m=2 integrand over the '>' kernel."""
    resp = prove(KIND, Fraction(1), ">", Fraction(8, 5))
    tex = render_equation(resp["parameters"], KIND, Fraction(1), ">", Fraction(8, 5))
    assert tex == (
        "\\psi_1\\left(1\\right) - \\dfrac{8}{5} = \\int_0^1 "
        "\\frac{x^{2} \\left(3 x + 2\\right) \\left(x - \\log{\\left(x \\right)} - 1\\right)}"
        "{5 \\left(1 - x\\right)} \\mathrm{d} x > 0"
    )


def test_render_equation_lt_third():
    """psi'(1/3) < 21/2 prints the singular x^{-2/3} '<' kernel factor."""
    resp = prove(KIND, Fraction(1, 3), "<", Fraction(21, 2))
    tex = render_equation(resp["parameters"], KIND, Fraction(1, 3), "<", Fraction(21, 2))
    assert tex == (
        "\\dfrac{21}{2} - \\psi_1\\left(\\dfrac{1}{3}\\right) = \\int_0^1 "
        "\\frac{\\left(8 x + 1\\right) \\left(x \\log{\\left(x \\right)} - x + 1\\right)}"
        "{9 x^{\\frac{2}{3}} \\left(1 - x\\right)} \\mathrm{d} x > 0"
    )


def test_verify_dispatch(monkeypatch):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, KIND, "trigamma")
    resp = prove(KIND, Fraction(1), ">", Fraction(8, 5))
    res = verify(KIND, Fraction(1), ">", Fraction(8, 5), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


def test_full_pipeline_when_registered():
    """Once kernels.EXACT_TYPES / solve.FAMILY / integrand.constant_mpf wire
    the type, solver.prove(exact=True) must emit self-checking proofs."""
    if KIND not in solver.FAMILY or KIND not in EXACT_TYPES:
        pytest.skip(f"{KIND} not yet registered (merge wiring is the leader's)")
    resp = solver.prove(KIND, "1", ">", "8/5", exact=True)
    res = check(KIND, Fraction(1), ">", Fraction(8, 5), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


def test_full_pipeline_false_claim_when_registered():
    """Through solve.prove(exact=True) a false claim is WrongDirection —
    certified_cmp decides before the kernel scan."""
    if KIND not in solver.FAMILY or KIND not in EXACT_TYPES:
        pytest.skip(f"{KIND} not yet registered (merge wiring is the leader's)")
    with pytest.raises(WrongDirection):
        solver.prove(KIND, "1", ">", "2", exact=True)
