"""Exact-mode tests for the beta_even family: beta4..beta10 (new W7 types).

Drives kernels.beta_even.prove + exact_check.beta_even.check directly; the
types are wired into solve.FAMILY/EXACT_TYPES (full-pipeline test below)
and have no site counterpart. Bounds are CF convergents of β(s) computed
with mpmath at 80dps (workdps, not bare mp.dps assignment):

    β(4) ≈ 0.98894455174110533611   β(6) ≈ 0.99868522221843813544
    β(8) ≈ 0.99984999024682965634   β(10) ≈ 0.99998316402619687741

Every emitted proof must pass check() exactly. False claims raise
WrongDirection mid-scan: the identity ∫f = ±(C−r) holds exactly for every
solved candidate, so a false claim's first sign-definite P is uniformly
non-positive. A true claim beyond the (m,n) ≤ 10 budget exhausts to
NoSolution (measured: β(4) > 1873596/1894541, err 7.7e-14, needs (11,11)).
"""

from fractions import Fraction

import pytest

from attention_calculator import solve as solver
from attention_calculator.certificate import verify_cert
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify
from attention_calculator.exact_check.beta_even import check
from attention_calculator.kernels import beta_even

KINDS = ("beta4", "beta6", "beta8", "beta10")


def run(kind, power, comp, bound):
    resp = beta_even.prove(kind, Fraction(power), comp, Fraction(bound))
    return check(kind, Fraction(power), comp, Fraction(bound), resp["parameters"])


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        # β(4) '<' upper bounds, loose -> 5.2e-12 tight
        ("beta4", "1", "<", "2"),  # (0,0)
        ("beta4", "1", "<", "1"),  # (0,0)
        ("beta4", "1", "<", "98309/99408"),  # err -5.2e-12 -> (9,9)
        # β(4) '>' lower bounds
        ("beta4", "1", ">", "0"),
        ("beta4", "1", ">", "1/2"),
        ("beta4", "1", ">", "9/10"),  # (1,0)
        ("beta4", "1", ">", "89/90"),  # err -5.6e-5 -> (3,1)
        ("beta4", "1", ">", "805/814"),  # err -1.1e-6 -> (4,4)
        ("beta4", "1", ">", "5725/5789"),  # err -1.7e-9 -> (6,7)
        # β(6) '<' upper bounds
        ("beta6", "1", "<", "2"),
        ("beta6", "1", "<", "1"),
        ("beta6", "1", "<", "559814/560551"),  # err -4.1e-13 -> (8,9)
        # β(6) '>' lower bounds
        ("beta6", "1", ">", "0"),
        ("beta6", "1", ">", "9/10"),  # (1,0)
        ("beta6", "1", ">", "129889/130060"),  # err -1.3e-11 -> (7,7)
        ("beta6", "1", ">", "4048587/4053917"),  # err -2.6e-14 -> (9,10)
        # β(8) '<' upper bounds
        ("beta8", "1", "<", "2"),
        ("beta8", "1", "<", "1"),
        ("beta8", "1", "<", "199957/199987"),  # err +2.5e-12 -> (6,5)
        # β(8) '>' lower bounds
        ("beta8", "1", ">", "0"),
        ("beta8", "1", ">", "9/10"),
        ("beta8", "1", ">", "6665/6666"),  # err -5.2e-9 -> (3,2)
        ("beta8", "1", ">", "86648/86661"),  # err -5.5e-11 -> (5,3)
        # β(10) '<' upper bounds — β(10) hugs 1, so even 3.8e-13 gaps are easy
        ("beta10", "1", "<", "2"),
        ("beta10", "1", "<", "59396/59397"),  # err +1.1e-10 -> (4,0)
        ("beta10", "1", "<", "475165/475173"),  # err +3.8e-13 -> (6,2)
        # β(10) '>' lower bounds
        ("beta10", "1", ">", "0"),
        ("beta10", "1", ">", "9/10"),
        ("beta10", "1", ">", "49999/50000"),  # err -3.2e-6 -> (1,0)
        ("beta10", "1", ">", "999983/1000000"),  # err -1.6e-7 -> (1,0)
        # power as coefficient: 2β(4) ≈ 1.9778891, β(4)/2 ≈ 0.4944723
        ("beta4", "2", ">", "1"),
        ("beta4", "2", "<", "2"),
        ("beta4", "1/2", ">", "0"),
        ("beta4", "1/2", "<", "1/2"),
        ("beta4", "1/2", ">", "89/180"),  # (89/90)/2: err -2.8e-5 -> (3,1)
        ("beta4", "3", "<", "3"),
        ("beta6", "2", ">", "1"),
        ("beta6", "2", "<", "2"),
        ("beta6", "1/2", "<", "1/2"),
        # negative coefficient flips the inequality's constant side only:
        # -β(4) < -89/90 ⇔ β(4) > 89/90 — same (3,1) proof, mirrored
        ("beta4", "-1", "<", "0"),
        ("beta4", "-1", "<", "-89/90"),
        # degenerate power=0: true constant identity, like zeta_odd's 0-vs-r
        ("beta4", "0", "<", "1/2"),
        ("beta4", "0", ">", "-1"),
        ("beta6", "0", "<", "1"),
        ("beta8", "2", ">", "1"),
        ("beta8", "0", "<", "1"),
        ("beta10", "0", ">", "-1"),
    ],
)
def test_emitted_proofs_verify(kind, power, comp, bound):
    res = run(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("beta4", "1", "<", "89/90"),  # 89/90 < β(4) by 5.6e-5
        ("beta4", "1", "<", "805/814"),  # below by 1.1e-6
        ("beta4", "1", "<", "5725/5789"),  # below by 1.7e-9
        ("beta4", "1", ">", "179/181"),  # 179/181 > β(4) by 5.7e-6
        ("beta4", "1", ">", "984/995"),  # above by 1.7e-7
        ("beta4", "2", "<", "89/45"),  # 2·(89/90) < 2β(4)
        ("beta6", "1", "<", "759/760"),  # below by 1.0e-6
        ("beta6", "1", "<", "4048587/4053917"),  # below by 2.6e-14
        ("beta6", "1", ">", "760/761"),  # above by 7.2e-7
        ("beta6", "1", ">", "559814/560551"),  # above by 4.1e-13
        ("beta8", "1", "<", "86648/86661"),  # below by 5.5e-11
        ("beta8", "1", ">", "199957/199987"),  # above by 2.5e-12
        ("beta8", "1", ">", "1"),  # β(8) < 1
        ("beta10", "1", "<", "49999/50000"),  # below by 3.2e-6
        ("beta10", "1", ">", "59396/59397"),  # above by 1.1e-10
        ("beta10", "1", ">", "1"),  # β(10) < 1
    ],
)
def test_false_claim_raises_wrong_direction(kind, power, comp, bound):
    # prove() itself does no direction screening; the mid-scan WrongDirection
    # is the exact-moment certificate that the flipped inequality holds
    with pytest.raises(WrongDirection):
        beta_even.prove(kind, Fraction(power), comp, Fraction(bound))


def test_true_claim_beyond_budget_reports_no_solution():
    # β(4) > 1873596/1894541 is true (err 7.7e-14) but first resolves at
    # (11,11), past the (m,n) ≤ 10 budget
    with pytest.raises(NoSolution):
        beta_even.prove("beta4", Fraction(1), ">", Fraction(1873596, 1894541))


def test_corrupted_params_fail_identity():
    resp = beta_even.prove("beta4", Fraction(1), ">", Fraction(89, 90))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check("beta4", Fraction(1), ">", Fraction(89, 90), params)
    assert not res["identity_ok"]


def test_proof_does_not_verify_flipped_claim():
    resp = beta_even.prove("beta4", Fraction(1), ">", Fraction(89, 90))
    res = check("beta4", Fraction(1), "<", Fraction(89, 90), resp["parameters"])
    assert not res["identity_ok"]


def test_render_equation_smoke():
    resp = beta_even.prove("beta4", Fraction(1), ">", Fraction(89, 90))
    tex = beta_even.render_equation(resp["parameters"], "beta4", Fraction(1), ">", Fraction(89, 90))
    assert tex.startswith("\\beta(4) - \\dfrac{89}{90} = \\int_0^1 ")
    # r=3 is odd: the weight must print as ln^3(1/x), not ln^3(x)
    assert tex.endswith("\\ln^{3}(1/x) \\mathrm{d} x > 0")


@pytest.mark.parametrize("kind", KINDS)
def test_verify_dispatch(monkeypatch, kind):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, kind, "beta_even")
    resp = beta_even.prove(kind, Fraction(1), ">", Fraction(9, 10))
    res = verify(kind, Fraction(1), ">", Fraction(9, 10), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize("kind", KINDS)
def test_full_pipeline(kind):
    """Through the registered family, solver.prove(exact=True) must emit
    self-checking proofs whose certificates verify_cert accepts."""
    resp = solver.prove(kind, "1", ">", "9/10", exact=True)
    res = check(kind, Fraction(1), ">", Fraction(9, 10), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
    assert verify_cert(resp["certificate"])
