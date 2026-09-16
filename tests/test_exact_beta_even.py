"""Exact-mode tests for the beta_even family: beta4, beta6 (new W7 types).

Drives kernels.beta_even.prove + exact_check.beta_even.check directly — the
types are not yet wired into solve.FAMILY/EXACT_TYPES and have no site
counterpart. Bounds are CF convergents of β(4)/β(6) computed with mpmath at
80dps (workdps, not bare mp.dps assignment):

    β(4) ≈ 0.98894455174110533611   β(6) ≈ 0.99868522221843813544

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
from attention_calculator.kernels import EXACT_TYPES, beta_even

KINDS = ("beta4", "beta6")


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
def test_full_pipeline_when_registered(kind):
    """Once kernels.EXACT_TYPES / solve.FAMILY / integrand.constant_mpf wire
    the family, solver.prove(exact=True) must emit self-checking proofs whose
    certificates verify_cert accepts."""
    if kind not in solver.FAMILY or kind not in EXACT_TYPES:
        pytest.skip(f"{kind} not yet registered (merge wiring is the leader's)")
    resp = solver.prove(kind, "1", ">", "9/10", exact=True)
    res = check(kind, Fraction(1), ">", Fraction(9, 10), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]
    assert verify_cert(resp["certificate"])
