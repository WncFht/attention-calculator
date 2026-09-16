"""Exact-mode tests for pi_sqrt2 (lemniscate-family completion, W3).

Claims ``q·S ⋚ r`` with S = π√2 = Γ(1/4)Γ(3/4) ≈ 4.4429: '>' proves via
K = (1-x⁴)^{-3/4} on x^{4m+2}(1-x)(a+bx⁴), '<' via K = (1-x⁴)^{-1/4} on
x^{4m+3}(1-x)(a+bx⁴). Every emitted proof is re-verified by
exact_check.pi_sqrt2.check — a dict-equal moment identity plus the
a >= 0 and a+b >= 0 nonneg rule on P.

The type predates its registry entries (solve.FAMILY, kernels.EXACT_TYPES
and integrand.constant_mpf are the leader's wiring); tests install the
three through monkeypatch so the full solve.prove(exact=True) ->
certificate path runs for real.
"""

from fractions import Fraction

import mpmath as mp
import pytest

from attention_calculator import certificate, integrand, solve
from attention_calculator.engine import EqualClaim, NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response
from attention_calculator.exact_check.pi_sqrt2 import check
from attention_calculator.kernels import EXACT_TYPES
from attention_calculator.kernels.pi_sqrt2 import prove, render_equation

mp.mp.dps = 60
S = mp.pi * mp.sqrt(2)


@pytest.fixture(autouse=True)
def registered(monkeypatch):
    """Wire pi_sqrt2 into the exact pipeline the way the leader will.

    constant_mpf is looked up lazily inside solve.certified_cmp, so
    patching the integrand attribute covers it; FAMILY is the shared dict
    exact_check/certificate read; EXACT_TYPES is solve's imported name.
    """
    orig = integrand.constant_mpf

    def const(kind: str, power: Fraction):
        if kind == "pi_sqrt2":
            q = mp.mpf(power.numerator) / power.denominator
            return q * mp.pi * mp.sqrt(2)
        return orig(kind, power)

    monkeypatch.setattr(integrand, "constant_mpf", const)
    monkeypatch.setitem(solve.FAMILY, "pi_sqrt2", "pi_sqrt2")
    monkeypatch.setattr(solve, "EXACT_TYPES", [*EXACT_TYPES, "pi_sqrt2"])


def tight_bound(comp: str, digs: int) -> Fraction:
    """``digs``-place floor of S for '>', ceil for '<'."""
    v = int(mp.floor(S * 10**digs))
    out = Fraction(v, 10**digs)
    return out + Fraction(1, 10**digs) if comp == "<" else out


def assert_verified(power: str, comp: str, bound: Fraction) -> dict:
    """prove() must emit and check() must confirm the exact identity."""
    resp = prove("pi_sqrt2", Fraction(power), comp, bound)
    res = check("pi_sqrt2", Fraction(power), comp, bound, resp["parameters"])
    assert res["identity_ok"], (power, comp, bound, res)
    assert res["nonneg"], (power, comp, bound, res)
    return resp


# ------------------------------------------------------- verified bounds


@pytest.mark.parametrize(
    ("comp", "bound", "m"),
    [
        (">", "4", 0),  # a=4, b=0
        (">", "22/5", 5),  # the research doc's emitted identity
        ("<", "5", 0),
        ("<", "9/2", 13),
    ],
)
def test_research_doc_bounds(comp, bound, m):
    """The four bounds verified at dps=100 in the W3 research note land at
    the documented m and pass the pipeline self-check."""
    resp = solve.prove("pi_sqrt2", "1", comp, bound, exact=True)
    p = resp["parameters"]
    assert int(p["m"]) == m
    res = verify_response("pi_sqrt2", Fraction(1), comp, Fraction(bound), resp)
    assert res["identity_ok"]
    assert res["nonneg"]


def test_gt_4_is_x2_kernel():
    """pi*sqrt(2) > 4: m=0 emits ∫4x²(1-x)(1-x⁴)^{-3/4}dx = π√2 - 4."""
    resp = solve.prove("pi_sqrt2", "1", ">", "4", exact=True)
    p = resp["parameters"]
    assert (p["a_val"], p["b_val"], p["u_val"]) == ("4", "0", "1")
    res = verify_response("pi_sqrt2", Fraction(1), ">", Fraction(4), resp)
    assert res["integrand"] == {"pi_sqrt2": Fraction(1), "1": Fraction(-4)}


def test_gt_225_emitted_identity():
    """pi*sqrt(2) > 22/5 emits the note's verified identity: m=5 with
    P = 2577291/5992448 + 5499385/749056·x⁴ on x^{22}(1-x)(1-x⁴)^{-3/4}."""
    resp = solve.prove("pi_sqrt2", "1", ">", "22/5", exact=True)
    p = resp["parameters"]
    assert p["a_val"] == "2577291/5992448"
    assert p["b_val"] == "5499385/749056"
    res = verify_response("pi_sqrt2", Fraction(1), ">", Fraction(22, 5), resp)
    assert res["integrand"] == {"pi_sqrt2": Fraction(1), "1": Fraction(-22, 5)}


def test_lt_5_emitted_identity():
    """pi*sqrt(2) < 5 at m=0: P = (13 + 56x⁴)/3 on x³(1-x)(1-x⁴)^{-1/4}."""
    resp = solve.prove("pi_sqrt2", "1", "<", "5", exact=True)
    p = resp["parameters"]
    assert (p["a_val"], p["b_val"], p["u_val"]) == ("13/3", "56/3", "3")
    res = verify_response("pi_sqrt2", Fraction(1), "<", Fraction(5), resp)
    assert res["integrand"] == {"pi_sqrt2": Fraction(-1), "1": Fraction(5)}


def test_n_is_direction_flag():
    """n encodes the kernel direction like beta_family: 0 '>', 1 '<'."""
    assert solve.prove("pi_sqrt2", "1", ">", "4", exact=True)["parameters"]["n"] == 0
    assert solve.prove("pi_sqrt2", "1", "<", "5", exact=True)["parameters"]["n"] == 1


def test_certificate_verifies():
    """prove_exact's embedded certificate independently re-verifies."""
    resp = solve.prove("pi_sqrt2", "1", ">", "22/5", exact=True)
    cert = resp["certificate"]
    assert cert["check"]["identity_ok"] is True
    assert cert["check"]["nonneg"] is True
    assert certificate.verify_cert(cert)


# ---------------------------------------------------------- coverage


@pytest.mark.parametrize(
    ("comp", "digs"),
    [(">", 2), ("<", 2)],
)
def test_tight_bounds(comp, digs):
    """1e-2 floor/ceil bounds prove (m = 11 '>' / 29 '<' observed)."""
    bound = tight_bound(comp, digs)
    assert_verified("1", comp, bound)


@pytest.mark.parametrize(
    ("power", "comp", "bound"),
    [
        ("2", ">", "8"),  # 2S ≈ 8.886 > 8 at m=0 (a=8, b=0)
        ("2", "<", "9"),  # 2S < 9 at m=13
        ("3", ">", "13"),  # 3S ≈ 13.33 > 13
        ("1/2", ">", "2"),  # S/2 ≈ 2.221 > 2
    ],
)
def test_power_coefficient(power, comp, bound):
    """The power slot scales the claim q·S ⋚ r directly."""
    assert_verified(power, comp, Fraction(bound))


def test_beyond_budget_is_honest_no_solution():
    """'<' bound 4443/1000 sits ~1.2e-4 off S — past the m<=256 budget."""
    with pytest.raises(NoSolution):
        prove("pi_sqrt2", Fraction(1), "<", Fraction(4443, 1000))


# ------------------------------------------------------- false claims


@pytest.mark.parametrize(
    ("power", "comp", "bound"),
    [
        ("1", ">", "9/2"),  # S ≈ 4.4429 < 4.5
        ("1", "<", "4"),
        ("1", "<", "22/5"),
        ("2", ">", "9"),  # 2S ≈ 8.886 < 9
    ],
)
def test_false_claim_wrong_direction(power, comp, bound):
    """False claims die in certified_cmp before the search runs."""
    with pytest.raises(WrongDirection):
        solve.prove("pi_sqrt2", power, comp, bound, exact=True)


def test_kernel_level_false_claim_is_no_solution():
    """Called directly, a false claim exhausts the scan: every solved P is
    sign-indefinite (a non-positive P would certify the true opposite
    inequality), so the kernel reports NoSolution, not WrongDirection —
    direction verdicts are upstream certified_cmp's job."""
    with pytest.raises(NoSolution):
        prove("pi_sqrt2", Fraction(1), ">", Fraction(9, 2))


# ----------------------------------------------------------- edge cases


@pytest.mark.parametrize(
    ("power", "comp", "bound", "exc"),
    [
        ("0", ">", "1", WrongDirection),  # 0 > 1 false
        ("0", ">", "-1", NoSolution),  # 0 > -1 true but unprovable
        ("0", "<", "1", NoSolution),  # 0 < 1 true but unprovable
        ("0", "<", "-1", WrongDirection),  # 0 < -1 false
        ("-1", ">", "-5", NoSolution),  # -S > -5 true, P sign-indefinite at all m
        ("-1", ">", "-4", WrongDirection),  # -S > -4 false
        ("-1", "<", "4", NoSolution),  # -S < 4 true, unprovable
        ("-1", "<", "-5", WrongDirection),  # -S < -5 false
    ],
)
def test_nonpositive_power(power, comp, bound, exc):
    """q <= 0 mirrors the varpi/gamma convention: certified_cmp rejects
    false claims; true ones are honest NoSolution since the moment machine
    cannot express a negative-coefficient claim with non-negative P."""
    with pytest.raises(exc):
        solve.prove("pi_sqrt2", power, comp, bound, exact=True)


def test_equal_claim():
    """0 vs 0 is the only reachable equality (S is irrational)."""
    with pytest.raises(EqualClaim):
        solve.prove("pi_sqrt2", "0", "<", "0", exact=True)


# ------------------------------------------------------------- checker


def test_corrupted_params_fail_identity():
    resp = solve.prove("pi_sqrt2", "1", ">", "22/5", exact=True)
    p = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = verify("pi_sqrt2", Fraction(1), ">", Fraction(22, 5), p)
    assert not res["identity_ok"]


def test_params_against_wrong_comp_fail():
    """'>' params checked as '<' moments (different kernel) fail."""
    resp = prove("pi_sqrt2", Fraction(1), ">", Fraction(4))
    res = check("pi_sqrt2", Fraction(1), "<", Fraction(4), resp["parameters"])
    assert not res["identity_ok"]


def test_zero_integrand_guard():
    """All-zero params: poly_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {"m": 0, "n": 0, "au_val": "0", "bu_val": "0", "u_val": "1"}
    for comp in (">", "<"):
        res = check("pi_sqrt2", Fraction(1), comp, Fraction(0), params)
        assert res["integrand"] == {}
        assert not res["nonneg"]
        assert not res["identity_ok"]


def test_target_drops_zero_keys():
    """q = 0 filters the S slot out of the target vector."""
    res = check(
        "pi_sqrt2",
        Fraction(0),
        "<",
        Fraction(1),
        {"m": 0, "n": 1, "au_val": "1", "bu_val": "0", "u_val": "1"},
    )
    assert res["target"] == {"1": Fraction(1)}
    assert "pi_sqrt2" not in res["target"]


def test_emitted_params_shape():
    resp = prove("pi_sqrt2", Fraction(1), "<", Fraction(5))
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
    assert p["c_val"] == "0" and p["cu_val"] == "0"
    assert resp["type"] == "pi_sqrt2"
    assert resp["solution"].startswith("a = ")
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u
    assert Fraction(p["b_val"]) == Fraction(p["bu_val"]) / u


# ------------------------------------------------------------ rendering


def test_render_equation_gt():
    resp = prove("pi_sqrt2", Fraction(1), ">", Fraction(22, 5))
    tex = render_equation(resp["parameters"], "pi_sqrt2", Fraction(1), ">", Fraction(22, 5))
    assert tex == (
        "\\pi\\sqrt{2} - \\dfrac{22}{5} = \\int_0^1 "
        "\\frac{x^{22} \\left(1 - x\\right) \\left(43995080 x^{4} + 2577291\\right)}"
        "{5992448 \\left(1 - x^{4}\\right)^{\\frac{3}{4}}} \\mathrm{d} x > 0"
    )


def test_render_equation_lt():
    resp = prove("pi_sqrt2", Fraction(1), "<", Fraction(5))
    tex = render_equation(resp["parameters"], "pi_sqrt2", Fraction(1), "<", Fraction(5))
    assert tex == (
        "5 - \\pi\\sqrt{2} = \\int_0^1 "
        "\\frac{x^{3} \\left(1 - x\\right) \\left(56 x^{4} + 13\\right)}"
        "{3 \\sqrt[4]{1 - x^{4}}} \\mathrm{d} x > 0"
    )


def test_render_equation_coefficient():
    """q != 1 prints the coefficient before \\pi\\sqrt{2}."""
    resp = prove("pi_sqrt2", Fraction(2), ">", Fraction(8))
    tex = render_equation(resp["parameters"], "pi_sqrt2", Fraction(2), ">", Fraction(8))
    assert tex.startswith("2\\pi\\sqrt{2} - 8 = ")
