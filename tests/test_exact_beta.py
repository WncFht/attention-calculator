"""Exact-mode tests for beta kernels (golden/varpi/gauss) and gamma.

Pins two site-bug clusters that exact_check must reject:
- gauss '<' m=5 records come from the site's transposed solve — their
  params satisfy a transposed system, not the real moment equations.
- power != 1 records printed by the site were flagged by bench/verify.py,
  but the printed integrand literally carries the power coefficient (gamma
  folds it into every term; beta kernels solve the scaled target), so the
  identities are true — only the transposed cluster is genuinely false.

And the exact-mode behavior: no transposed params, honest NoSolution when
the '<' search budget is exhausted, and self-checked proofs otherwise.
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.engine import EqualClaim, NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response

GOLDEN = Path(__file__).parent.parent / "bench" / "data" / "golden.jsonl"

# the site's transposed '<' shots (all gauss, m=5 n=1) — genuinely false
TRANSPOSED = {
    ("1", "1398/1675"),
    ("1", "21/25"),
    ("1", "167/200"),
    ("2", "2796/1675"),
    ("2", "17/10"),
    ("2", "167/100"),
    ("3", "21293/8504"),
    ("3", "251/100"),
    ("3", "313/125"),
    ("1/2", "699/1675"),
    ("1/2", "21/50"),
    ("1/2", "209/500"),
}

# zero-integrand degenerate prints: ``0 = int 0 dx > 0`` asserts 0 > 0
ZERO_INTEGRAND = {
    ("golden", "0", "<", "0"),
    ("golden", "0", ">", "0"),
    ("varpi", "0", "<", "0"),
    ("varpi", "0", ">", "0"),
    ("gauss", "0", "<", "0"),
    ("gauss", "0", ">", "0"),
}


def corpus():
    for line in GOLDEN.open():
        r = json.loads(line)
        if r.get("type") in ("golden", "varpi", "gauss", "gamma") and r.get("success"):
            yield r


def test_corpus_sweep():
    """Every golden success record for the 4 types is checked; only the
    transposed cluster may fail identity_ok, only the zero-integrand
    degenerates may fail nonneg."""
    false_identity, bad_nonneg, total = [], [], 0
    for r in corpus():
        total += 1
        res = verify(
            r["type"],
            Fraction(r["power"]),
            r["comparison"],
            Fraction(r["rational"]),
            r["parameters"],
        )
        if not res["identity_ok"]:
            false_identity.append((r["type"], r["power"], r["comparison"], r["rational"]))
        if not res["nonneg"]:
            bad_nonneg.append((r["type"], r["power"], r["comparison"], r["rational"]))
    assert total == 123
    assert set(false_identity) == {("gauss", p, "<", b) for p, b in TRANSPOSED}
    assert set(bad_nonneg) == ZERO_INTEGRAND


@pytest.mark.parametrize(("power", "rational"), sorted(TRANSPOSED))
def test_transposed_records_fail(power, rational):
    params = next(
        r["parameters"]
        for r in corpus()
        if r["type"] == "gauss"
        and r["power"] == power
        and r["comparison"] == "<"
        and r["rational"] == rational
    )
    res = verify("gauss", Fraction(power), "<", Fraction(rational), params)
    assert not res["identity_ok"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("golden", "1", ">", "3/2"),
        ("golden", "2", ">", "3"),
        ("golden", "1", "<", "2"),
        ("varpi", "1", ">", "5/2"),
        ("varpi", "2", "<", "6"),
        ("gauss", "1", ">", "4/5"),
        ("gauss", "1", "<", "9/10"),
        ("gamma", "1", ">", "0"),
        ("gamma", "1", "<", "1"),
        ("gamma", "1", ">", "57/100"),
        ("gamma", "2", ">", "1"),  # power != 1 stays consistent
        ("gamma", "2", "<", "7/6"),
        ("gamma", "1/2", "<", "1/3"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    res = verify_response(kind, Fraction(power), comp, Fraction(bound), resp)
    assert res["identity_ok"]
    assert res["nonneg"]


def test_gauss_window_never_emits_transposed():
    """The (G, ~0.855] window: exact mode either finds a true proof at
    deeper m or reports NoSolution — never the transposed params."""
    # the two tightest claims need m ~ 10^6 (the d_m/q_m feasibility ratio
    # converges to G^-1 like 1 - c/m) — honest NoSolution within budget
    for bound in ("1398/1675", "167/200"):
        with pytest.raises(NoSolution):
            solve.prove("gauss", "1", "<", bound, exact=True)
    # the rest of the window resolves at deeper m with real proofs
    for bound, m in (("21/25", 18),):
        resp = solve.prove("gauss", "1", "<", bound, exact=True)
        assert int(resp["parameters"]["m"]) == m
        res = verify_response("gauss", Fraction(1), "<", Fraction(bound), resp)
        assert res["identity_ok"] and res["nonneg"]
    resp = solve.prove("gauss", "2", "<", "17/10", exact=True)
    assert int(resp["parameters"]["m"]) == 6
    res = verify_response("gauss", Fraction(2), "<", Fraction(17, 10), resp)
    assert res["identity_ok"] and res["nonneg"]


def test_gamma_nonpositive_power_honest_failure():
    # power <= 0 '<' claims are true but unprovable in the printed form —
    # the coefficient scales every printed term, flipping the sign
    with pytest.raises(NoSolution):
        solve.prove("gamma", "-2", "<", "1", exact=True)
    with pytest.raises(NoSolution):
        solve.prove("gamma", "0", "<", "1", exact=True)
    with pytest.raises(WrongDirection):
        solve.prove("gamma", "-2", ">", "1", exact=True)


def test_equal_claim():
    with pytest.raises(EqualClaim):
        solve.prove("gauss", "0", "<", "0", exact=True)


def test_zero_integrand_site_params_fail_nonneg():
    """The site emits ``0 = int 0 dx > 0`` for 0-vs-0; the equality holds
    vacuously but the strict inequality is false."""
    params = next(
        r["parameters"]
        for r in corpus()
        if r["type"] == "gauss" and r["power"] == "0" and r["rational"] == "0"
    )
    res = verify("gauss", Fraction(0), ">", Fraction(0), params)
    assert res["identity_ok"]  # 0 = int 0 dx is a true equation
    assert not res["nonneg"]  # but it proves nothing about >
