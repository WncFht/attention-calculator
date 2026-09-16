"""Exact-mode tests for the quadlog family.

Covers pi, pi_n, catalan, zeta3, arctan_q, arccot_q. mode=exact must emit
proofs that pass exact_check; false claims are rejected by the certified
direction comparison before any search. The corpus sweep runs check() on
every golden success record — the family solves the printed claim
s*(power*C - bound) directly, so all site params satisfy identity_ok,
including power=0 (printed 0*C vs bound degenerates to a true constant
identity) and power!=1 records (bench/out/verify.jsonl's false-identity
flags there were verify-side artifacts: claimed_lhs ignored power).
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.engine import EqualClaim, NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response

GOLDEN = Path(__file__).resolve().parent.parent / "bench" / "data" / "golden.jsonl"
TYPES = ("pi", "pi_n", "catalan", "zeta3", "arctan_q", "arccot_q")


def check(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    return verify_response(kind, Fraction(power), comp, Fraction(bound), resp)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("pi", "1", ">", "3"),
        ("pi", "1", "<", "22/7"),
        ("pi", "3", ">", "9"),
        ("pi", "1/2", "<", "8/5"),
        ("pi", "8", ">", "25"),
        ("pi_n", "1", ">", "3"),
        ("pi_n", "2", ">", "9"),
        ("pi_n", "2", "<", "10"),
        ("pi_n", "3", ">", "31"),
        ("pi_n", "5/2", ">", "17"),
        ("pi_n", "5/2", "<", "35/2"),
        ("pi_n", "8", "<", "9489"),
        ("catalan", "1", ">", "9/10"),
        ("catalan", "1", "<", "11/12"),
        ("catalan", "2", "<", "2"),
        ("catalan", "1/2", ">", "2/5"),
        ("catalan", "0", "<", "1/2"),
        ("zeta3", "1", "<", "5/4"),
        ("zeta3", "2", ">", "2"),
        ("zeta3", "1/2", ">", "3/5"),
        ("zeta3", "0", "<", "1"),
        ("arctan_q", "1", ">", "7/9"),
        ("arctan_q", "1", "<", "4/5"),
        ("arctan_q", "2", ">", "1"),
        ("arctan_q", "1/2", "<", "1/2"),
        ("arccot_q", "1", ">", "7/9"),
        ("arccot_q", "2", "<", "1/2"),
        ("arccot_q", "1/2", ">", "1"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    res = check(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("pi", "1", "<", "3"),  # pi > 3
        ("pi_n", "2", "<", "9"),  # pi^2 > 9
        ("pi_n", "5/2", "<", "17"),  # pi^(5/2) > 17
        ("catalan", "2", ">", "2"),  # 2C < 2
        ("zeta3", "1", ">", "5/4"),  # zeta(3) < 5/4
        ("arctan_q", "1", ">", "4/5"),  # arctan 1 < 4/5
        ("arccot_q", "1/2", "<", "1"),  # arccot(1/2) > 1
    ],
)
def test_false_claim_rejected_before_search(kind, power, comp, bound):
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound, exact=True)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("catalan", "0", "<", "0"),
        ("catalan", "0", ">", "0"),
        ("zeta3", "0", "<", "0"),
    ],
)
def test_zero_power_equality_claim(kind, power, comp, bound):
    # power=0 makes the claimed constant 0*C = 0; bound 0 is exact equality.
    with pytest.raises(EqualClaim):
        solve.prove(kind, power, comp, bound, exact=True)


def test_true_claim_budget_exhaustion_reports_no_solution():
    # arctan(3) < 5/4 is true but exceeds the (m,n) budget — the failure must
    # be NoSolution (honest "not found"), never a fake proof
    try:
        resp = solve.prove("arctan_q", "3", "<", "5/4", exact=True)
    except NoSolution:
        return
    res = verify_response("arctan_q", Fraction(3), "<", Fraction(5, 4), resp)
    assert res["identity_ok"]


def test_corrupted_params_fail_identity():
    resp = solve.prove("pi", "1", ">", "3", exact=True)
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = verify("pi", Fraction(1), ">", Fraction(3), params)
    assert not res["identity_ok"]


def test_params_checked_against_wrong_direction_fail():
    # a proof of pi - 3 = ∫f > 0 does not prove the flipped claim 3 - pi
    resp = solve.prove("pi", "1", ">", "3", exact=True)
    res = verify("pi", Fraction(1), "<", Fraction(3), resp["parameters"])
    assert not res["identity_ok"]


def test_zero_power_emitted_identity_is_true():
    # site prints e.g. '1/2 - 0*C = ∫(1+x^2)ln(1/x)/(2(1+x^2)) > 0': the claim
    # is degenerate but the printed equation is a true identity (params solve
    # {C: 0, 1: bound})
    resp = solve.prove("catalan", "0", "<", "1/2", exact=True)
    res = verify_response("catalan", Fraction(0), "<", Fraction(1, 2), resp)
    assert res["identity_ok"]
    assert res["integrand"] == {"1": Fraction(1, 2)}


@pytest.mark.skipif(not GOLDEN.exists(), reason="bench/data 语料不入库（rsync 同步）")
def test_golden_corpus_sweep():
    """Every golden success record in the six types must verify exactly —
    except the four vacuous-equality emissions (0⋚0 claims answered with
    `0 = ∫0 dx > 0`): their identity holds but the integrand is identically
    zero, so nonneg must flag them."""
    vacuous = {
        ("catalan", "0", "<", "0"),
        ("catalan", "0", ">", "0"),
        ("zeta3", "0", "<", "0"),
        ("zeta3", "0", ">", "0"),
    }
    failures = []
    total = 0
    for line in GOLDEN.open():
        rec = json.loads(line)
        if rec["type"] not in TYPES or not rec.get("success"):
            continue
        total += 1
        params = json.loads(rec["raw_calculate"])["parameters"]
        res = verify(
            rec["type"],
            Fraction(rec["power"]),
            rec["comparison"],
            Fraction(rec["rational"]),
            params,
        )
        key = (rec["type"], rec["power"], rec["comparison"], rec["rational"])
        if key in vacuous:
            assert res["identity_ok"] and not res["nonneg"], f"vacuous record changed: {key}"
            continue
        if not (res["identity_ok"] and res["nonneg"]):
            failures.append(key)
    assert not failures, f"{len(failures)}/{total} records failed: {failures[:10]}"
