"""Exact-mode tests for the log family.

Covers ln_q, ln_q_square, artanh_q, arcoth_q. mode=exact must emit proofs
that pass exact_check — including ln_q_square at q in {5, 7}, where the
site 500s on a singular moment system (the (1,0) determinant is ∝ c-4,
(1,1) ∝ c-6) but the true system solves fine once singular candidates are
skipped. False claims are rejected by the certified direction comparison
before the search — so q in {5, 7} false claims surface as WrongDirection,
never the site's InternalError.
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response

GOLDEN = Path(__file__).resolve().parent.parent / "bench" / "data" / "golden.jsonl"
TYPES = ("ln_q", "ln_q_square", "artanh_q", "arcoth_q")


def check(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    return verify_response(kind, Fraction(power), comp, Fraction(bound), resp)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("ln_q", "3", "<", "11/10"),
        ("ln_q", "2", ">", "11/16"),
        ("ln_q", "5", "<", "375/233"),
        ("ln_q", "3/2", ">", "2/5"),
        ("ln_q_square", "2", "<", "1/2"),
        ("ln_q_square", "3", "<", "121/100"),
        ("ln_q_square", "3/2", "<", "17/100"),
        # q in {5, 7}: site 500s on the singular moment system; exact mode solves
        ("ln_q_square", "5", ">", "2"),
        ("ln_q_square", "5", "<", "3"),
        ("ln_q_square", "7", ">", "3"),
        ("ln_q_square", "7", "<", "4"),
        ("artanh_q", "1/2", "<", "11/20"),
        ("artanh_q", "1/2", ">", "1/2"),
        ("artanh_q", "9/10", "<", "3/2"),
        ("arcoth_q", "2", "<", "11/20"),
        ("arcoth_q", "2", ">", "1/2"),
        ("arcoth_q", "3/2", "<", "9/10"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    res = check(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("ln_q", "3", ">", "11/10"),  # ln 3 < 11/10
        ("ln_q", "3/2", "<", "2/5"),  # ln(3/2) > 2/5
        ("ln_q_square", "2", ">", "1/2"),  # ln^2 2 < 1/2
        # false claims at the site's crash values: WrongDirection, not 500 —
        # certified_cmp decides direction before the kernel ever runs
        ("ln_q_square", "5", ">", "3"),
        ("ln_q_square", "5", "<", "5709/2204"),
        ("ln_q_square", "7", ">", "4"),
        ("artanh_q", "1/2", ">", "11/20"),  # artanh(1/2) < 11/20
        ("arcoth_q", "2", ">", "11/20"),  # arcoth 2 < 11/20
    ],
)
def test_false_claim_rejected_before_search(kind, power, comp, bound):
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound, exact=True)


def test_tight_true_claim_beyond_budget_reports_no_solution():
    # ln^2 5 > 5709/2204 is true by ~1.3e-8; the first proof sits at (6, 15),
    # outside the m,n <= 10 budget — honest "not found", never a fake proof
    # and never the site's 500.
    with pytest.raises(NoSolution):
        solve.prove("ln_q_square", "5", ">", "5709/2204", exact=True)


def test_true_claim_budget_exhaustion_reports_no_solution():
    # ln 22 < 2105/681 is true but needs exponents > 10 (site-mode test pins
    # NoSolution already); exact mode must fail the same honest way
    try:
        resp = solve.prove("ln_q", "22", "<", "2105/681", exact=True)
    except NoSolution:
        return
    res = verify_response("ln_q", Fraction(22), "<", Fraction(2105, 681), resp)
    assert res["identity_ok"]


def test_corrupted_params_fail_identity():
    resp = solve.prove("ln_q", "3", "<", "11/10", exact=True)
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = verify("ln_q", Fraction(3), "<", Fraction(11, 10), params)
    assert not res["identity_ok"]


def test_corrupted_qtilde_fails_identity():
    # artanh params smuggle q~ through c_val — a wrong c_val builds the wrong
    # printed denominator and the identity must fail
    resp = solve.prove("artanh_q", "1/2", "<", "11/20", exact=True)
    assert resp["parameters"]["c_val"] == "3"  # (1+q)/(1-q) at q = 1/2
    params = dict(resp["parameters"], c_val="4")
    res = verify("artanh_q", Fraction(1, 2), "<", Fraction(11, 20), params)
    assert not res["identity_ok"]


def test_params_checked_against_wrong_direction_fail():
    resp = solve.prove("ln_q", "3", "<", "11/10", exact=True)
    res = verify("ln_q", Fraction(3), ">", Fraction(11, 10), resp["parameters"])
    assert not res["identity_ok"]


@pytest.mark.skipif(not GOLDEN.exists(), reason="bench/data 语料不入库（rsync 同步）")
def test_golden_corpus_sweep():
    """Every golden success record in the four types must verify exactly.

    The log family has no reproduced site distortion (the only quirk is the
    q in {5, 7} crash, which emits no success records), so the sweep should
    be unanimous; a failure here means a site-emitted false identity worth
    pinning as a regression test.
    """
    failures = []
    total = 0
    for line in GOLDEN.open():
        rec = json.loads(line)
        if rec["type"] not in TYPES or not rec.get("success"):
            continue
        total += 1
        res = verify(
            rec["type"],
            Fraction(rec["power"]),
            rec["comparison"],
            Fraction(rec["rational"]),
            rec["parameters"],
        )
        if not (res["identity_ok"] and res["nonneg"]):
            failures.append((rec["type"], rec["power"], rec["comparison"], rec["rational"]))
    assert not failures, f"{len(failures)}/{total} records failed: {failures[:10]}"
