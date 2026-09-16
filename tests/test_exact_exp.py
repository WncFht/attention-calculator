"""Exact-mode smoke tests for the exp family (reference implementation).

mode=exact must emit proofs that pass exact_check; false claims must be
rejected by the certified direction comparison before any search.
"""

from fractions import Fraction

import pytest

from attention_calculator import solve
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify_response


def check(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    return verify_response(kind, Fraction(power), comp, Fraction(bound), resp)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("e", "1", ">", "8/3"),
        ("e", "1", ">", "5/2"),
        ("e", "2", "<", "6"),
        ("e_q", "2", ">", "7"),
        ("e_pi", "1", ">", "23"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    res = check(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


def test_false_claim_rejected_before_search():
    with pytest.raises(WrongDirection):
        solve.prove("e", "1", "<", "8/3", exact=True)


def test_true_claim_budget_exhaustion_reports_no_solution():
    # e < 3 is true but may exceed the budget — the failure mode must be
    # NoSolution (honest "not found"), never a fake proof
    try:
        resp = solve.prove("e", "1", "<", "3", exact=True)
    except NoSolution:
        return
    res = verify_response("e", Fraction(1), "<", Fraction(3), resp)
    assert res["identity_ok"]
