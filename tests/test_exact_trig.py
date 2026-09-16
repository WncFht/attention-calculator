"""Exact-mode tests for trig kernels: trig_q (sin_q/cos_q/tan_q/cot_q)
and trig_pi (sin_pi_q/cos_pi_q/sin_q_degree/cos_q_degree).

The golden corpus doubles as the regression set: every site-emitted trig
proof outside the known trig_pi (m,n)=(1,8) bias cluster must pass the
exact checker, and every (1,8) record must fail its identity -- those
params encode the site's corrupted stored formula (BIAS18). mode=exact
re-solves those inputs against the true moments and must never re-emit
the biased coefficients.
"""

import json
from collections.abc import Iterator
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.engine import EqualClaim, NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response

GOLDEN = Path(__file__).resolve().parent.parent / "bench" / "data" / "golden.jsonl"
TRIG_Q = ("sin_q", "cos_q", "tan_q", "cot_q")
TRIG_PI = ("sin_pi_q", "cos_pi_q", "sin_q_degree", "cos_q_degree")


def golden_records() -> Iterator[dict]:
    """Successful golden /calculate records for the eight trig types."""
    with GOLDEN.open() as f:
        for line in f:
            r = json.loads(line)
            if r.get("type") in TRIG_Q + TRIG_PI and r.get("success"):
                yield r


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("sin_q", "1/2", ">", "11/23"),
        ("sin_q", "1/2", "<", "1/2"),
        ("cos_q", "1/2", ">", "4/5"),
        ("tan_q", "1/2", ">", "1/2"),
        ("tan_q", "1/2", "<", "3/5"),
        ("cot_q", "1/2", ">", "3/2"),
        ("sin_pi_q", "1/5", ">", "1/2"),
        ("cos_pi_q", "1/10", ">", "9/10"),
        ("sin_q_degree", "10", ">", "1/6"),
        ("cos_q_degree", "18", ">", "9/10"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    res = verify_response(kind, Fraction(power), comp, Fraction(bound), resp)
    assert res["identity_ok"]
    assert res["nonneg"]


def test_golden_corpus_sweep():
    """Sweep every golden trig record: only the 30 known (1,8) biased
    emissions may fail, and only on the identity (their P stays nonneg)."""
    failures = {}
    total = 0
    for r in golden_records():
        total += 1
        res = verify(
            r["type"],
            Fraction(r["power"]),
            r["comparison"],
            Fraction(r["rational"]),
            r["parameters"],
        )
        key = (r["type"], r["power"], r["comparison"], r["rational"])
        if not res["identity_ok"]:
            assert res["nonneg"], (
                f"sign-definite but wrong identity expected, got nonneg=False: {key}"
            )
            failures[key] = (r["parameters"]["m"], r["parameters"]["n"])
        assert res["nonneg"], f"indefinite P in emitted proof: {key}"
    assert total == 521
    # every failing record is the corrupted stored (1,8) formula
    assert set(failures.values()) == {(1, 8)}
    assert len(failures) == 30


BIAS18_PINS = [
    ("sin_q_degree", "18", "<", "305/987"),
    ("cos_q_degree", "72", "<", "305/987"),
    ("sin_pi_q", "1/10", "<", "305/987"),
    ("cos_pi_q", "2/5", "<", "305/987"),
]


@pytest.mark.parametrize(("kind", "power", "comp", "bound"), BIAS18_PINS)
def test_biased_site_params_fail_identity(kind, power, comp, bound):
    """The site's (1,8) emissions are false identities under true moments."""
    site = next(
        r["parameters"]
        for r in golden_records()
        if (r["type"], r["power"], r["comparison"], r["rational"]) == (kind, power, comp, bound)
    )
    assert (site["m"], site["n"]) == (1, 8)
    res = verify(kind, Fraction(power), comp, Fraction(bound), site)
    assert not res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(("kind", "power", "comp", "bound"), BIAS18_PINS)
def test_exact_mode_never_emits_biased_params(kind, power, comp, bound):
    """Re-solving a known false-identity input in exact mode must yield a
    check-passing proof or an honest failure -- never the biased params."""
    try:
        resp = solve.prove(kind, power, comp, bound, exact=True)
    except (NoSolution, WrongDirection):
        return
    res = verify_response(kind, Fraction(power), comp, Fraction(bound), resp)
    assert res["identity_ok"]
    assert res["nonneg"]


def test_false_claim_rejected_before_search():
    with pytest.raises(WrongDirection):
        solve.prove("sin_pi_q", "1/5", ">", "3/5", exact=True)
    with pytest.raises(WrongDirection):
        solve.prove("sin_q", "1/2", "<", "0", exact=True)


def test_equal_claim_reports_equal():
    # sin(pi/6) = 1/2 exactly (Niven point)
    with pytest.raises(EqualClaim):
        solve.prove("sin_pi_q", "1/6", ">", "1/2", exact=True)
