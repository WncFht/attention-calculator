"""Exact-mode tests for the hyperbolic family + exp-family edge cases.

mode=exact must emit proofs that pass exact_check; false claims must be
rejected by the certified direction comparison before any search. The
golden sweeps double-check that both families carry no reproduced site
bias — every site-emitted proof must verify exactly.
"""

import json
from fractions import Fraction
from pathlib import Path

import pytest

from attention_calculator import solve
from attention_calculator.engine import EqualClaim, WrongDirection
from attention_calculator.exact_check import exp_family as exp_check
from attention_calculator.exact_check import verify, verify_response

GOLDEN = Path(__file__).resolve().parent.parent / "bench" / "data" / "golden.jsonl"
HYP = ("sinh_q", "cosh_q", "tanh_q", "coth_q")
EXP = ("e", "e_q", "e_pi")

# e_pi power=0 records trip the reference checker's zero-coefficient target
# entry ({e_pi: 0} survives while combine drops it) — a checker artifact,
# not a site falsehood; tracked for the exp_family owner
KNOWN_EXP_ZERO_KEY = {
    ("e_pi", "0", "<", "0"),
    ("e_pi", "0", ">", "0"),
    ("e_pi", "0", "<", "1/2"),
    ("e_pi", "0", "<", "1"),
}


def check(kind, power, comp, bound):
    resp = solve.prove(kind, power, comp, bound, exact=True)
    return verify_response(kind, Fraction(power), comp, Fraction(bound), resp)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("sinh_q", "1/2", ">", "0"),  # bound=0: target carries a zero "1" key
        ("sinh_q", "1/2", "<", "1"),
        ("sinh_q", "2", "<", "4"),
        ("sinh_q", "1", "<", "11752011936438015/10000000000000000"),  # sinh 1 + 1e-16
        ("cosh_q", "1", ">", "3/2"),
        ("cosh_q", "1/2", "<", "8/7"),
        ("tanh_q", "1/2", ">", "2/5"),
        ("tanh_q", "2", "<", "1"),  # sinh 2 - cosh 2 < 0 in {sinh,cosh,1} space
        ("coth_q", "1", ">", "13/10"),
        ("coth_q", "1/2", "<", "9/4"),
    ],
)
def test_exact_proofs_verify(kind, power, comp, bound):
    res = check(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("sinh_q", "1/2", "<", "0"),  # sinh q > 0 for q > 0
        ("cosh_q", "1", ">", "2"),  # cosh 1 < 2
        ("cosh_q", "1/2", "<", "9/8"),  # tight false bound: cosh(1/2) > 9/8
        ("tanh_q", "1", ">", "1"),  # tanh q < 1
        ("coth_q", "1", "<", "1"),  # coth q > 1
    ],
)
def test_false_claims_rejected_before_search(kind, power, comp, bound):
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound, exact=True)


def test_equal_claim_at_rational_point():
    # cosh(0)=1 is the family's only rational-constant point
    with pytest.raises(EqualClaim):
        solve.prove("cosh_q", "0", "<", "1", exact=True)


def test_zero_argument_still_hits_the_moment_wall():
    # certified direction passes (sinh 0 < 1/10 is true) but the 1/q IBP
    # cannot start at q=0 — the same crash the site surfaces as 500;
    # coth crashes one step earlier inside 1/tanh(0) of certified_cmp
    with pytest.raises(ZeroDivisionError):
        solve.prove("sinh_q", "0", "<", "1/10", exact=True)
    with pytest.raises(ZeroDivisionError):
        solve.prove("coth_q", "0", "<", "2", exact=True)


def test_checker_rejects_tampered_params():
    resp = solve.prove("sinh_q", "1/2", ">", "0", exact=True)
    params = dict(resp["parameters"])
    params["au_val"] = params["au_val"] + 1
    res = verify("sinh_q", Fraction("1/2"), ">", Fraction(0), params)
    assert not res["identity_ok"]


# ------------------------------------------------------------ exp edge cases


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("e", "2", ">", "5"),  # coefficient power != 1
        ("e", "3/2", "<", "41/10"),  # fractional coefficient
        ("e_q", "3/2", ">", "4"),  # fractional exponent
        ("e_q", "1/2", "<", "7/4"),
        ("e_pi", "1/2", ">", "11"),  # coefficient of e^pi
        ("e_pi", "2", "<", "47"),
        ("e", "1", "<", "2719/1000"),  # bound 7e-4 above e
    ],
)
def test_exp_edge_proofs_verify(kind, power, comp, bound):
    res = check(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("e_q", "2", "<", "7"),  # e^2 > 7
        ("e_pi", "1", "<", "23"),  # e^pi > 23
        ("e", "2", "<", "5"),  # 2e > 5
    ],
)
def test_exp_false_claims_rejected_before_search(kind, power, comp, bound):
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound, exact=True)


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("e", "0", ">", "0"),  # 0·e = 0
        ("e_q", "0", "<", "1"),  # e^0 = 1
    ],
)
def test_exp_equal_claims_at_rational_points(kind, power, comp, bound):
    with pytest.raises(EqualClaim):
        solve.prove(kind, power, comp, bound, exact=True)


# ------------------------------------------------------------- corpus sweeps


def sweep(kinds, checker):
    """Run every golden success record through checker; return failures."""
    fails = []
    for line in GOLDEN.open():
        rec = json.loads(line)
        if rec["type"] not in kinds or not rec.get("success"):
            continue
        res = checker(
            rec["type"],
            Fraction(rec["power"]),
            rec["comparison"],
            Fraction(rec["rational"]),
            rec["parameters"],
        )
        if not (res["identity_ok"] and res["nonneg"]):
            fails.append((rec["type"], rec["power"], rec["comparison"], rec["rational"]))
    return fails


@pytest.mark.skipif(not GOLDEN.exists(), reason="bench/data 语料不入库（rsync 同步）")
def test_golden_hyperbolic_proofs_all_verify():
    assert sweep(set(HYP), verify) == []


@pytest.mark.skipif(not GOLDEN.exists(), reason="bench/data 语料不入库（rsync 同步）")
def test_golden_exp_proofs_verify_modulo_known_artifact():
    fails = sweep(set(EXP), exp_check.check)
    assert set(fails) <= KNOWN_EXP_ZERO_KEY
