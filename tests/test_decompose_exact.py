"""Tests for decompose_exact: provable-by-construction composite decomposition.

Every successful run must satisfy the two invariants the module is built
around: each emitted sub-claim actually passed ``solve.prove(...,
exact=True)`` inside the run (``all_proved``), and the conjunction of the
rational term bounds certifies the original claim over QQ (``certifies``).
The golden problems are the two site fixtures; the harder cases tighten the
slack toward the measured (m,n)-search frontier, exercise reciprocal
factors, and pin honest failure modes.
"""

from fractions import Fraction

import pytest

from attention_calculator.decompose_exact import decompose_exact


def check_certifies(r, comp):
    """Independent re-check: the emitted rational term bounds imply the claim."""
    bounds = [Fraction(b) for b in r["term_bounds"]]
    total = sum(bounds)
    assert Fraction(r["bound_sum"]) == total
    if comp == ">":
        assert total >= Fraction(r["rhs"])
    else:
        assert total <= Fraction(r["rhs"])


@pytest.mark.parametrize(
    "problem",
    [
        "pi^2+8*pi>35",  # golden fixture 000
        "e*pi+phi+sin(1)<11",  # golden fixture 001 (product + func term)
        "pi^2+8*pi>35002/1000",  # golden 1 with ~7x tighter slack
        "e*pi<4443/520",  # two-factor product, slack ~3e-4
        "pi^2/e>36/10",  # reciprocal factor: 1/e needs an upper bound on e
        "2/pi+pi<4",  # negative-power term among the summands
        "varpi+gauss<35/10",  # hard constants: lemniscate + Gauss frontiers
        "sin(1)+cos(1)<17/12",  # two trig units
        "8*pi-pi^2>15",  # negative-coefficient term (upper bound on pi^2)
        "gamma+pi<19/5",  # gamma composite kernel inside a decomposition
    ],
)
def test_decomposition_proves(problem):
    """Each sub-claim proves and the rational bound sum certifies the claim."""
    r = decompose_exact(problem)
    assert r["all_proved"], r["failures"]
    assert r["certifies"]
    assert r["success"]
    check_certifies(r, r["comparison"])
    for step in r["steps"]:
        assert Fraction(step["bound"]) is not None
        assert step.get("error") is None


def test_single_atom_direct():
    """A lone atom proves against R itself — the loosest valid bound."""
    r = decompose_exact("e^pi>2314069/100000")
    assert r["success"], r["failures"]
    assert r["steps"][0]["bound"] == "2314069/100000"


def test_pair_mode():
    """atom-vs-atom claim chains two proofs through a rational midpoint."""
    r = decompose_exact("pi>e")
    assert r["success"]
    assert len(r["steps"]) == 2
    mid = Fraction(r["steps"][0]["bound"])
    assert Fraction(r["steps"][1]["bound"]) == mid


def test_false_claim_certified():
    """A strictly-false claim is rejected by the certified slack check."""
    r = decompose_exact("pi^2+8*pi<35")
    assert not r["success"]
    assert r["failures"] == [{"reason": "claim_certified_false"}]


def test_below_frontier_honest_failure():
    """Slack under e_pi's provability floor reports NoSolution, not a proof.

    e^pi ~ 23.1406926328; the bound leaves ~3e-9 of slack while the
    (m,n)<=30 search frontier for e_pi '>' sits near 1e-7, so no rational
    bound covering R is reachable — the run must surface that honestly.
    """
    r = decompose_exact("e^pi>2314069263/100000000")
    assert not r["success"]
    assert not r["all_proved"]
    assert any(f.get("error") == "NoSolution" for f in r["failures"])
