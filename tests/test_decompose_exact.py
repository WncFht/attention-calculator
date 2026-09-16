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

from attention_calculator import decompose
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


# ---------------------------------------------------------------- exact kinds
# mode=exact admits every kernel-registered type as an atom (W7-decomp-ext);
# decompose.py stays byte-frozen, so these spellings are exact-mode only.


@pytest.mark.parametrize(
    "problem, kind, power",
    [
        ("zeta5 > 1", "zeta5", "1"),
        ("zeta7 > 1", "zeta7", "1"),
        ("zeta9 > 1", "zeta9", "1"),
        ("zeta11 < 2", "zeta11", "1"),
        ("zeta(7) > 1", "zeta7", "1"),  # builtin spelling
        ("beta4 > 9/10", "beta4", "1"),
        ("beta6 > 9/10", "beta6", "1"),
        ("beta8 > 9/10", "beta8", "1"),
        ("beta10 > 9/10", "beta10", "1"),
        ("dirichlet_beta(8) > 9/10", "beta8", "1"),  # builtin spelling
        ("pi_sqrt2 > 4", "pi_sqrt2", "1"),
        ("pi3 > 5", "pi3", "1"),
        ("pi3_u > 3/5", "pi3_u", "1"),
        ("pi3_a > 7/2", "pi3_a", "1"),
        ("gamma14 > 7/2", "gamma14", "1"),
        ("gamma34 > 1", "gamma34", "1"),
        ("gamma12 < 2", "gamma12", "1"),
        ("Gamma(1/4) > 7/2", "gamma14", "1"),  # builtin spelling
        ("ln(3)^3 > 1", "ln_q_cube", "3"),
        ("ln(2)^2 > 9/20", "ln_q_square", "2"),
        ("ln(3)^4 > 7/5", "ln_q_quad", "3"),
        ("arcsin(1/2) > 1/2", "arcsin_q", "1/2"),
        ("arsinh(1) > 4/5", "arsinh_q", "1"),
        ("gaussint(1) > 7/10", "gaussint_q", "1"),
        ("dawson(1) > 1/2", "dawson_q", "1"),
        ("erfiint(1) > 7/5", "erfiint_q", "1"),
        ("li2(1/2) > 1/2", "li2_q", "1/2"),
        ("psi1(1) > 3/2", "psi1_q", "1"),
        ("si(2) > 3/2", "si_q", "2"),
        ("cin(1) > 1/5", "cin_q", "1"),
        ("cot(1) > 3/5", "cot_q", "1"),
        ("cosh(1) > 3/2", "cosh_q", "1"),
        ("coth(1) > 13/10", "coth_q", "1"),
        ("arccot(1) > 3/4", "arccot_q", "1"),
        ("artanh(1/2) > 1/2", "artanh_q", "1/2"),
        ("arcoth(2) > 1/2", "arcoth_q", "2"),
        ("sinpi(1/4) > 3/5", "sin_pi_q", "1/4"),
        ("cospi(1/4) > 3/5", "cos_pi_q", "1/4"),
        ("sind(45) > 3/5", "sin_q_degree", "45"),
        ("cosd(45) > 3/5", "cos_q_degree", "45"),
    ],
)
def test_exact_kind_atoms(problem, kind, power):
    """Each registered exact type works as a single-atom claim (direct path)."""
    r = decompose_exact(problem)
    assert r["all_proved"], r["failures"]
    assert r["certifies"]
    check_certifies(r, r["comparison"])
    (step,) = r["steps"]
    assert step["type"] == kind
    assert step["power"] == power


@pytest.mark.parametrize(
    "problem",
    [
        "si(-2) < -3/2",  # Si is odd: kernel folds sigma=sign(q)
        "arsinh(-1) < -4/5",  # asinh odd
        "li2(-1) < 0",  # signed q inside the li2 kernel's domain
        "sinh(-1) < -1",  # parity-folding site kind, exact path
        "tanh(-1) < -3/4",
        "arctan(-1) < -7/10",
        "cosh(-1) < 8/5",  # even: cosh(-1) = cosh(1) ~ 1.543 < 1.6
    ],
)
def test_signed_argument_atoms(problem):
    """Signed-domain kinds pass a negative q to the power slot — never a
    reciprocal factor (the site path used to mis-model these as 1/U)."""
    r = decompose_exact(problem)
    assert r["all_proved"], r["failures"]
    assert r["certifies"]
    check_certifies(r, r["comparison"])


def test_function_atom_in_sum():
    """A function kind composes into a multi-term allocation."""
    r = decompose_exact("ln(3)^3 + zeta5 > 2")
    assert r["all_proved"], r["failures"]
    assert r["certifies"]
    check_certifies(r, r["comparison"])
    kinds = {s["type"] for s in r["steps"]}
    assert kinds == {"ln_q_cube", "zeta5"}


def test_exact_coef_product():
    """A product of exact coefficient kinds proves (pi*zeta5 ~ 3.2573)."""
    r = decompose_exact("pi*zeta5 > 16/5")
    assert r["all_proved"], r["failures"]
    assert r["certifies"]
    check_certifies(r, r["comparison"])


def test_coef_folds_into_power():
    """2*zeta5 folds the coefficient into the sub-claim's power slot."""
    r = decompose_exact("2*zeta5 + e > 4")
    assert r["all_proved"], r["failures"]
    zeta_step = next(s for s in r["steps"] if s["type"] == "zeta5")
    assert zeta_step["power"] == "2"
    assert zeta_step["coefficient"] == "2"


def test_exact_kind_pair_mode():
    """atom-vs-atom across new kinds chains through a rational midpoint."""
    r = decompose_exact("si(1) > cin(1)")
    assert r["success"], r["failures"]
    assert len(r["steps"]) == 2
    assert {s["type"] for s in r["steps"]} == {"si_q", "cin_q"}
    mid = Fraction(r["steps"][0]["bound"])
    assert Fraction(r["steps"][1]["bound"]) == mid


@pytest.mark.parametrize(
    "problem",
    [
        "zeta13 + e > 3",  # unregistered symbol -> site ERR_ATOM
        "li3(2) > 0",  # unregistered function kind -> ERR_ATOM
        "zeta5^2 + pi > 9",  # exact constants reject powers -> ERR_POW_SPECIAL
        "ln(2)^5 + e > 2",  # (ln q)^5 has no registered kernel -> ERR_ATOM
        "li2(1/2)*pi > 1",  # function atoms still banned inside products
        "li2(3/2) > 0",  # domain: Li2 needs q<1
        "arcsin(2) > 1",  # domain: (0,1)
        "si(0) > 1",  # domain: q!=0
        "psi1(-1) > 0",  # domain: q>0 (psi' pole)
        "ln(1/2) > 0",  # domain: q>1 (site rule, kept in exact mode)
        "coth(-1) < -1",  # domain: q>0 (kernel proves positive-q only)
        "pi3_u > 7/10",  # certified-false claim (U ~ 0.6845), not a parse error
    ],
)
def test_exact_rejections(problem):
    """Unknown kinds, domain-invalid args, and false claims reject cleanly."""
    if problem == "pi3_u > 7/10":
        r = decompose_exact(problem)
        assert not r["success"]
        assert r["failures"] == [{"reason": "claim_certified_false"}]
        return
    with pytest.raises(ValueError):
        decompose_exact(problem)


def test_site_path_still_rejects_exact_spellings():
    """Guard: the byte-frozen site decomposer must NOT accept new spellings."""
    for problem in ("zeta5 > 1", "li2(1/2) > 1/2", "Gamma(1/4) > 3", "ln(3)^3 > 1"):
        with pytest.raises(ValueError):
            decompose.decompose_inequality(problem)
