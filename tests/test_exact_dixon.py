"""Exact-mode tests for the Dixon family: pi3, pi3_u, pi3_a (W3, exact-only).

The Γ(1/3)-lattice constants — π₃ = B(1/3,1/3), U = A/π₃, A = 2π/√3 — share
one moment machine on the two cubic kernels (1−x³)^{−1/3} / (1−x³)^{−2/3},
with √3/(2π) as the varpi-1/π analogue flipping the reciprocal-direction
spans (docs/2026-09-16-w3-research-gamma-third.md). Every emitted proof is
re-verified by exact_check.dixon.check: dict-equal moments plus the P sign
certificate.

Verified-hit coverage: the research note's six end-to-end proofs land at the
documented m (5, 3, 7, 2, 19, 14); the edge ratio closes ~1/m, so ~1e-2-tight
bounds resolve inside m ~ 140 while the 8e-5 gap (pi3 < 53/10) honestly
exhausts the m <= 512 budget.
"""

from fractions import Fraction

import pytest

from attention_calculator import solve as solver
from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check import verify, verify_response
from attention_calculator.exact_check.dixon import check
from attention_calculator.kernels.dixon import prove, render_equation

KINDS = ("pi3", "pi3_u", "pi3_a")


def assert_verified(kind: str, q: str, comp: str, bound: Fraction) -> dict:
    """prove() must emit and check() must confirm the exact identity."""
    resp = prove(kind, Fraction(q), comp, bound)
    res = check(kind, Fraction(q), comp, bound, resp["parameters"])
    assert res["identity_ok"], (kind, q, comp, bound, res)
    assert res["nonneg"], (kind, q, comp, bound, res)
    return resp


# the research note's six verified end-to-end hits, with their m values
@pytest.mark.parametrize(
    ("kind", "comp", "bound", "m"),
    [
        ("pi3", "<", "27/5", 5),
        ("pi3", ">", "5", 3),  # proved as 1 - 5·pi3^{-1} > 0
        ("pi3_u", ">", "2/3", 7),
        ("pi3_u", "<", "5/7", 2),  # proved as (5/7)·U^{-1} - 1 > 0
        ("pi3_a", "<", "11/3", 19),
        ("pi3_a", ">", "18/5", 14),
    ],
)
def test_doc_hits(kind, comp, bound, m):
    resp = assert_verified(kind, "1", comp, Fraction(bound))
    assert int(resp["parameters"]["m"]) == m
    assert int(resp["parameters"]["cu_val"]) == 0


# ~1e-2-tight true bounds per direction, with the resolving m they hit
@pytest.mark.parametrize(
    ("kind", "comp", "bound", "m"),
    [
        ("pi3", ">", "132/25", 58),
        ("pi3", "<", "531/100", 57),
        ("pi3_u", ">", "17/25", 33),
        ("pi3_u", "<", "137/200", 141),
        ("pi3_a", ">", "181/50", 52),
        ("pi3_a", "<", "91/25", 64),
    ],
)
def test_tight_bounds_both_directions(kind, comp, bound, m):
    """~1e-2-tight bounds prove in both directions inside the m budget."""
    resp = assert_verified(kind, "1", comp, Fraction(bound))
    assert int(resp["parameters"]["m"]) == m


@pytest.mark.parametrize(
    ("kind", "q", "comp", "bound"),
    [
        ("pi3", "2", "<", "11"),  # 2·pi3 = 10.5998 < 11
        ("pi3", "2", ">", "21/2"),
        ("pi3", "1/2", "<", "3"),
        ("pi3_u", "3", ">", "2"),  # 3U = 2.0534 > 2
        ("pi3_u", "1/3", "<", "1/4"),
        ("pi3_a", "1/2", "<", "2"),
        ("pi3_a", "3", ">", "10"),
    ],
)
def test_nonunit_coefficients(kind, q, comp, bound):
    assert_verified(kind, q, comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "comp", "bound"),
    [
        ("pi3", "<", "5"),  # pi3 > 5
        ("pi3", ">", "27/5"),
        ("pi3_u", ">", "7/10"),  # U < 0.7
        ("pi3_u", "<", "13/20"),
        ("pi3_a", "<", "18/5"),
        ("pi3_a", ">", "11/3"),
    ],
)
def test_false_claim_rejected(kind, comp, bound):
    """False claims exhaust the search: every solved P is sign-indefinite
    (a non-positive one would certify the opposite inequality), so the kernel
    reports honest NoSolution. Through the mode=exact pipeline certified_cmp
    rejects them earlier with WrongDirection — covered by the gated test."""
    with pytest.raises((NoSolution, WrongDirection)):
        prove(kind, Fraction(1), comp, Fraction(bound))


@pytest.mark.parametrize(
    ("kind", "comp", "bound"),
    [
        ("pi3", "<", "100"),
        ("pi3", ">", "1/2"),
        ("pi3", ">", "0"),  # bound = 0: the t = 0 emission, claim is its rational part
        ("pi3_u", ">", "1/10"),
        ("pi3_u", "<", "2"),
        ("pi3_a", "<", "50"),
        ("pi3_a", ">", "1/10"),
    ],
)
def test_loose_bound_fallback(kind, comp, bound):
    """Bounds beyond the 2x2 interval edges emit the +b template: a single
    t·x^{3m+res}(1-x)·K·outer moment plus a positive rational leftover."""
    resp = assert_verified(kind, "1", comp, Fraction(bound))
    p = resp["parameters"]
    assert p["cu_val"] == "1"
    assert p["bu_val"] == "0" and p["a_val"] == "0" and p["c_val"] == "0"
    assert Fraction(p["b_val"]) > 0
    assert resp["solution"] == "a = 0, b = 0"


def test_beyond_budget_is_honest_no_solution():
    """pi3 < 53/10 is true by ~8e-5 but needs m ~ 7000 at the ~1/m rate."""
    with pytest.raises(NoSolution):
        prove("pi3", Fraction(1), "<", Fraction(53, 10))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("comp", [">", "<"])
def test_zero_coefficient_rejected(kind, comp):
    """q = 0 degenerates the claim to a rational comparison."""
    with pytest.raises(ValueError):
        prove(kind, Fraction(0), comp, Fraction(1))


@pytest.mark.parametrize("kind", KINDS)
def test_zero_integrand_guard(kind):
    """All-zero params: poly_nonneg vacuously accepts P = 0, but the bool
    guard on the integrand moment rejects the vacuous 0 = int 0 dx > 0."""
    params = {"m": 0, "n": 0, "au_val": "0", "bu_val": "0", "cu_val": "0", "u_val": "1"}
    res = check(kind, Fraction(1), ">", Fraction(0), params)
    assert res["integrand"] == {}
    assert not res["nonneg"]
    assert not res["identity_ok"]


@pytest.mark.parametrize("kind", KINDS)
def test_corrupted_params_fail_identity(kind):
    resp = prove(kind, Fraction(1), ">", Fraction(1, 10))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check(kind, Fraction(1), ">", Fraction(1, 10), params)
    assert not res["identity_ok"]


def test_params_against_wrong_comp_fail():
    resp = prove("pi3_u", Fraction(1), ">", Fraction(2, 3))
    res = check("pi3_u", Fraction(1), "<", Fraction(2, 3), resp["parameters"])
    assert not res["identity_ok"]


def test_fallback_corrupted_b_fails():
    resp = prove("pi3", Fraction(1), "<", Fraction(100))
    params = dict(resp["parameters"], b_val="99")
    res = check("pi3", Fraction(1), "<", Fraction(100), params)
    assert not res["identity_ok"]


@pytest.mark.parametrize("kind", KINDS)
def test_emitted_params_shape(kind):
    bound = {"pi3": "5", "pi3_u": "2/3", "pi3_a": "18/5"}[kind]
    resp = prove(kind, Fraction(1), ">", Fraction(bound))
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
    assert resp["type"] == kind
    assert resp["solution"].startswith("a = ")
    u = Fraction(p["u_val"])
    assert Fraction(p["a_val"]) == Fraction(p["au_val"]) / u
    assert Fraction(p["b_val"]) == Fraction(p["bu_val"]) / u
    # n is the direction flag, matching beta_family/pi_sqrt2 encoding
    assert int(p["n"]) == 0


def test_render_equation_doc_identities():
    """Pinned LaTeX for the note's verified identities (note normalizes the
    pi3 '>' identity by the bound — ours keeps the varpi '>' coefficient
    convention, 5x the note's printed scale)."""
    resp = prove("pi3", Fraction(1), "<", Fraction(27, 5))
    tex = render_equation(resp["parameters"], "pi3", Fraction(1), "<", Fraction(27, 5))
    assert tex == (
        "\\dfrac{27}{5} - \\pi_{3} = \\int_0^1 \\frac{x^{17} \\left(1 - x\\right) "
        "\\left(2185 x^{3} + 16858\\right)}{1755 \\left(1 - x^{3}\\right)^{\\frac{2}{3}}} "
        "\\mathrm{d} x > 0"
    )
    resp = prove("pi3_u", Fraction(1), ">", Fraction(2, 3))
    tex = render_equation(resp["parameters"], "pi3_u", Fraction(1), ">", Fraction(2, 3))
    assert tex == (
        "U - \\dfrac{2}{3} = \\int_0^1 \\frac{x^{22} \\left(1 - x\\right) "
        "\\left(4978610 x^{3} + 1253720\\right)}{1003833 \\sqrt[3]{1 - x^{3}}} "
        "\\mathrm{d} x > 0"
    )


def test_render_equation_reciprocal_lhs():
    """The reciprocal-direction claims print the reciprocal LHS."""
    resp = prove("pi3", Fraction(1), ">", Fraction(5))
    tex = render_equation(resp["parameters"], "pi3", Fraction(1), ">", Fraction(5))
    assert tex.startswith("1 - 5\\pi_{3}^{-1} = \\int_0^1")
    assert "\\sqrt{3}" in tex and "\\pi \\sqrt[3]{1 - x^{3}}" in tex
    resp = prove("pi3_u", Fraction(1), "<", Fraction(5, 7))
    tex = render_equation(resp["parameters"], "pi3_u", Fraction(1), "<", Fraction(5, 7))
    assert tex.startswith("\\dfrac{5}{7}U^{-1} - 1 = \\int_0^1")


def test_render_equation_fallback():
    resp = prove("pi3", Fraction(1), "<", Fraction(100))
    tex = render_equation(resp["parameters"], "pi3", Fraction(1), "<", Fraction(100))
    assert tex.endswith(" > 0")
    assert " + " not in tex.split("\\mathrm{d} x")[0]  # the +b sits after dx
    assert "\\mathrm{d} x+" in tex


@pytest.mark.parametrize("kind", KINDS)
def test_verify_dispatch(monkeypatch, kind):
    """exact_check.verify routes via solve.FAMILY (merge wiring is the
    leader's); setitem on the shared dict exercises the real dispatch."""
    monkeypatch.setitem(solver.FAMILY, kind, "dixon")
    bound = {"pi3": "5", "pi3_u": "2/3", "pi3_a": "18/5"}[kind]
    comp = {"pi3": ">", "pi3_u": ">", "pi3_a": ">"}[kind]
    resp = prove(kind, Fraction(1), comp, Fraction(bound))
    res = verify(kind, Fraction(1), comp, Fraction(bound), resp["parameters"])
    assert res["identity_ok"]
    assert res["nonneg"]


TRUE_BOUNDS = {
    "pi3": ((">", "5"), ("<", "27/5")),
    "pi3_u": ((">", "2/3"), ("<", "7/10")),
    "pi3_a": ((">", "18/5"), ("<", "11/3")),
}
FALSE_BOUNDS = {
    "pi3": (">", "27/5"),
    "pi3_u": ("<", "2/3"),
    "pi3_a": ("<", "18/5"),
}


@pytest.mark.parametrize("kind", KINDS)
def test_full_pipeline(kind):
    """Through the registered family, solver.prove(exact=True) must emit
    self-checking proofs — and certified_cmp must reject false claims with
    WrongDirection."""
    for comp, bound in TRUE_BOUNDS[kind]:
        resp = solver.prove(kind, "1", comp, bound, exact=True)
        res = verify_response(kind, Fraction(1), comp, Fraction(bound), resp)
        assert res["identity_ok"]
        assert res["nonneg"]
    comp, bound = FALSE_BOUNDS[kind]
    with pytest.raises(WrongDirection):
        solver.prove(kind, "1", comp, bound, exact=True)
