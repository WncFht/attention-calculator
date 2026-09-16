"""Exact-mode tests for the zeta_odd family: zeta5, zeta7 (new W3 types).

Drives kernels.zeta_odd.prove + exact_check.zeta_odd.check directly — the
types are not yet wired into solve.FAMILY/TYPES and have no site counterpart.
Bounds are CF convergents of ζ(5)/ζ(7) chosen with mpmath:

    ζ(5) ≈ 1.03692775514336992633   ζ(7) ≈ 1.00834927738192282684

Every emitted proof must pass check() exactly. False claims raise
WrongDirection mid-scan: the identity ∫f = ±(C−r) holds exactly for every
solved candidate, so a false claim's first sign-definite P is uniformly
non-positive (observed on all probes below — never a passing fake proof,
never silent exhaustion). A true claim beyond the (m,n) ≤ 10 budget
exhausts to NoSolution.
"""

from fractions import Fraction

import pytest

from attention_calculator.engine import NoSolution, WrongDirection
from attention_calculator.exact_check.zeta_odd import check
from attention_calculator.kernels import zeta_odd


def run(kind, power, comp, bound):
    resp = zeta_odd.prove(kind, Fraction(power), comp, Fraction(bound))
    return check(kind, Fraction(power), comp, Fraction(bound), resp["parameters"])


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        # ζ(5) '<' upper bounds, loose -> 7e-9 tight
        ("zeta5", "1", "<", "2"),
        ("zeta5", "1", "<", "28/27"),  # err +1.1e-4 -> (2,2)
        ("zeta5", "1", "<", "365/352"),  # err +4.1e-6 -> (4,3)
        ("zeta5", "1", "<", "10895/10507"),  # err +7.3e-9 -> (6,7)
        # ζ(5) '>' lower bounds
        ("zeta5", "1", ">", "1"),
        ("zeta5", "1", ">", "337/325"),  # err -4.7e-6 -> (3,4)
        ("zeta5", "1", ">", "702/677"),  # err -1.3e-7 -> (5,5)
        ("zeta5", "1", ">", "11597/11184"),  # err -1.2e-9 -> (7,7)
        # ζ(7) '<' upper bounds
        ("zeta7", "1", "<", "2"),
        ("zeta7", "1", "<", "120/119"),  # err +5.4e-5 -> (2,0)
        ("zeta7", "1", "<", "483/479"),  # err +1.5e-6 -> (4,1)
        ("zeta7", "1", "<", "1570/1557"),  # err +1.1e-7 -> (6,1)
        ("zeta7", "1", "<", "466417/462555"),  # err +1.3e-12 -> (8,9)
        # ζ(7) '>' lower bounds
        ("zeta7", "1", ">", "1"),
        ("zeta7", "1", ">", "121/120"),  # err -1.6e-5 -> (3,0)
        ("zeta7", "1", ">", "1087/1078"),  # err -4.8e-7 -> (5,1)
        ("zeta7", "1", ">", "4227/4192"),  # err -4.1e-8 -> (5,3)
        # power as coefficient: 2ζ(5) ≈ 2.0738555, ζ(5)/2 ≈ 0.5184639,
        # 3ζ(7) ≈ 3.0250478, ζ(7)/2 ≈ 0.5041746
        ("zeta5", "2", ">", "2"),
        ("zeta5", "2", "<", "21/10"),
        ("zeta5", "2", "<", "365/176"),  # 2·(365/352): err +8.1e-6
        ("zeta5", "1/2", ">", "1/2"),
        ("zeta5", "1/2", "<", "13/25"),
        ("zeta5", "3", ">", "3"),
        ("zeta5", "3", "<", "16/5"),
        ("zeta7", "3", ">", "3"),
        ("zeta7", "3", "<", "61/20"),
        ("zeta7", "3", "<", "1570/519"),  # 3·(1570/1557): err +3.4e-7
        ("zeta7", "2", ">", "2"),
        ("zeta7", "1/2", ">", "1/2"),
        ("zeta7", "1/2", "<", "51/100"),
        # ζ(9) ≈ 1.00200839282608221442, ζ(11) ≈ 1.00049418860411946456
        ("zeta9", "1", ">", "1"),
        ("zeta9", "1", ">", "499/498"),  # err -3.6e-7
        ("zeta9", "1", "<", "2"),
        ("zeta9", "1", "<", "5488/5477"),  # err +5.9e-9
        ("zeta9", "3", "<", "4"),
        ("zeta9", "1/2", ">", "1/2"),
        ("zeta11", "1", ">", "1"),
        ("zeta11", "1", "<", "101/100"),
        ("zeta11", "1", "<", "4049/4047"),  # err +4.6e-9
        ("zeta11", "2", ">", "2"),
        # negative coefficient flips the inequality's constant side only
        ("zeta5", "-1", "<", "-1"),
        # degenerate power=0: true constant identity, like zeta3's 0-vs-r
        ("zeta5", "0", "<", "1/2"),
        ("zeta7", "0", ">", "-1"),
    ],
)
def test_emitted_proofs_verify(kind, power, comp, bound):
    res = run(kind, power, comp, bound)
    assert res["identity_ok"]
    assert res["nonneg"]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "bound"),
    [
        ("zeta5", "1", ">", "28/27"),  # 28/27 > ζ(5)
        ("zeta5", "1", "<", "337/325"),  # 337/325 < ζ(5)
        ("zeta7", "1", ">", "5797/5749"),  # 5797/5749 > ζ(7) by 7.5e-10
        ("zeta7", "1", "<", "121/120"),  # 121/120 < ζ(7)
        ("zeta7", "2", "<", "2"),  # 2 < 2ζ(7)
    ],
)
def test_false_claim_raises_wrong_direction(kind, power, comp, bound):
    # prove() itself does no direction screening; the mid-scan WrongDirection
    # is the exact-moment certificate that the flipped inequality holds
    with pytest.raises(WrongDirection):
        zeta_odd.prove(kind, Fraction(power), comp, Fraction(bound))


def test_true_claim_beyond_budget_reports_no_solution():
    # ζ(5) > 229834/221649 is true (err 1e-12) but needs more than (m,n) ≤ 10
    with pytest.raises(NoSolution):
        zeta_odd.prove("zeta5", Fraction(1), ">", Fraction(229834, 221649))


def test_corrupted_params_fail_identity():
    resp = zeta_odd.prove("zeta5", Fraction(1), "<", Fraction(28, 27))
    params = dict(resp["parameters"], au_val=str(int(resp["parameters"]["au_val"]) + 1))
    res = check("zeta5", Fraction(1), "<", Fraction(28, 27), params)
    assert not res["identity_ok"]


def test_proof_does_not_verify_flipped_claim():
    resp = zeta_odd.prove("zeta7", Fraction(1), "<", Fraction(120, 119))
    res = check("zeta7", Fraction(1), ">", Fraction(120, 119), resp["parameters"])
    assert not res["identity_ok"]


def test_render_equation_smoke():
    resp = zeta_odd.prove("zeta5", Fraction(1), "<", Fraction(28, 27))
    tex = zeta_odd.render_equation(resp["parameters"], "zeta5", Fraction(1), "<", Fraction(28, 27))
    assert tex.startswith("\\dfrac{28}{27} - \\zeta(5) = \\int_0^1 ")
    assert tex.endswith("\\ln^4(x) \\mathrm{d} x > 0")
