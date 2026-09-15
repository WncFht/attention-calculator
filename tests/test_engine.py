"""Engine primitive tests, anchored on the article's printed moment tables."""

from fractions import Fraction as F

import pytest

from attention_calculator.engine import (
    gauss_solve,
    mn_order,
    poly_nonneg,
    poly_nonpos,
    solve_moment,
)


def test_mn_order_shape():
    seq = list(mn_order(3))
    # 站点约定：同和、同 |m-n| 时 m 小者先（pi>8/3 实测取 (0,1) 而非 (1,0)）
    assert seq[:4] == [(0, 0), (0, 1), (1, 0), (1, 1)]
    # m+n non-decreasing
    assert [a + b for a, b in seq] == sorted(a + b for a, b in seq)
    # within each sum, |m-n| ascending
    by_sum = {}
    for m, n in seq:
        by_sum.setdefault(m + n, []).append(abs(m - n))
    assert all(v == sorted(v) for v in by_sum.values())


def test_gauss_solve():
    assert gauss_solve(
        [[F(1), F(2)], [F(3), F(-1)]], [F(5), F(1)]
    ) == [F(1), F(2)]
    with pytest.raises(ValueError):
        gauss_solve([[F(1), F(1)], [F(2), F(2)]], [F(1), F(2)])


def test_poly_nonneg_rules():
    assert poly_nonneg([F(10), F(-16), F(10)])  # vertex inside, disc>0
    assert not poly_nonneg([F(1), F(-5), F(1)])  # dips below 0
    assert poly_nonneg([F(3), F(-2)])  # 3-2x >= 1 > 0
    assert not poly_nonneg([F(1), F(-3)])  # 1-3x < 0 near 1
    assert poly_nonpos([F(-1), F(-1)])  # -1-x <= 0


def test_solve_moment_article_vector():
    """Article type-1 n=2: solve 8π>25 -> a=10, b=-16, c=10."""

    def m2(a, b, c):
        return {
            "pi": -b / F(2),
            "ln2": a - c,
            "1": -2 * a / F(3) + 19 * b / F(12) + 7 * c / F(10),
        }

    basis = [m2(F(1), F(0), F(0)), m2(F(0), F(1), F(0)), m2(F(0), F(0), F(1))]
    assert solve_moment(basis, {"pi": F(8), "ln2": F(0), "1": F(-25)}) == [
        F(10),
        F(-16),
        F(10),
    ]
