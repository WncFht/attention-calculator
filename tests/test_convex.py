"""Tests for the /convex sibling app clone (POST /convex/prove).

Deterministic parts pinned against docs/sibling-apps.md §1 samples:
parser error strings (ast.dump verbatim), normalization text/latex,
curvature labels, status/reason texts, response key shape, and the
documented example's proof fields. Numeric minima are asserted
approximately — exact float digits await the probe agent's pinning.
"""

from pathlib import Path

import pytest

from attention_calculator import server

MISSING = "请输入一个不等式。"
REASON_INCONCLUSIVE = (
    "整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。"
)
REASON_FAILED = "数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。"
REASON_NO_LINE = "不等式数值上已通过，但当前情形没有生成中间直线证明。"
REASON_SEARCH_MISS = "不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。"

SAMPLE = "exp(x)-log(x)-261/112>0"


@pytest.fixture
def client():
    return server.app.test_client()


def prove(client, **form):
    """POST /convex/prove with urlencoded fields."""
    return client.post("/convex/prove", data=form)


def test_convex_page_and_routes(client):
    """GET /convex/ 与 /convex/en 返回 site-convex.html 原字节；错误带 ok:false。"""
    expected = Path("bench/data/site-convex.html").read_bytes()
    for path in ("/convex/", "/convex/en"):
        resp = client.get(path)
        assert resp.status_code == 200
        assert resp.get_data() == expected, path
    assert client.get("/convex/static/bg1.png").status_code == 200

    # 姊妹应用统一 {"error", "ok":false} 包络；域内未匹配路径与错方法
    # 一律 500（probe: /health/xyz、/convex/Prove 实测站端行为）
    resp = client.get("/convex/nope")
    assert resp.status_code == 500
    assert resp.get_json() == {"error": "服务器内部错误，请稍后再试。", "ok": False}
    resp = client.get("/convex/prove")
    assert resp.status_code == 500
    assert resp.get_json()["ok"] is False


def test_missing_inequality(client):
    """缺 inequality 字段 -> 400 原文案。"""
    resp = client.post("/convex/prove", data={})
    assert resp.status_code == 400
    assert resp.get_json() == {"error": MISSING, "ok": False}


def test_parser_error_surface(client):
    """ast.dump 泄漏的错误串逐字复现。"""
    resp = prove(client, inequality="foo(x)>0")
    assert resp.status_code == 400
    assert resp.get_json()["error"] == (
        "unsupported atom Call(func=Name(id='foo', ctx=Load()), "
        "args=[Name(id='x', ctx=Load())], keywords=[])"
    )
    resp = prove(client, inequality="x+y>0")
    assert resp.get_json()["error"] == "unsupported atom Name(id='y', ctx=Load())"
    # log(17/30) 在 line 字段上触发实测错误
    resp = prove(client, inequality=SAMPLE, line="log(17/30)")
    assert resp.get_json()["error"] == "log only supports argument x"


def test_caret_equivalent(client):
    """`^` 与 `**` 响应字节相同。"""
    a = prove(client, inequality="x^2+1>log(x)").get_data()
    b = prove(client, inequality="x**2+1>log(x)").get_data()
    assert a == b


def test_sample_proved_full_shape(client):
    """docs/sibling-apps.md §1 的完整样例响应：确定性字段逐字，数值近似。"""
    resp = prove(client, inequality=SAMPLE)
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["ok"] is True
    r = data["result"]
    assert set(r) == {
        "curvature",
        "minimum",
        "normalized",
        "ok",
        "proof",
        "provided_line",
        "reason",
        "status",
    }
    assert r["curvature"] == {"difference": "convex", "left": "convex", "right": "concave"}
    assert r["normalized"] == {
        "difference": "e^x - ln x - 261/112",
        "difference_latex": "e^x - \\ln x - \\frac{261}{112}",
        "left": "e^x",
        "left_latex": "e^x",
        "right": "ln x + 261/112",
        "right_latex": "\\ln x + \\frac{261}{112}",
    }
    assert r["status"] == "proved" and r["ok"] is True and r["reason"] == ""
    assert r["provided_line"] is None

    m = r["minimum"]
    assert set(m) == {"value", "value_text", "x", "x_text"}
    assert m["x"] == pytest.approx(0.567143290409784, abs=1e-9)
    assert m["value"] == pytest.approx(8.981904537908036e-06, abs=1e-12)
    assert m["x_text"] == f"{m['x']:.12g}" and m["value_text"] == f"{m['value']:.12g}"

    p = r["proof"]
    assert set(p) == {
        "formula_latex",
        "left_gap_min",
        "left_gap_min_x",
        "line_latex",
        "line_text",
        "tangent_at",
        "tangent_at_latex",
    }
    assert p["tangent_at"] == "17/30"
    assert p["tangent_at_latex"] == "\\frac{17}{30}"
    assert p["line_text"] == "(30/17*x - 1 + ln(17/30)) + 261/112"
    assert p["line_latex"] == ("(\\frac{30}{17}x - 1 + \\ln \\frac{17}{30}) + \\frac{261}{112}")
    assert p["formula_latex"] == (
        "e^x > (\\frac{30}{17}x - 1 + \\ln \\frac{17}{30}) + \\frac{261}{112} "
        "\\ge \\ln x + \\frac{261}{112}"
    )
    assert p["left_gap_min"] == pytest.approx(8.004855962417956e-06, abs=1e-10)
    assert p["left_gap_min_x"] == pytest.approx(0.5679840376059393, abs=1e-9)


def test_lt_direction_normalized(client):
    """`<`/`<=` 取反归一，与 `>` 形式的 normalized 完全一致。"""
    r = prove(client, inequality="log(x)+261/112<exp(x)").get_json()["result"]
    assert r["normalized"]["difference"] == "e^x - ln x - 261/112"
    assert r["normalized"]["right"] == "ln x + 261/112"


def test_normalize_rules(client):
    """正常数留左、负常数移右；负系数凹项取负移右。"""
    r = prove(client, inequality="x^2+1>log(x)").get_json()["result"]
    assert r["normalized"]["left"] == "x^2 + 1"
    assert r["normalized"]["right"] == "ln x"
    assert r["normalized"]["difference"] == "x^2 + 1 - ln x"


def test_status_inconclusive(client):
    """左侧凹 -> inconclusive，minimum/proof 为 null。"""
    r = prove(client, inequality="log(x)>0").get_json()["result"]
    assert r["status"] == "inconclusive" and r["ok"] is False
    assert r["reason"] == REASON_INCONCLUSIVE
    assert r["minimum"] is None and r["proof"] is None
    assert r["curvature"]["left"] == "concave"


def test_status_failed(client):
    """数值最小值为负 -> failed，minimum 仍携带。"""
    r = prove(client, inequality="exp(x)-log(x)-3>0").get_json()["result"]
    assert r["status"] == "failed" and r["ok"] is False
    assert r["reason"] == REASON_FAILED
    assert r["minimum"]["value"] < 0 and r["proof"] is None


def test_proved_no_line_when_left_affine(client):
    """left 仿射跳过搜索 -> proved + REASON_NO_LINE + proof null。"""
    r = prove(client, inequality="x>log(x)").get_json()["result"]
    assert r["status"] == "proved" and r["ok"] is True
    assert r["proof"] is None and r["reason"] == REASON_NO_LINE
    assert r["curvature"]["left"] == "affine"


def test_proved_search_miss_zero_gap(client):
    """e^x>=x+1: 切线与右侧重合零 gap 被拒 -> 没找到漂亮切点（实测怪癖）。"""
    r = prove(client, inequality="exp(x)>=x+1").get_json()["result"]
    assert r["status"] == "proved" and r["ok"] is True
    assert r["proof"] is None and r["reason"] == REASON_SEARCH_MISS


def test_provided_line(client):
    """line=1 对 x^2+1>0 验证通过但 proof 仍为 null（不回灌）。"""
    r = prove(client, inequality="x^2+1>0", line="1").get_json()["result"]
    assert r["status"] == "proved"
    assert r["proof"] is None
    pl = r["provided_line"]
    assert set(pl) == {
        "b",
        "left_gap_min",
        "left_gap_min_x",
        "m",
        "ok",
        "right_gap_min",
        "right_gap_min_x",
    }
    assert pl["m"] == 0.0 and pl["b"] == 1.0 and pl["ok"] is True
    assert pl["left_gap_min"] == pytest.approx(0.0, abs=1e-9)
    assert pl["right_gap_min"] == pytest.approx(1.0, abs=1e-9)


def test_line_parse_errors_400(client):
    """line 惰性解析：只在进入证明路径时解析（probe: x>0+foo(x) 边界失败→200；
    x^2+1>0+foo(x) 可证路径→400 unsupported atom）。"""
    resp = prove(client, inequality="x>0", line="foo(x)")
    assert resp.status_code == 200
    resp = prove(client, inequality="x^2+1>0", line="foo(x)")
    assert resp.status_code == 400
    assert "unsupported atom" in resp.get_json()["error"]
