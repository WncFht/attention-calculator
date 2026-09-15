"""Endpoint shape tests for the Flask app.

Kernel families and the decompose module are stubbed: /calculate is patched at
solve.prove, /get_integral_image injects a fake family module into sys.modules,
/decompose_inequality injects a fake attention_calculator.decompose.

Statuses/messages follow the probed site contract: format errors -> 400,
domain and search failures -> 404, render failures -> 500 generic text.
"""

import sys
import types
from fractions import Fraction

import pytest

import attention_calculator
from attention_calculator import engine, server, solve

PARAMS = {
    "m": 3,
    "n": 3,
    "a_val": "47/120",
    "b_val": "-13/120",
    "c_val": "0",
    "au_val": 47,
    "bu_val": -13,
    "cu_val": 0,
    "u_val": 120,
}


@pytest.fixture
def client():
    return server.app.test_client()


def post_calc(client, **form):
    """POST /calculate with the given urlencoded fields."""
    return client.post("/calculate", data=form)


def test_index_served(client):
    """GET / serves the calculator page."""
    resp = client.get("/")
    assert resp.status_code == 200
    assert "注意力计算器" in resp.get_data(as_text=True)


def test_routes_and_error_envelope(client, monkeypatch):
    """Site-surface details: JSON error bodies, /en/ 404, /attention mount."""
    def fake(*a):
        return {"parameters": dict(PARAMS), "solution": "a = 47/120, b = -13/120"}

    monkeypatch.setattr(solve, "prove", fake)
    # absent comparison defaults to '>' (present-but-empty stays 400)
    resp = post_calc(client, type="pi", power="1", rational="22/7")
    assert resp.status_code == 200

    assert client.get("/en").status_code == 200
    resp = client.get("/en/")
    assert resp.status_code == 404
    assert resp.get_json() == {"error": "请求的页面不存在"}
    for path in ("/attention", "/attention/", "/attention/en"):
        assert client.get(path).status_code == 200, path
    resp = client.get("/attention/static/title.png")
    assert resp.status_code == 200

    # every non-2xx is JSON: unknown path 404, wrong method -> 500
    assert client.get("/no-such-path").get_json() == {"error": "请求的页面不存在"}
    for method, path in [("get", "/calculate"), ("post", "/get_integral_image"),
                         ("post", "/"), ("post", "/en")]:
        resp = getattr(client, method)(path)
        assert resp.status_code == 500, (method, path)
        assert resp.get_json() == {"error": INTERNAL_ERROR}

    resp = client.get("/favicon.ico")
    assert resp.status_code == 204 and resp.get_data() == b""


def test_calculate_success(client, monkeypatch):
    """Successful prove -> site envelope; type echoes the request verbatim."""
    def fake(*a):
        return {"parameters": dict(PARAMS), "solution": "a = 47/120, b = -13/120"}

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison="<", rational="22/7")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["success"] is True
    assert data["type"] == "pi"
    # site wire types: m/n int, au/bu/cu/u_val strings, unified_form always present
    assert data["parameters"]["m"] == 3
    assert data["parameters"]["a_val"] == "47/120"
    assert data["parameters"]["au_val"] == "47"
    assert data["parameters"]["u_val"] == "120"
    assert data["parameters"]["unified_form"] == {}
    assert data["equations"]["solution"] == "a = 47/120, b = -13/120"


def test_calculate_wrong_direction(client, monkeypatch):
    """WrongDirection -> 404 with the site's '方向反了' message."""
    def fake(*a):
        raise engine.WrongDirection

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison=">", rational="22/7")
    assert resp.status_code == 404
    assert resp.get_json() == {"error": "要证明的式子不等号方向反了"}


def test_calculate_no_solution_limits(client, monkeypatch):
    """NoSolution -> 404, message carries the per-type exponent budget."""
    def fake(*a):
        raise engine.NoSolution

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison="<", rational="1/1")
    assert resp.status_code == 404
    assert resp.get_json()["error"] == "在指数不超过30的范围内未找到<方向的解"
    resp = post_calc(client, type="sin_q", power="1", comparison=">", rational="1/1")
    assert resp.get_json()["error"] == "在指数不超过10的范围内未找到>方向的解"


BAD_TYPE = "无效的证明类型"
BAD_COMP = "无效的不等号方向"
BAD_POWER = "左侧系数格式无效"
BAD_POWER_DENOM = "左侧系数分母不能为0"
BAD_POWER_CAP = "左侧系数请输入小于10^16的整数或分数"
BAD_RATIONAL = "右侧有理数格式无效"
ZERO_DENOM = "右侧有理数分母不能为0"
RATIONAL_TOO_BIG = "右侧有理数请输入小于10^16的整数或分数"
INTERNAL_ERROR = "服务器内部错误，请稍后再试"


def test_calculate_format_errors(client):
    """Malformed input -> 400 with the site's probed texts."""
    def form(power="1", rational="1/2", **kw):
        return {"type": "pi", "power": power, "comparison": ">", "rational": rational, **kw}

    cases = [
        (form(type="nope"), BAD_TYPE),
        (form(comparison="="), BAD_COMP),
        (form(power=""), BAD_POWER),
        (form(power="x"), BAD_POWER),
        (form(power="-1"), BAD_POWER),
        (form(power="1/0"), BAD_POWER_DENOM),
        (form(power="1" + "0" * 17), BAD_POWER_CAP),
        (form(rational=""), BAD_RATIONAL),
        (form(rational="abc"), BAD_RATIONAL),
        (form(rational="3.14"), BAD_RATIONAL),
        (form(rational="-1/10"), BAD_RATIONAL),
        (form(rational="1/-2"), BAD_RATIONAL),
        (form(rational="1/0"), ZERO_DENOM),
        (form(rational="1" + "0" * 17 + "/1"), RATIONAL_TOO_BIG),
        (form(rational="1/" + "1" + "0" * 16), RATIONAL_TOO_BIG),
        # 站端把整个右侧字段校验置于左侧之前（成对探针实测）
        (form(power="abc", rational="xyz"), BAD_RATIONAL),
        (form(power="abc", rational="1/0"), ZERO_DENOM),
        (form(power="1/0", rational="1/0"), ZERO_DENOM),
        (form(power="1" + "0" * 17, rational="1/0"), ZERO_DENOM),
        (form(power="1/0", rational="1" + "0" * 17), RATIONAL_TOO_BIG),
        # type 缺省回退 pi；显式空串仍是无效类型
        ({"power": "1", "comparison": "<", "rational": "22/7"}, None),
        (form(type=""), BAD_TYPE),
    ]
    for form, error in cases:
        resp = post_calc(client, **form)
        if error is None:
            assert resp.status_code == 200 and resp.get_json()["type"] == "pi"
            continue
        assert resp.status_code == 400, (form, resp.get_json())
        assert resp.get_json() == {"error": error}, form


def test_calculate_domain_errors(client):
    """Out-of-domain parameters -> 404 with the site's probed texts."""
    pi_frac_msg = "请在输入一个在(0,1/2)内的分数，本情况不支持整数"
    cases = [
        ({"type": "ln_q", "power": "1"}, "请在ln后输入一个大于1的数"),
        ({"type": "ln_q", "power": "1/2"}, "请在ln后输入一个大于1的数"),
        ({"type": "ln_q_square", "power": "1"}, "请在ln后输入一个大于1的数"),
        ({"type": "sin_q", "power": "4"}, "请在sin后输入一个在(0,π)内的数"),
        ({"type": "cos_q", "power": "4"}, "请在cos后输入一个在(0,π/2)内的数"),
        ({"type": "tan_q", "power": "2"}, "请在tan后输入一个在(0,π/2)内的数"),
        ({"type": "cot_q", "power": "2"}, "请在cot后输入一个在(0,π/2)内的数"),
        ({"type": "sin_pi_q", "power": "1/2"}, pi_frac_msg),
        ({"type": "sin_pi_q", "power": "1"}, pi_frac_msg),
        ({"type": "cos_pi_q", "power": "1"}, pi_frac_msg),
    ]
    for form, error in cases:
        resp = post_calc(client, comparison=">", rational="1/2", **form)
        assert resp.status_code == 404, (form, resp.get_json())
        assert resp.get_json() == {"error": error}, form


def test_calculate_kernel_value_error(client, monkeypatch):
    """Kernel-side ValueError -> 400 with the kernel's message."""
    def fake(*a):
        raise ValueError("左侧系数格式无效")

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="e", power="0", comparison=">", rational="8/3")
    assert resp.status_code == 400
    assert resp.get_json()["error"] == "左侧系数格式无效"


def test_get_integral_image(client, monkeypatch):
    """Query params are coerced and forwarded to the family render_equation."""
    seen = {}

    def fake_render(params, kind, power, comp, bound):
        seen.update(
            {"params": params, "kind": kind, "power": power, "comp": comp, "bound": bound}
        )
        return "\\dfrac{22}{7} - \\pi = \\int_0^1 f(x) \\mathrm{d} x > 0"

    fake = types.ModuleType("attention_calculator.kernels.quadlog")
    fake.render_equation = fake_render
    monkeypatch.setitem(sys.modules, "attention_calculator.kernels.quadlog", fake)

    resp = client.get(
        "/get_integral_image",
        query_string={
            "type": "pi",
            "comparison": "<",
            "rational": "\\frac{22}{7}",
            "coef": "1",
            **{k: str(PARAMS[k]) for k in server.IMG_INT_KEYS + server.IMG_FRAC_KEYS},
        },
    )
    assert resp.status_code == 200
    assert resp.get_json() == {
        "equation": "\\dfrac{22}{7} - \\pi = \\int_0^1 f(x) \\mathrm{d} x > 0"
    }
    assert seen["kind"] == "pi"
    assert seen["params"]["m"] == 3
    assert seen["params"]["u_val"] == 120
    assert seen["params"]["a_val"] == Fraction(47, 120)
    # coef/rational 原文传给渲染层（站端原样回显，不约分）
    assert seen["power"] == "1"
    assert seen["bound"] == "\\frac{22}{7}"


def test_get_integral_image_echoes_raw(client):
    """The site echoes request digits unreduced, even the \\frac macro verbatim."""
    base = {
        "type": "pi",
        "comparison": "<",
        "coef": "1",
        **{k: str(PARAMS[k]) for k in server.IMG_INT_KEYS + server.IMG_FRAC_KEYS},
    }
    resp = client.get("/get_integral_image",
                      query_string={**base, "rational": "3140/1000"})
    assert "\\dfrac{3140}{1000}" in resp.get_json()["equation"]

    resp = client.get("/get_integral_image",
                      query_string={**base, "rational": "\\frac{3140}{1000}"})
    eq = resp.get_json()["equation"]
    assert "\\frac{3140}{1000}" in eq and "dfrac{3140" not in eq

    resp = client.get("/get_integral_image",
                      query_string={**base, "rational": "3",
                                    "coef": "2/4", "comparison": ">"})
    assert "\\dfrac{2}{4}\\pi" in resp.get_json()["equation"]


def test_get_integral_image_errors(client, monkeypatch):
    """Ordered per-field 400s like the site; render failures -> generic 500."""
    base = {
        "type": "pi",
        "comparison": "<",
        "rational": "22/7",
        "coef": "1",
        **{k: str(PARAMS[k]) for k in server.IMG_INT_KEYS + server.IMG_FRAC_KEYS},
    }

    def get(**over):
        q = dict(base)
        for k, v in over.items():
            if v is None:
                q.pop(k)
            else:
                q[k] = v
        return client.get("/get_integral_image", query_string=q)

    # type/comparison first: absent defaults, empty/invalid 400
    assert get(type=None).status_code == 200
    assert get(comparison=None).status_code == 200
    for over, err in [
        ({"type": "bogus"}, "无效的证明类型"),
        ({"type": ""}, "无效的证明类型"),
        ({"comparison": "="}, "无效的不等号方向"),
        ({"comparison": ""}, "无效的不等号方向"),
    ]:
        resp = get(**over)
        assert resp.status_code == 400 and resp.get_json() == {"error": err}, over

    # integer fields in probed order; m,n bounded [0,30], u_val >= 1,
    # au/bu/cu unbounded (negatives and huge values pass)
    for over, err in [
        ({"m": "abc"}, "m必须是整数"),
        ({"m": "3.5"}, "m必须是整数"),
        ({"m": "3/2"}, "m必须是整数"),
        ({"m": "+3"}, "m必须是整数"),  # int() 会吃 '+'，站端 strip+regex 拒绝
        ({"m": "3 0"}, "m必须是整数"),  # 内部空白拒绝
        ({"m": None}, "m必须是整数"),
        ({"n": None}, "n必须是整数"),
        ({"u_val": "x"}, "u_val必须是整数"),
        ({"m": "-1"}, "m过小"),
        ({"m": "31"}, "m过大"),
        ({"n": "-1"}, "n过小"),
        ({"n": "31"}, "n过大"),
        ({"u_val": "0"}, "u_val过小"),
        ({"m": "x", "a_val": "zz"}, "m必须是整数"),  # 整数组先于分数组
        ({"n": "x", "m": "y"}, "m必须是整数"),  # m 先于 n
    ]:
        resp = get(**over)
        assert resp.status_code == 400 and resp.get_json() == {"error": err}, over
    for ok in ({"m": "30"}, {"n": "30"}, {"m": " 3 "}, {"m": "\t3\n"},
               {"au_val": "-5"}, {"cu_val": "9"}, {"u_val": "10" + "0" * 16}):
        assert get(**ok).status_code == 200, ok

    # fraction fields: format + denominator-0, no numerator cap
    for over, err in [
        ({"a_val": "x"}, "a_val格式无效"),
        ({"a_val": "1/0"}, "a_val分母不能为0"),
        ({"b_val": None}, "b_val格式无效"),
        ({"a_val": "x", "b_val": "y"}, "a_val格式无效"),  # a 先于 b
    ]:
        resp = get(**over)
        assert resp.status_code == 400 and resp.get_json() == {"error": err}, over
    assert get(a_val="-1/2").status_code == 200
    assert get(b_val="9" * 20 + "/7").status_code == 200

    # coef/rational unvalidated, echoed raw (probed: coef=x -> x\pi,
    # rational absent -> "0 - \pi")
    eq = get(coef="x").get_json()["equation"]
    assert "x\\pi" in eq
    eq = get(rational=None).get_json()["equation"]
    assert eq.startswith("0 - \\pi")

    # a kernel crash still lands in the site's generic 500
    monkeypatch.setitem(server.render.FAMILY, "pi", "no_such_family")
    resp = get()
    assert resp.status_code == 500
    assert resp.get_json() == {"error": INTERNAL_ERROR}


def fake_decompose_module(payload):
    """Build a stub module exposing decompose_inequality(problem)."""
    mod = types.ModuleType("attention_calculator.decompose")
    mod.decompose_inequality = lambda problem: payload
    return mod


def install_decompose(monkeypatch, mod):
    """Register a fake attention_calculator.decompose for the lazy import."""
    monkeypatch.setitem(sys.modules, "attention_calculator.decompose", mod)
    monkeypatch.setattr(attention_calculator, "decompose", mod, raising=False)


def test_decompose_inequality(client, monkeypatch):
    """The decompose module's dict is passed through verbatim."""
    payload = {
        "success": True,
        "problem": "pi^2+8*pi>35",
        "normalized_latex": "\\pi^{2}+8\\pi>35",
        "direct_basic": False,
        "basic_count": 2,
        "decomposition_latex": "\\pi^{2}>\\dfrac{227}{23}",
        "steps": [{"type": "pi_n", "equation": "..."}],
    }
    install_decompose(monkeypatch, fake_decompose_module(payload))
    resp = client.post("/decompose_inequality", data={"problem": "pi^2+8*pi>35"})
    assert resp.status_code == 200
    assert resp.get_json() == payload


def test_decompose_inequality_errors(client, monkeypatch):
    """Empty problem -> 400; ValueError from the module -> 400 error body."""
    resp = client.post("/decompose_inequality", data={"problem": "  "})
    assert resp.status_code == 400
    assert resp.get_json() == {"error": "请输入一个只包含一个 > 或 < 的不等式"}

    def boom(problem):
        raise ValueError("无法解析该不等式")

    mod = types.ModuleType("attention_calculator.decompose")
    mod.decompose_inequality = boom
    install_decompose(monkeypatch, mod)
    resp = client.post("/decompose_inequality", data={"problem": "???"})
    assert resp.status_code == 400
    assert resp.get_json() == {"error": "无法解析该不等式"}
