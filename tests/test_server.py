"""Endpoint shape tests for the Flask app.

Kernel families and the decompose module are stubbed: /calculate is patched at
solve.prove, /get_integral_image injects a fake family module into sys.modules,
/decompose_inequality injects a fake attention_calculator.decompose.
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
    "au_val": "47",
    "bu_val": "-13",
    "cu_val": "0",
    "u_val": "120",
    "unified_form": {},
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


def test_calculate_success(client, monkeypatch):
    """Successful prove -> site success envelope with parameters + equations."""
    def fake(*a):
        return {"parameters": dict(PARAMS), "solution": "a = 47/120, b = -13/120"}

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison="<", rational="22/7")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["success"] is True
    assert data["type"] == "pi"
    assert data["parameters"]["m"] == 3
    assert data["parameters"]["a_val"] == "47/120"
    assert data["equations"]["solution"] == "a = 47/120, b = -13/120"


def test_calculate_wrong_direction(client, monkeypatch):
    """WrongDirection -> the site's '方向反了' message."""
    def fake(*a):
        raise engine.WrongDirection

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison=">", rational="22/7")
    assert resp.get_json() == {"success": False, "error": "要证明的式子不等号方向反了"}


def test_calculate_no_solution_limits(client, monkeypatch):
    """NoSolution message carries the per-type exponent budget."""
    def fake(*a):
        raise engine.NoSolution

    monkeypatch.setattr(solve, "prove", fake)
    resp = post_calc(client, type="pi", power="1", comparison="<", rational="1/1")
    assert resp.get_json()["error"] == "在指数不超过30的范围内未找到<方向的解"
    resp = post_calc(client, type="sin_q", power="1", comparison=">", rational="1/1")
    assert resp.get_json()["error"] == "在指数不超过10的范围内未找到>方向的解"


def test_calculate_validation(client):
    """Malformed input -> 400 with the original page's Chinese prompts."""
    resp = post_calc(client, type="nope", power="1", comparison=">", rational="1/2")
    assert resp.status_code == 400
    assert resp.get_json()["success"] is False

    resp = post_calc(client, type="sin_q", power="", comparison=">", rational="1/2")
    assert resp.get_json()["error"] == "请输入sin后的值"

    resp = post_calc(client, type="pi", power="abc", comparison=">", rational="1/2")
    assert resp.get_json()["error"] == "请输入π的系数"

    resp = post_calc(client, type="pi_n", power="11", comparison=">", rational="1/2")
    assert resp.get_json()["error"] == "π的次数不要超过10"

    resp = post_calc(client, type="ln_q", power="1", comparison=">", rational="1/2")
    assert resp.get_json()["error"] == "请输入大于1的值"

    resp = post_calc(client, type="pi", power="1", comparison=">", rational="xyz")
    assert resp.get_json()["error"] == "请输入分子和分母"


def test_get_integral_image(client, monkeypatch):
    """Query params are forwarded to the family render_equation; latex returned."""
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
            **{k: PARAMS[k] for k in server.IMAGE_KEYS},
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
    assert seen["power"] == Fraction(1)
    assert seen["bound"] == Fraction(22, 7)


def test_get_integral_image_errors(client, monkeypatch):
    """Bad type -> 400; missing family module -> 502."""
    resp = client.get("/get_integral_image", query_string={"type": "bogus"})
    assert resp.status_code == 400
    assert "error" in resp.get_json()

    monkeypatch.setitem(server.render.FAMILY, "pi", "no_such_family")
    resp = client.get(
        "/get_integral_image",
        query_string={
            "type": "pi",
            "comparison": "<",
            "rational": "22/7",
            "coef": "1",
            **{k: PARAMS[k] for k in server.IMAGE_KEYS},
        },
    )
    assert resp.status_code == 502
    assert "error" in resp.get_json()


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
    """Empty problem -> 400; ValueError from the module -> error body."""
    resp = client.post("/decompose_inequality", data={"problem": "  "})
    assert resp.status_code == 400
    assert resp.get_json()["success"] is False

    def boom(problem):
        raise ValueError("无法解析该不等式")

    mod = types.ModuleType("attention_calculator.decompose")
    mod.decompose_inequality = boom
    install_decompose(monkeypatch, mod)
    resp = client.post("/decompose_inequality", data={"problem": "???"})
    assert resp.get_json() == {"success": False, "error": "无法解析该不等式"}
