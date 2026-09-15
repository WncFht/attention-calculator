"""Flask app reproducing zhuyidao.net (注意力计算器).

Endpoints mirror the live site: POST /calculate, GET /get_integral_image,
POST /decompose_inequality, plus GET / and /en serving the single-page UI.
Run with ``python -m attention_calculator.server`` or
``waitress-serve --port=8080 attention_calculator.server:app``.
"""

from __future__ import annotations

import re
from fractions import Fraction

from flask import Flask, jsonify, render_template, request

from . import engine, render, solve
from .kernels import TYPES

app = Flask(__name__)
app.json.ensure_ascii = False  # 错误文案为中文, 直接输出 UTF-8

# 搜索预算: e、pi 两类型指数上限 30, 其余 10 (见 docs/kernel-spec.md 搜索顺序)
EXPONENT_LIMIT = {"pi": 30, "e": 30}

# 缺系数时的提示, 与原页 JS 的校验文案一致
POWER_PROMPT = {
    "pi": "请输入π的系数",
    "e": "请输入e的系数",
    "pi_n": "请输入π的次数",
    "e_q": "请输入e的次数",
    "ln_q": "请输入ln后的值",
    "ln_q_square": "请输入ln后的值",
    "sin_q": "请输入sin后的值",
    "cos_q": "请输入cos后的值",
    "tan_q": "请输入tan后的值",
    "cot_q": "请输入cot后的值",
    "sin_q_degree": "请输入sin后的度数",
    "cos_q_degree": "请输入cos后的度数",
    "sin_pi_q": "请输入sin(qπ)内的值",
    "cos_pi_q": "请输入cos(qπ)内的值",
    "arctan_q": "请输入arctan后的值",
    "arccot_q": "请输入arccot后的值",
    "sinh_q": "请输入sinh后的值",
    "cosh_q": "请输入cosh后的值",
    "tanh_q": "请输入tanh后的值",
    "coth_q": "请输入coth后的值",
    "artanh_q": "请输入artanh后的值",
    "arcoth_q": "请输入arcoth后的值",
    "gamma": "请输入欧拉常数γ的系数",
    "golden": "请输入黄金分割率φ的系数",
    "catalan": "请输入卡塔兰常数C的系数",
    "zeta3": "请输入阿培里常数ζ(3)的系数",
    "e_pi": "请输入盖尔方德常数e^π的系数",
    "varpi": "请输入双纽线周率ϖ的系数",
    "gauss": "请输入高斯常数G的系数",
}

# 前端把有理数显示成 LaTeX \frac{n}{d} 再回传给 /get_integral_image
FRAC_RE = re.compile(r"^\\d?frac\{\s*(-?\d+)\s*\}\{\s*(-?\d+)\s*\}$")


def parse_bound(text: str) -> Fraction:
    """Parse 'n', 'n/d', or LaTeX '\\frac{n}{d}'/'\\dfrac{n}{d}' into a Fraction."""
    text = text.strip()
    m = FRAC_RE.match(text)
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    return Fraction(text)


def invalid(message: str):
    """A failed-input JSON body (HTTP 400), same shape as the site's errors."""
    return jsonify({"success": False, "error": message}), 400


@app.get("/")
@app.get("/en")
@app.get("/en/")
def index():
    """Serve the calculator page; /en only flips the frontend's initial language."""
    return render_template("index.html")


@app.post("/calculate")
def calculate():
    """Run a proof search; wraps solve.prove into the site's response shape."""
    kind = request.form.get("type", "")
    power = request.form.get("power", "").strip()
    comp = request.form.get("comparison", "")
    rational = request.form.get("rational", "").strip()

    if kind not in TYPES:
        return invalid("请选择要证明的不等式类型")
    prompt = POWER_PROMPT[kind]
    if not power:
        return invalid(prompt)
    try:
        power_val = parse_bound(power)
    except (ValueError, ZeroDivisionError):
        return invalid(prompt)
    if kind == "pi_n" and power_val.numerator > 10:
        return invalid("π的次数的分子不要超过10" if "/" in power else "π的次数不要超过10")
    if kind in ("ln_q", "ln_q_square") and power_val <= 1:
        return invalid("请输入大于1的值")
    if comp not in (">", "<"):
        return invalid("请选择 > 或者 <")
    if not rational:
        return invalid("请输入分子和分母")
    try:
        parse_bound(rational)
    except (ValueError, ZeroDivisionError):
        return invalid("请输入分子和分母")

    try:
        result = solve.prove(kind, power, comp, rational)
    except ModuleNotFoundError:
        return jsonify(
            {"success": False, "error": f"kernel family for type {kind!r} not available"}
        ), 502
    except engine.WrongDirection:
        return jsonify({"success": False, "error": "要证明的式子不等号方向反了"})
    except ValueError as exc:
        return invalid(str(exc))
    except engine.NoSolution:
        limit = EXPONENT_LIMIT.get(kind, 10)
        return jsonify(
            {"success": False, "error": f"在指数不超过{limit}的范围内未找到{comp}方向的解"}
        )
    return jsonify(
        {
            "success": True,
            "type": kind,
            "parameters": result["parameters"],
            "equations": {"solution": result["solution"]},
        }
    )


IMAGE_KEYS = ("m", "n", "a_val", "b_val", "c_val", "u_val", "au_val", "bu_val", "cu_val")


@app.get("/get_integral_image")
def get_integral_image():
    """Render the LaTeX proof equation for solved parameters via the kernel family."""
    kind = request.args.get("type", "")
    if kind not in TYPES:
        return jsonify({"error": "unsupported type"}), 400
    comp = request.args.get("comparison", "")
    if comp not in (">", "<"):
        return jsonify({"error": "bad comparison"}), 400
    try:
        params = render.coerce_params({k: request.args.get(k, "") for k in IMAGE_KEYS})
        power = parse_bound(request.args.get("coef", "1"))
        bound = parse_bound(request.args.get("rational", ""))
    except (ValueError, ZeroDivisionError):
        return jsonify({"error": "bad parameters"}), 400
    try:
        equation = render.render_equation(params, kind, power, comp, bound)
    except ModuleNotFoundError:
        return jsonify({"error": f"kernel family for type {kind!r} not available"}), 502
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    return jsonify({"equation": equation})


@app.post("/decompose_inequality")
def decompose_inequality():
    """Decompose a composite inequality into basic-type sub-proofs."""
    problem = request.form.get("problem", "").strip()
    if not problem:
        return invalid("请输入一个组合不等式。")
    from . import decompose
    try:
        result = decompose.decompose_inequality(problem)
    except ValueError as exc:
        return jsonify({"success": False, "error": str(exc)}), 400
    return jsonify(result)


def main():
    """Serve the app on 0.0.0.0:8080 via waitress."""
    from waitress import serve

    serve(app, listen="0.0.0.0:8080")


if __name__ == "__main__":
    main()
