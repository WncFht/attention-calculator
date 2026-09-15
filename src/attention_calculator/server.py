"""Flask app reproducing zhuyidao.net (注意力计算器).

Endpoints mirror the live site: POST /calculate, GET /get_integral_image,
POST /decompose_inequality, plus GET / and /en serving the single-page UI.
Error statuses/texts follow bench probes: format errors -> 400, domain and
search failures -> 404, internal render failures -> 500 (site's catch-all).
Run with ``python -m attention_calculator.server`` or
``waitress-serve --port=8080 attention_calculator.server:app``.
"""

from __future__ import annotations

import math
import re
from fractions import Fraction

from flask import Flask, jsonify, render_template, request

from . import engine, render, solve
from .kernels import TYPES

app = Flask(__name__)
app.json.ensure_ascii = False  # 错误文案为中文, 直接输出 UTF-8

# 搜索预算: e、pi 两类型指数上限 30, 其余 10 (见 docs/kernel-spec.md 搜索顺序)
EXPONENT_LIMIT = {"pi": 30, "e": 30}

# 站端只接受 "n" 或 "n/d" (非负整数组成); 小数、负数、其它写法都算格式错误
NUM_RE = re.compile(r"^\d+(/\d+)?$")

# 右侧有理数分子/分母必须 < 10^16 (probe: cap:*)
RATIONAL_CAP = 10**16

INTERNAL_ERROR = "服务器内部错误，请稍后再试"  # noqa: RUF001 -- 站端原文


def fail(message: str, status: int):
    """A failure JSON body; the site sends just ``{"error": ...}``."""
    return jsonify({"error": message}), status


def split_num(text: str) -> tuple[int, int] | None:
    """Split the site's 'n'/'n/d' wire format into (numerator, denominator).

    Returns None on format violation (empty, non-digit, negative, decimal);
    denominator 0 is returned as-is so the caller can pick the right message.
    """
    text = text.strip()
    if not NUM_RE.match(text):
        return None
    num, _, den = text.partition("/")
    return int(num), int(den) if den else 1


def domain_error(kind: str, power: Fraction) -> str | None:
    """Site's 404-level domain prechecks, with probed messages verbatim."""
    if kind in ("ln_q", "ln_q_square") and power <= 1:
        return "请在ln后输入一个大于1的数"
    if kind == "sin_q" and not 0 < float(power) < math.pi:
        return "请在sin后输入一个在(0,π)内的数"
    if kind in ("cos_q", "tan_q", "cot_q") and not 0 < float(power) < math.pi / 2:
        return f"请在{kind.split('_')[0]}后输入一个在(0,π/2)内的数"
    if kind in ("sin_pi_q", "cos_pi_q") and (
        power.denominator == 1 or not 0 < power < Fraction(1, 2)
    ):
        return "请在输入一个在(0,1/2)内的分数，本情况不支持整数"  # noqa: RUF001 -- 站端原文(含"在"字笔误)
    if kind == "artanh_q" and (power.denominator == 1 or not 0 < power < 1):
        return "请在输入一个在(0,1)内的分数，本情况不支持整数"  # noqa: RUF001 -- 站端原文(含"在"字笔误)
    if kind == "arcoth_q" and power <= 1:
        return "请在输入一个大于1的数"  # noqa: RUF001 -- 站端原文(含"在"字笔误)
    if kind in ("sin_q_degree", "cos_q_degree") and not 0 < power < 90:
        return "请在输入一个在(0,90)内的数"  # 站端原文(含"在"字笔误)
    return None


@app.get("/")
@app.get("/en")
@app.get("/en/")
def index():
    """Serve the calculator page; /en only flips the frontend's initial language."""
    return render_template("index.html")


@app.post("/calculate")
def calculate():
    """Run a proof search; wraps solve.prove into the site's response shape."""
    # type 缺省回退 pi（实测：不发 type 字段照常出解；发 type="" 报无效类型）
    kind = request.form.get("type", "pi")
    power = request.form.get("power", "").strip()
    comp = request.form.get("comparison", "")
    rational = request.form.get("rational", "").strip()

    if kind not in TYPES:
        return fail("无效的证明类型", 400)
    if comp not in (">", "<"):
        return fail("无效的不等号方向", 400)

    # 校验顺序（成对探针钉死）：右侧 format→分母→上限 全部先于左侧一切检查
    bound_parts = split_num(rational)
    if bound_parts is None:
        return fail("右侧有理数格式无效", 400)
    if bound_parts[1] == 0:
        return fail("右侧有理数分母不能为0", 400)
    if bound_parts[0] >= RATIONAL_CAP or bound_parts[1] >= RATIONAL_CAP:
        return fail("右侧有理数请输入小于10^16的整数或分数", 400)

    power_parts = split_num(power)
    if power_parts is None:
        return fail("左侧系数格式无效", 400)
    if power_parts[1] == 0:
        return fail("左侧系数分母不能为0", 400)
    if power_parts[0] >= RATIONAL_CAP or power_parts[1] >= RATIONAL_CAP:
        return fail("左侧系数请输入小于10^16的整数或分数", 400)
    power_val = Fraction(*power_parts)

    err = domain_error(kind, power_val)
    if err:
        return fail(err, 404)

    try:
        result = solve.prove(kind, power, comp, rational)
    except engine.WrongDirection:
        return fail("要证明的式子不等号方向反了", 404)
    except engine.NoSolution:
        limit = EXPONENT_LIMIT.get(kind, 10)
        return fail(f"在指数不超过{limit}的范围内未找到{comp}方向的解", 404)
    except ValueError as exc:
        return fail(str(exc), 400)
    except ModuleNotFoundError:
        return fail(INTERNAL_ERROR, 500)

    params = dict(result["parameters"])
    params.setdefault("unified_form", {})
    # 线上响应里 m/n 是 int, au/bu/cu/u_val 是字符串
    params["m"], params["n"] = int(params["m"]), int(params["n"])
    for k in ("au_val", "bu_val", "cu_val", "u_val"):
        params[k] = str(params[k])
    return jsonify(
        {
            "success": True,
            "type": kind,
            "parameters": params,
            "equations": {"solution": result["solution"]},
        }
    )


IMAGE_KEYS = ("m", "n", "a_val", "b_val", "c_val", "u_val", "au_val", "bu_val", "cu_val")


@app.get("/get_integral_image")
def get_integral_image():
    """Render the LaTeX proof equation for solved parameters via the kernel family.

    The site answers every render failure with a generic 500; keep that.
    """
    try:
        kind = request.args.get("type", "")
        comp = request.args.get("comparison", "")
        params = render.coerce_params({k: request.args.get(k, "") for k in IMAGE_KEYS})
        # coef/rational 原文交给渲染层回显（站端不约分）；先解析一遍，
        # 不可解析的输入和站端一样走 catch-all 500
        power = request.args.get("coef", "1").strip()
        bound = request.args.get("rational", "").strip()
        render.wire_fraction(power)
        render.wire_fraction(bound)
        equation = render.render_equation(params, kind, power, comp, bound)
    except Exception:
        return jsonify({"error": INTERNAL_ERROR}), 500
    return jsonify({"equation": equation})


@app.post("/decompose_inequality")
def decompose_inequality():
    """Decompose a composite inequality into basic-type sub-proofs."""
    problem = request.form.get("problem", "").strip()
    if not problem:
        return fail("请输入一个组合不等式。", 400)
    from . import decompose

    try:
        result = decompose.decompose_inequality(problem)
    except ValueError as exc:
        return fail(str(exc), 400)
    return jsonify(result)


def main():
    """Serve the app on 0.0.0.0:8080 via waitress."""
    from waitress import serve

    serve(app, listen="0.0.0.0:8080")


if __name__ == "__main__":
    main()
