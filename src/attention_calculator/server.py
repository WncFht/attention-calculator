"""Flask app reproducing zhuyidao.net (注意力计算器).

Endpoints mirror the live site: POST /calculate, GET /get_integral_image,
POST /decompose_inequality, plus GET / and /en serving the single-page UI.
Error statuses/texts follow bench probes: format errors -> 400, domain and
search failures -> 404, internal render failures -> 500 (site's catch-all).
Run with ``python -m attention_calculator.server`` or
``waitress-serve --port=8080 attention_calculator.server:app``.
"""

from __future__ import annotations

import json
import math
import re
from fractions import Fraction

from flask import Flask, Response, render_template, request

from . import engine, render, solve
from .kernels import TYPES

app = Flask(__name__)
# Jinja 默认吞掉模板尾部换行；站端页面以 \n 结尾（byte diff 实测）
app.jinja_env.keep_trailing_newline = True


def respond(payload: dict, status: int = 200) -> Response:
    """JSON body matching the live site byte-for-byte: compact separators,
    \\uXXXX escapes (ensure_ascii), alphabetically sorted keys, and the
    trailing newline old Flask's jsonify appended."""
    body = json.dumps(payload, ensure_ascii=True, sort_keys=True,
                      separators=(",", ":")) + "\n"
    return Response(body, status=status, mimetype="application/json")

# 搜索预算: e、pi 两类型指数上限 30, 其余 10 (见 docs/kernel-spec.md 搜索顺序)
EXPONENT_LIMIT = {"pi": 30, "e": 30}

# 站端只接受 "n" 或 "n/d" (非负整数组成); 小数、负数、其它写法都算格式错误
NUM_RE = re.compile(r"^\d+(/\d+)?$")

# 右侧有理数分子/分母必须 < 10^16 (probe: cap:*)
RATIONAL_CAP = 10**16

INTERNAL_ERROR = "服务器内部错误，请稍后再试"  # noqa: RUF001 -- 站端原文

# 站端把一切非 2xx 包成 JSON：未知路径 404 文案固定，方法不对走 catch-all 500
NOT_FOUND = "请求的页面不存在"


def fail(message: str, status: int):
    """A failure JSON body; the site sends just ``{"error": ...}``."""
    return respond({"error": message}, status)


@app.errorhandler(404)
def not_found(_):
    """Site answers unknown paths with a JSON body, not Flask's HTML page."""
    return fail(NOT_FOUND, 404)


@app.errorhandler(405)
def method_not_allowed(_):
    """Wrong method -> the site's generic 500 (e.g. GET /calculate)."""
    return fail(INTERNAL_ERROR, 500)


@app.errorhandler(Exception)
def unhandled(_):
    """Everything else non-2xx is also JSON-wrapped server-side."""
    return fail(INTERNAL_ERROR, 500)


def split_num(text: str) -> tuple[int, int] | None:
    """Split the site's 'n'/'n/d' wire format into (numerator, denominator).

    站端把空格和换行当透明字符全局删掉再校验（"2 2/7"→22/7、"1\n"→1 均实测），
    tab 等其它空白不赦免（"\\t22/7"→400）。Returns None on format violation;
    denominator 0 is returned as-is so the caller can pick the right message.
    """
    text = text.replace(" ", "").replace("\n", "")
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
def index():
    """Serve the calculator page; /en is byte-identical (language is client-side)."""
    return render_template("index.html")


@app.get("/en/")
def en_slash():
    """The site 404s on /en/ even though /en works (no slash redirect)."""
    return fail(NOT_FOUND, 404)


# 同一应用还挂在 /attention 前缀下（页面相同，静态资源路径前缀改写）
@app.get("/attention")
@app.get("/attention/")
@app.get("/attention/en")
def attention():
    """Serve the calculator under the /attention mount the live site exposes."""
    return render_template("attention.html")


@app.get("/attention/static/<path:name>")
def attention_static(name: str):
    """Mirror of /static under the /attention prefix."""
    from flask import send_from_directory

    return send_from_directory(app.static_folder, name)


@app.get("/favicon.ico")
def favicon():
    """The site answers the favicon with an empty 204, not a file."""
    return "", 204


@app.post("/calculate")
def calculate():
    """Run a proof search; wraps solve.prove into the site's response shape."""
    # type 缺省回退 pi（实测：不发 type 字段照常出解；发 type="" 报无效类型）
    kind = request.form.get("type", "pi")
    # power/rational 不做 strip——空白字符规则由 split_num 统一实现
    # （站端删掉 " "/"\\n" 但不赦免 tab，strip() 会误吃两端 tab）
    power = request.form.get("power", "")
    # 缺席 -> 站端默认 '>'；存在但为空 -> 400（缺席与空串不同待遇）
    comp = request.form.get("comparison", ">")
    rational = request.form.get("rational", "")

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

    # 求解器拿到的是清理后的 n/d 文本（站端校验即清理；原始串可能含空格）
    power_wire = "/".join(map(str, power_parts))
    rational_wire = "/".join(map(str, bound_parts))
    try:
        result = solve.prove(kind, power_wire, comp, rational_wire)
    except engine.WrongDirection:
        return fail("要证明的式子不等号方向反了", 404)
    except engine.NoSolution:
        limit = EXPONENT_LIMIT.get(kind, 10)
        return fail(f"在指数不超过{limit}的范围内未找到{comp}方向的解", 404)
    except engine.EqualClaim as exc:
        return fail(str(exc), 404)
    except ValueError as exc:
        return fail(str(exc), 400)
    except Exception:
        # 站端对一切内核异常走 catch-all 500（probe: e_q 0、ln_q_square 5/7）
        return fail(INTERNAL_ERROR, 500)

    params = dict(result["parameters"])
    params.setdefault("unified_form", {})
    # 线上响应里 m/n 是 int, au/bu/cu/u_val 是字符串
    params["m"], params["n"] = int(params["m"]), int(params["n"])
    for k in ("au_val", "bu_val", "cu_val", "u_val"):
        params[k] = str(params[k])
    return respond(
        {
            "success": True,
            # kernels may normalize the wire type (site echoes sin_q_degree
            # requests back as "sin_pi_q"); default to the request's type
            "type": result.get("type", kind),
            "parameters": params,
            "equations": {"solution": result["solution"]},
        }
    )


# 图像端点逐字段校验（probe 钉死）：
# 顺序 type -> comparison -> 整数字段 m,n,u,au,bu,cu -> 分数字段 a,b,c；
# coef/rational 完全不校验、原文回显（coef=x -> x\pi，rational 缺席 -> 0）。
IMG_FRAC_RE = re.compile(r"^-?\d+(/\d+)?$")
# 整数字段语法 = strip() 后 ^-?\d+$：两端空白（含 tab）可、"+3"/内部空白/下划线拒
# （probe: m=+3 → m必须是整数——站端不是裸 int()，int() 会吃掉 '+'）
IMG_INT_RE = re.compile(r"^-?\d+$")
IMG_INT_KEYS = ("m", "n", "u_val", "au_val", "bu_val", "cu_val")
IMG_FRAC_KEYS = ("a_val", "b_val", "c_val")
IMG_INT_BOUNDS = {"m": (0, 30), "n": (0, 30)}
# u_val 下限按类型：gamma 两段式核允许 u=0（golden 实测），其余要 >=1
IMG_U_MIN = {"gamma": 0}


@app.get("/get_integral_image")
def get_integral_image():
    """Render the LaTeX proof equation for solved parameters via the kernel family."""
    kind = request.args.get("type", "pi")
    comp = request.args.get("comparison", ">")
    if not kind or kind not in TYPES:
        return fail("无效的证明类型", 400)
    if comp not in (">", "<"):
        return fail("无效的不等号方向", 400)

    params = {}
    for k in IMG_INT_KEYS:
        v = request.args.get(k, "")
        # 先 strip() 再过 ^-?\d+$（见 IMG_INT_RE 注释）
        if not IMG_INT_RE.match(v.strip()):
            return fail(f"{k}必须是整数", 400)
        params[k] = int(v)
        lo, hi = IMG_INT_BOUNDS.get(k, (None, None))
        if k == "u_val":
            lo = IMG_U_MIN.get(kind, 1)
        if lo is not None and params[k] < lo:
            return fail(f"{k}过小", 400)
        if hi is not None and params[k] > hi:
            return fail(f"{k}过大", 400)
    for k in IMG_FRAC_KEYS:
        # 分数字段同 /calculate 语法：先去 " "/"\\n" 再卡格式
        v = request.args.get(k, "").replace(" ", "").replace("\n", "")
        if not IMG_FRAC_RE.match(v):
            return fail(f"{k}格式无效", 400)
        _, _, den = v.partition("/")
        if den and int(den) == 0:
            return fail(f"{k}分母不能为0", 400)
        params[k] = Fraction(v)

    power = request.args.get("coef", "1")
    bound = request.args.get("rational", "0")
    try:
        equation = render.render_equation(params, kind, power, comp, bound)
    except Exception:
        return fail(INTERNAL_ERROR, 500)
    return respond({"equation": equation})


@app.post("/decompose_inequality")
def decompose_inequality():
    """Decompose a composite inequality into basic-type sub-proofs."""
    problem = request.form.get("problem", "").strip()
    if not problem:
        return fail("请输入一个只包含一个 > 或 < 的不等式", 400)
    from . import decompose

    try:
        result = decompose.decompose_inequality(problem)
    except ValueError as exc:
        return fail(str(exc), 400)
    return respond(result)


def main():
    """Serve the app on 0.0.0.0:8080 via waitress."""
    from waitress import serve

    serve(app, listen="0.0.0.0:8080")


if __name__ == "__main__":
    main()
