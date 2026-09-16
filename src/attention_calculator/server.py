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
from .kernels import EXACT_TYPES, TYPES

app = Flask(__name__)
# Jinja 默认吞掉模板尾部换行；站端页面以 \n 结尾（byte diff 实测）
app.jinja_env.keep_trailing_newline = True


def respond(payload: dict, status: int = 200) -> Response:
    """JSON body matching the live site byte-for-byte: compact separators,
    \\uXXXX escapes (ensure_ascii), alphabetically sorted keys, and the
    trailing newline old Flask's jsonify appended."""
    body = json.dumps(payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")) + "\n"
    return Response(body, status=status, mimetype="application/json")


# 站端只接受 "n" 或 "n/d" (非负整数组成); 小数、负数、其它写法都算格式错误
NUM_RE = re.compile(r"^\d+(/\d+)?$")

# 右侧有理数分子/分母必须 < 10^16 (probe: cap:*)
RATIONAL_CAP = 10**16

INTERNAL_ERROR = "服务器内部错误，请稍后再试"  # 站端原文
# 姊妹应用（/convex、/health）的内部错误文案带句号，且所有未匹配/错方法
# 一律 500（probe: /health/xyz、/health/en/、POST /health/en、/convex/Prove）
SIBLING_INTERNAL_ERROR = "服务器内部错误，请稍后再试。"  # 站端原文（带句号）

# 站端把一切非 2xx 包成 JSON：未知路径 404 文案固定，方法不对走 catch-all 500
NOT_FOUND = "请求的页面不存在"


def in_sibling(path: str) -> bool:
    """路径是否落在姊妹应用挂载点内；'/healthxyz' 这类前缀兄弟不算
    （probe: GET /healthxyz -> 主站 404 裸包络）。"""
    return any(path == p or path.startswith(p + "/") for p in ("/convex", "/health"))


def fail(message: str, status: int):
    """A failure JSON body; the site sends just ``{"error": ...}``.

    /convex 与 /health 下的姊妹应用用新版包络 ``{"error": ..., "ok": false}``
    （见 docs/sibling-apps.md §0）。
    """
    payload = {"error": message}
    if in_sibling(request.path):
        payload["ok"] = False
    return respond(payload, status)


@app.errorhandler(404)
def not_found(_):
    """Site answers unknown paths with a JSON body, not Flask's HTML page.
    姊妹应用域内未匹配路径不走 404，一律 500 内部错误（probe: /health/xyz）。"""
    if in_sibling(request.path):
        return fail(SIBLING_INTERNAL_ERROR, 500)
    return fail(NOT_FOUND, 404)


@app.errorhandler(405)
def method_not_allowed(_):
    """Wrong method -> the site's generic 500 (e.g. GET /calculate)."""
    if in_sibling(request.path):
        return fail(SIBLING_INTERNAL_ERROR, 500)
    return fail(INTERNAL_ERROR, 500)


@app.errorhandler(Exception)
def unhandled(_):
    """Everything else non-2xx is also JSON-wrapped server-side."""
    if in_sibling(request.path):
        return fail(SIBLING_INTERNAL_ERROR, 500)
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
    if kind in ("ln_q", "ln_q_square", "ln_q_cube") and power <= 1:
        return "请在ln后输入一个大于1的数"
    if kind == "arcsin_q" and not 0 < power < 1:
        return "arcsin后的值只能在(0,1)内，请输入一个在(0,1)内的分数"
    if kind == "arsinh_q" and power == 0:
        return "请在arsinh后输入一个非0的数"
    if kind == "sin_q" and not 0 < float(power) < math.pi:
        return "请在sin后输入一个在(0,π)内的数"
    if kind in ("cos_q", "tan_q", "cot_q") and not 0 < float(power) < math.pi / 2:
        return f"请在{kind.split('_')[0]}后输入一个在(0,π/2)内的数"
    if kind in ("sin_pi_q", "cos_pi_q") and (
        power.denominator == 1 or not 0 < power < Fraction(1, 2)
    ):
        return "请在输入一个在(0,1/2)内的分数，本情况不支持整数"  # 站端原文(含"在"字笔误)
    if kind == "artanh_q" and (power.denominator == 1 or not 0 < power < 1):
        return "请在输入一个在(0,1)内的分数，本情况不支持整数"  # 站端原文(含"在"字笔误)
    if kind == "arcoth_q" and power <= 1:
        return "请在输入一个大于1的数"  # 站端原文(含"在"字笔误)
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


# 本地演示页：站端页面引用 jsdelivr/hertzen 三个 CDN 依赖（MathJax/KaTeX/
# html2canvas），离线或 CDN 不可达时公式不排版。/demo 把 URL 改写为
# static/vendor 下的本地副本；/ 与 /en 保持站端逐字节不动（parity 夹具）。
DEMO_CDN = (
    ("https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js", "mathjax/tex-mml-chtml.js"),
    ("https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css", "katex/katex.min.css"),
    ("https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js", "katex/katex.min.js"),
    ("https://html2canvas.hertzen.com/dist/html2canvas.min.js", "html2canvas.min.js"),
)
DEMO_ASSETS = {cdn: f"/static/vendor/{name}" for cdn, name in DEMO_CDN}


@app.get("/demo")
def demo():
    """Same page as / with the three CDN URLs rewritten to vendored copies."""
    html = render_template("index.html")
    for cdn, local in DEMO_ASSETS.items():
        html = html.replace(cdn, local)
    return html


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

    # mode=exact 走数学正确性路径（docs/2026-09-16-math-correctness-plan.md）；
    # 缺省/其他值保持站端逐字节行为。先读 mode：EXACT_TYPES 只在 exact 下有效
    exact = request.form.get("mode") == "exact"
    if kind not in TYPES and not (exact and kind in EXACT_TYPES):
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
        result = solve.prove(kind, power_wire, comp, rational_wire, exact=exact)
    except (engine.WrongDirection, engine.NoSolution) as exc:
        return fail(solve.failure_text(exc, kind, comp), 404)
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
    payload = {
        "success": True,
        # kernels may normalize the wire type (site echoes sin_q_degree
        # requests back as "sin_pi_q"); default to the request's type
        "type": result.get("type", kind),
        "parameters": params,
        "equations": {"solution": result["solution"]},
    }
    # mode=exact 的证明附机器可检证书（W5）；site 模式无此键，逐字节不变
    if "certificate" in result:
        payload["certificate"] = result["certificate"]
    return respond(payload)


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
    except Exception:
        return fail("组合证明生成失败，请稍后再试", 500)
    return respond(result)


# ===== 姊妹应用 /convex（凹凸不等式计算器，docs/sibling-apps.md §1） =====


@app.get("/convex")
@app.get("/convex/")
@app.get("/convex/en")
def convex_page():
    """站端 /convex、/convex/ 与 /convex/en 返回同一份字节（语言由前端 JS 切换；
    无斜杠路径直接 200 不做跳转，probe page2:convex 实测）。"""
    return render_template("convex.html")


@app.get("/convex/static/<path:name>")
def convex_static(name: str):
    """/convex/static/* 与主站共用 static 目录（bg1.png sha256 实测一致）。"""
    from flask import send_from_directory

    return send_from_directory(app.static_folder, name)


@app.post("/convex/prove")
def convex_prove():
    """凹凸不等式证明端点；convex 模块惰性导入（仅在命中路由时加载）。"""
    from . import convex

    # 前置校验按 strip 后文本：空 -> 400；strip 后 >500 -> 输入过长（恰好 500 通过）。
    # line/domain 无长度限制（probe: line 3001 字符照常进解析）
    inequality = request.form.get("inequality", "").strip()
    if not inequality:
        return fail("请输入一个不等式。", 400)
    if len(inequality) > 500:
        return fail("输入过长，请输入一个较短的单变量不等式", 400)
    line = request.form.get("line") or None  # 隐藏参数，缺席与空串同等处理
    domain = request.form.get("domain") or None  # '0,inf'/'0,10' 形；缺席走默认域
    try:
        result = convex.prove(inequality, line, domain)
    except Exception as exc:
        # 站端 prove 包装器把一切内核异常（不限 ValueError）包成 400 str(exc)：
        # "x/0>0"→"float division by zero"、x^1e309→"cannot convert Infinity
        # to integer ratio"、line 超长→"maximum recursion depth exceeded"
        return fail(str(exc), 400)
    return respond({"ok": True, "result": result})


# ===== 姊妹应用 /health（健康计算器，docs/sibling-apps.md §2） =====


@app.get("/health")
@app.get("/health/")
def health_page():
    """/health 与 /health/ 均 200 返回同一份中文页面（两路并存，无斜杠跳转）。"""
    return render_template("health.html")


@app.get("/health/en")
def health_page_en():
    """英文版是独立模板；/health/en/ 无路由（probe: 500 内部错误）。"""
    return render_template("health-en.html")


@app.post("/health/calculate")
def health_calculate():
    """JSON-only 计算端点；校验错误 400、成功 200，新版 ``ok`` 包络。"""
    from . import health

    # 非 JSON Content-Type、坏 JSON、非 dict 顶层一律同一条错（probe: tr:*）
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        return fail("提交内容格式不正确。", 400)
    fields, err = health.validate(data)
    if err:
        return fail(err, 400)
    return respond(
        {"ok": True, "record_id": health.next_record_id(), "results": health.calculate(fields)}
    )


def main():
    """Serve the app on 0.0.0.0:8080 via waitress."""
    from waitress import serve

    serve(app, listen="0.0.0.0:8080")


if __name__ == "__main__":
    main()
