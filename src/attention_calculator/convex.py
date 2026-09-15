"""POST /convex/prove 求解器 —— zhuyidao.net「凹凸不等式计算器」克隆。

管线（契约见 docs/sibling-apps.md §1，全部行为按该文档实测记录实现）:
ast 解析受限文法 → 项的线性组合 → 归一 left > right → 二阶导符号分类 →
定义域数值最小值 → (proved|inconclusive|failed) → 有理切点切线搜索。

错误面与站端逐字一致: ast.dump(node, show_empty=True) 直接进错误串
（`unsupported atom Call(func=Name(id='foo', ctx=Load()), ...)`）、
`log only supports argument x`、缺字段 `请输入一个不等式。`。

待探测钉死的细节全部收在 ProbeConfig 与候选生成器接缝中——数值对齐时
只改配置默认值，不动结构:

* 定义域: 含 log/sqrt/分数幂时 x>1e-8（实测"从 1e-08 起扫"）；无受限
  原子时是否退回全体实数未钉死（domain_pos_only / domain_neg_lo）。
* 扫描上界 domain_hi、网格数 grid_n、网格对数排布——均为假设值。
* 数值最小值机制: 网格取优 + f' 变号单元 brentq 求根（站端 15 位
  x 指向导数求根；bounded Brent 只到 ~1e-8）。bracket 选区与
  root_xtol 影响最后几个 ulp，样本最小值差 1 ulp。
* 切点候选集: 默认假设 = 数值最小值点 xmin 的连分数渐近分数列（样本
  tangent_at=17/30 正是 xmin≈0.567143290409784 的第 5 个渐近分数）。
  真实候选集与排除规则待探测。
* 阈值: 数值最小值判负 min_tol、切线 gap 判正 gap_tol（实测 e^x>=x+1
  的 ~0 gap 被拒 → 存在正阈值或边界判定）。
* 未实测的错误文案: 非比较输入 / 非法比较符 / 非仿射 line / 空侧渲染
  / 首项负号形态等，见各 PROBE 注释。
"""

from __future__ import annotations

import ast
import math
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise

import sympy as sp
from scipy.optimize import brentq

X = sp.symbols("x")

# ---- 站端原文文案（实测） ----
MISSING_INEQUALITY = "请输入一个不等式。"
REASON_INCONCLUSIVE = "整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。"  # noqa: E501 -- 站端原文
REASON_FAILED = "数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。"  # 站端原文
REASON_NO_LINE = "不等式数值上已通过，但当前情形没有生成中间直线证明。"  # 站端原文
REASON_SEARCH_MISS = "不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。"  # 站端原文

# ---- 未实测文案（PROBE） ----
NOT_INEQUALITY = MISSING_INEQUALITY  # PROBE: 非比较/语法错误输入的文案未实测
LINE_NOT_AFFINE = "line 只支持 m*x+b 形式的直线。"  # PROBE: 文案未实测

# Term = (coef: Fraction, atom)；atom = ("exp"|"log"|"sqrt"|"x"|"const")
#        或 ("pow", Fraction)；求值原子 evalatom = ("logat"|"sqrtat"|"expat", q)
#        或 ("powat", q, r)，只用于切线截距的字面显示。
Term = tuple[Fraction, tuple]
Atom = tuple
EvalAtom = tuple


# ============================ 解析 ============================

def is_const_expr(node: ast.expr) -> bool:
    """True if the AST subtree contains no x and no calls (pure arithmetic)."""
    if isinstance(node, ast.Constant):
        return isinstance(node.value, (int, float))
    if isinstance(node, ast.UnaryOp):
        return isinstance(node.op, (ast.UAdd, ast.USub)) and is_const_expr(node.operand)
    if isinstance(node, ast.BinOp):
        return (isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow))
                and is_const_expr(node.left) and is_const_expr(node.right))
    return False


def eval_const(node: ast.expr) -> Fraction:
    """Evaluate an is_const_expr subtree to a Fraction.

    小数按 str() 十进制串收纳（0.1 -> 1/10）；PROBE: 站端也可能用
    Fraction(0.1) 的二进制有理形态。
    """
    if isinstance(node, ast.Constant):
        v = node.value
        return Fraction(v) if isinstance(v, int) else Fraction(str(v))
    if isinstance(node, ast.UnaryOp):
        val = eval_const(node.operand)
        return -val if isinstance(node.op, ast.USub) else val
    a, b = eval_const(node.left), eval_const(node.right)
    op = node.op
    if isinstance(op, ast.Add):
        return a + b
    if isinstance(op, ast.Sub):
        return a - b
    if isinstance(op, ast.Mult):
        return a * b
    if isinstance(op, ast.Div):
        return a / b
    # Pow：整数指数精确，分数指数经 float 再按十进制串回收
    res = a ** b
    return res if isinstance(res, Fraction) else Fraction(str(res))


def parse_atom(node: ast.expr) -> Atom:
    """Term atom: exp/log/sqrt call on bare x, x^rational, or x itself."""
    if isinstance(node, ast.Name):
        if node.id == "x":
            return ("x",)
        raise ValueError(f"unsupported atom {ast.dump(node, show_empty=True)}")
    if isinstance(node, ast.Call):
        fname = node.func.id if isinstance(node.func, ast.Name) else ""
        if fname not in ("exp", "log", "sqrt"):
            raise ValueError(f"unsupported atom {ast.dump(node, show_empty=True)}")
        if (len(node.args) != 1 or node.keywords
                or not isinstance(node.args[0], ast.Name) or node.args[0].id != "x"):
            raise ValueError(f"{fname} only supports argument x")
        return (fname,)
    if (isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow)
            and isinstance(node.left, ast.Name) and node.left.id == "x"
            and is_const_expr(node.right)):
        return ("pow", eval_const(node.right))
    raise ValueError(f"unsupported atom {ast.dump(node, show_empty=True)}")


def mul_factors(node: ast.expr) -> list[ast.expr]:
    """Flatten a Mult chain into its factor list."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Mult):
        return mul_factors(node.left) + mul_factors(node.right)
    return [node]


def parse_term(node: ast.expr) -> Term:
    """One additive term -> (coefficient, atom).

    `系数*原子` / `原子/常数` 剥出系数；纯常数项归 ("const",)。多个非常数
    因子（x*x）或非常数分母 -> unsupported atom（项整体 dump）。
    """
    if is_const_expr(node):
        return eval_const(node), ("const",)
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
        if isinstance(node.op, ast.Div):
            numer, denom = [node.left], [node.right]
        else:
            numer, denom = mul_factors(node), []
        coef = Fraction(1)
        atoms = []
        for f in numer:
            if is_const_expr(f):
                coef *= eval_const(f)
            else:
                atoms.append(parse_atom(f))  # 各因子自己的错误串（Name('y') 等）
        for f in denom:
            if not is_const_expr(f):
                raise ValueError(f"unsupported atom {ast.dump(node, show_empty=True)}")
            coef /= eval_const(f)
        if len(atoms) != 1:
            raise ValueError(f"unsupported atom {ast.dump(node, show_empty=True)}")
        return coef, atoms[0]
    return Fraction(1), parse_atom(node)


def parse_sum(node: ast.expr) -> list[Term]:
    """An expression -> signed (coef, atom) term list, input order preserved."""
    out: list[Term] = []

    def walk(n: ast.expr, neg: bool) -> None:
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            walk(n.left, neg)
            walk(n.right, neg)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            walk(n.left, neg)
            walk(n.right, not neg)
        elif isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            walk(n.operand, neg != isinstance(n.op, ast.USub))
        else:
            c, a = parse_term(n)
            out.append((-c if neg else c, a))

    walk(node, False)
    return out


def parse_expression(text: str) -> ast.expr:
    """共享入口: `^`→`**` 后 ast.parse（实测二者响应字节相同）。"""
    try:
        return ast.parse(text.replace("^", "**"), mode="eval").body
    except SyntaxError:
        raise ValueError(NOT_INEQUALITY) from None


def parse_inequality(text: str) -> list[Term]:
    """`lhs > rhs`（含 >=/< /<=）归一为 difference 项列表 lhs - rhs。"""
    node = parse_expression(text)
    if not isinstance(node, ast.Compare) or len(node.ops) != 1:
        raise ValueError(NOT_INEQUALITY)
    op = node.ops[0]
    if isinstance(op, (ast.Gt, ast.GtE)):
        lhs, rhs = node.left, node.comparators[0]
    elif isinstance(op, (ast.Lt, ast.LtE)):
        lhs, rhs = node.comparators[0], node.left
    else:
        raise ValueError(NOT_INEQUALITY)  # PROBE: ==/!= 等比较符文案未实测
    terms = parse_sum(lhs) + [(-c, a) for c, a in parse_sum(rhs)]
    return [(c, a) for c, a in terms if c != 0]


# ============================ 显示 ============================

def fmt_frac(q: Fraction, tex: bool) -> str:
    """`261/112` / `\\frac{261}{112}`；调用方传绝对值，符号由 join 处理。"""
    if q.denominator == 1:
        return str(q.numerator)
    if tex:
        return f"\\frac{{{q.numerator}}}{{{q.denominator}}}"
    return f"{q.numerator}/{q.denominator}"


def atom_render(atom: Atom, tex: bool) -> str:
    """`e^x` `ln x` `√x` `x^1.5` / latex `\\ln x` `\\sqrt{x}` `x^{\\frac{3}{2}}`。"""
    kind = atom[0]
    if kind == "exp":
        return "e^x"
    if kind == "log":
        return "\\ln x" if tex else "ln x"
    if kind == "sqrt":
        return "\\sqrt{x}" if tex else "√x"
    if kind == "x":
        return "x"
    r = atom[1]  # pow: 文本态指数按小数（x^(3/2) -> x^1.5），latex 态按分数
    if r.denominator == 1:
        e = str(r.numerator)
    else:
        e = f"\\frac{{{r.numerator}}}{{{r.denominator}}}" if tex else str(float(r))
    return f"x^{{{e}}}" if tex else f"x^{e}"


def term_render(term: Term, tex: bool) -> str:
    """|coef|·atom；系数 1 裸原子。text 用 `*` 连接，latex 直接相邻。"""
    c, atom = abs(term[0]), term[1]
    if atom[0] == "const":
        return fmt_frac(c, tex)
    body = atom_render(atom, tex)
    if c == 1:
        return body
    return f"{fmt_frac(c, tex)}{body}" if tex else f"{fmt_frac(c, tex)}*{body}"


def render_terms(terms: list[Term], tex: bool) -> str:
    """` + `/` - ` 连接的项列表；空侧显示 0（PROBE: 未实测）；首项负号 `- `。"""
    if not terms:
        return "0"
    out = []
    for i, term in enumerate(terms):
        body = term_render(term, tex)
        if i == 0:
            out.append(f"- {body}" if term[0] < 0 else body)
        else:
            out.append((" - " if term[0] < 0 else " + ") + body)
    return "".join(out)


def evalatom_render(ea: EvalAtom, tex: bool) -> str:
    """切线截距里的"原子在 x0 取值"字面显示: `ln(17/30)` / `\\ln \\frac{17}{30}`。"""
    kind, q = ea[0], ea[1]
    if kind == "logat":
        return f"\\ln {fmt_frac(q, True)}" if tex else f"ln({q})"
    if kind == "sqrtat":
        return f"\\sqrt{{{fmt_frac(q, True)}}}" if tex else f"√({q})"
    if kind == "expat":
        return f"e^{{{fmt_frac(q, True)}}}" if tex else f"e^({q})"
    r = ea[2]  # powat
    if r.denominator == 1:
        e = str(r.numerator)
    else:
        e = f"\\frac{{{r.numerator}}}{{{r.denominator}}}" if tex else str(float(r))
    return f"\\left({fmt_frac(q, True)}\\right)^{{{e}}}" if tex else f"({q})^{e}"


def evalterm_render(coef: Fraction, ea: EvalAtom, tex: bool) -> str:
    """|coef|·evalatom（截距显示项）。"""
    body = evalatom_render(ea, tex)
    if abs(coef) == 1:
        return body
    c = fmt_frac(abs(coef), tex)
    return f"{c}{body}" if tex else f"{c}*{body}"


def line_render(m: Fraction, num: Fraction, items: list[tuple[Fraction, EvalAtom]],
                const: Fraction, tex: bool) -> str:
    """分离直线: `(m*x + 截距) + 右侧常数`，截距 = -m·x0 + Σ c_i·atom_i(x0)。

    实测怪癖: 斜率 1 印 `1x`/`1*x`（不约简）、`\\ln 1` 不化简、截距常数项
    排在求值项之前（`- 1 + ln(17/30)`）。
    """
    mx = f"{fmt_frac(m, True)}x" if tex else f"{fmt_frac(m, False)}*x"
    seq: list[tuple[Fraction, str]] = []
    if num or not items:
        seq.append((num, fmt_frac(abs(num), tex)))
    seq += [(c, evalterm_render(c, ea, tex)) for c, ea in items]
    inner = mx + "".join((" - " if c < 0 else " + ") + body for c, body in seq)
    line = f"({inner})"
    if const:
        line += (" - " if const < 0 else " + ") + fmt_frac(abs(const), tex)
    return line


# ============================ 归一与分类 ============================

def atom_curv(atom: Atom) -> int:
    """正系数原子的二阶导符号: +1 凸 / -1 凹 / 0 仿射。"""
    kind = atom[0]
    if kind == "exp":
        return 1
    if kind in ("log", "sqrt"):
        return -1
    if kind == "pow":
        r = atom[1]
        return 1 if (r > 1 or r < 0) else (-1 if 0 < r < 1 else 0)
    return 0  # x / const


def normalize(terms: list[Term]) -> tuple[list[Term], list[Term]]:
    """difference 拆成 left > right（实测规则，sibling-apps.md §1）:

    正系数项留左侧（凹凸留给分类器裁决）；负系数项取负后若为凹/仿射则
    移到右侧（右侧项系数恒正）；常数按符号归属（正留左、负移右）。
    """
    left, right = [], []
    for c, a in terms:
        if a[0] == "const":
            (left if c > 0 else right).append((abs(c), a))
        elif c > 0 or atom_curv(a) > 0:
            left.append((c, a))
        else:
            right.append((-c, a))
    return left, right


def atom_expr(atom: Atom) -> sp.Expr:
    """atom -> sympy 表达式。"""
    kind = atom[0]
    if kind == "exp":
        return sp.exp(X)
    if kind == "log":
        return sp.log(X)
    if kind == "sqrt":
        return sp.sqrt(X)
    if kind == "pow":
        return X ** sp.Rational(atom[1].numerator, atom[1].denominator)
    if kind == "x":
        return X
    return sp.Integer(1)  # const


def sum_expr(terms: list[Term]) -> sp.Expr:
    """项列表 -> sympy 表达式（系数走 Rational 保持精确）。"""
    total = sp.Integer(0)
    for c, a in terms:
        total += sp.Rational(c.numerator, c.denominator) * atom_expr(a)
    return total


def safe_func(f: Callable[[float], float]) -> Callable[[float], float]:
    """float 函数包装: 域外/奇异点返回 inf（真实数值边界）。"""

    def g(v: float) -> float:
        try:
            out = f(v)
        except (ArithmeticError, ValueError, TypeError):
            return math.inf
        return out if math.isfinite(out) else math.inf

    return g


def domain_of(terms: list[Term], cfg: ProbeConfig) -> tuple[float, float]:
    """定义域: 含 log/sqrt/分数幂/负幂 -> (1e-8, hi)（实测下界 1e-08）。"""
    need_pos = any(
        a[0] in ("log", "sqrt") or (a[0] == "pow" and a[1].denominator != 1)
        or (a[0] == "pow" and a[1] < 0)
        for _, a in terms)
    if need_pos or cfg.domain_pos_only:
        return cfg.domain_lo, cfg.domain_hi
    return cfg.domain_neg_lo, cfg.domain_hi


def make_grid(lo: float, hi: float, n: int) -> list[float]:
    """扫描网格：正域用对数等分（1e-8 起步时线性网格几乎全压在大端），
    跨零域用线性等分。PROBE: 站端网格形态未钉死。"""
    if lo > 0:
        step = math.log(hi / lo) / (n - 1)
        return [lo * math.exp(step * i) for i in range(n)]
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def numeric_min(expr: sp.Expr, domain: tuple[float, float],
                cfg: ProbeConfig) -> tuple[float, float]:
    """网格取优 + f' 变号单元 brentq 求根 -> (argmin, min)。

    站端 x 给到 ~15 位有效数字，指向导数求根而非直接极小化；
    候选 = 两端点 + 网格最优点 + 每个 f' 变号网格单元的 brentq 根。
    """
    f = safe_func(sp.lambdify(X, expr, "math"))
    fp = safe_func(sp.lambdify(X, sp.diff(expr, X), "math"))
    xs = make_grid(*domain, cfg.grid_n)
    cands = [xs[0], xs[-1], min(xs, key=f)]
    for a, b in pairwise(xs):
        fa, fb = fp(a), fp(b)
        if math.isfinite(fa) and math.isfinite(fb) and fa * fb < 0:
            cands.append(brentq(fp, a, b, xtol=cfg.root_xtol))
    x = min(cands, key=f)
    return float(x), float(f(x))


def classify(terms: list[Term], domain: tuple[float, float], cfg: ProbeConfig) -> str:
    """convex|concave|affine —— 项级二阶导符号一致即可定；混合符号时按
    f'' 数值采样再判（PROBE: 站端是否做数值兜底、不定串如何标记均未实测）。"""
    pos = any(atom_curv(a) * (1 if c > 0 else -1) > 0 for c, a in terms)
    neg = any(atom_curv(a) * (1 if c > 0 else -1) < 0 for c, a in terms)
    if pos and not neg:
        return "convex"
    if neg and not pos:
        return "concave"
    if not pos and not neg:
        return "affine"
    d2 = sp.diff(sum_expr(terms), X, 2)
    _, mn = numeric_min(d2, domain, cfg)
    _, mx_neg = numeric_min(-d2, domain, cfg)
    if mn >= -cfg.curv_tol:
        return "convex" if -mx_neg > cfg.curv_tol else "affine"
    if -mx_neg <= cfg.curv_tol:
        return "concave"
    return "indefinite"  # PROBE: 不定情形站端如何标记未实测


# ============================ 切线搜索 ============================

def convergent_candidates(xmin: float, domain: tuple[float, float],
                          cfg: ProbeConfig) -> Iterator[Fraction]:
    """默认候选集: xmin 的连分数渐近分数（分母 ≤ cfg.candidate_den）。

    样本钉死: xmin≈0.567143290409784 的渐近分数依次为
    0, 1, 1/2, 4/7, 17/30, 38/67 —— 站端取到 17/30；分母更小的可证
    中间分数（如 13/23）未被采用，排除 Farey/limit_denominator 全枚举。
    PROBE: 真实候选集与排除规则待探测代理钉死，替换 cfg.candidates 即可。
    """
    lo, hi = domain
    rem = Fraction(xmin)
    p0, p1, q0, q1 = 0, 1, 1, 0
    while rem:
        a = rem.numerator // rem.denominator
        p, q = a * p1 + p0, a * q1 + q0
        if q > cfg.candidate_den:
            return
        x0 = Fraction(p, q)
        if lo < x0 < hi:
            yield x0
        p0, p1, q0, q1 = p1, p, q1, q
        rem -= a
        if rem:
            rem = 1 / rem


def slope_at(g_terms: list[Term], x0: Fraction) -> Fraction | None:
    """g'(x0) 精确值；非有理数（sqrt/pow 在一般切点无理）返回 None。"""
    total = Fraction(0)
    xs = sp.Rational(x0.numerator, x0.denominator)
    for c, a in g_terms:
        kind = a[0]
        if kind == "x":
            d = Fraction(1)
        elif kind == "log":
            d = 1 / x0
        else:
            v = sp.diff(atom_expr(a), X).subs(X, xs)
            if not v.is_Rational:
                return None
            d = Fraction(v.p, v.q)
        total += c * d
    return total


def atom_value_at(atom: Atom, x0: Fraction) -> Fraction | EvalAtom:
    """atom 在 x0 的取值：有理数折成 Fraction，否则返回字面求值原子。

    log 永不折叠（实测 `\\ln 1` 不化简 → 按字面保留）；sqrt/pow 取值为
    有理数时折叠进截距常数项（PROBE: 站端是否同样按字面显示未实测）。
    """
    kind = atom[0]
    if kind == "x":
        return x0
    if kind == "log":
        return ("logat", x0)
    if kind == "exp":
        return ("expat", x0)
    q = sp.Rational(x0.numerator, x0.denominator)
    v = sp.sqrt(q) if kind == "sqrt" else sp.Pow(
        q, sp.Rational(atom[1].numerator, atom[1].denominator))
    return Fraction(v.p, v.q) if v.is_Rational else (
        ("sqrtat", x0) if kind == "sqrt" else ("powat", x0, atom[1]))


def intercept_parts(g_terms: list[Term], m: Fraction, x0: Fraction
                    ) -> tuple[Fraction, list[tuple[Fraction, EvalAtom]]]:
    """切线截距 = -m·x0 + Σ c_i·atom_i(x0) -> (常数项, 字面显示项)。"""
    num = -m * x0
    items = []
    for c, a in g_terms:
        v = atom_value_at(a, x0)
        if isinstance(v, Fraction):
            num += c * v
        else:
            items.append((c, v))
    return num, items


def evalatom_sym(ea: EvalAtom) -> sp.Expr:
    """求值原子 -> sympy 精确表达式（gap 检验走符号导数求根）。"""
    kind = ea[0]
    q = sp.Rational(ea[1].numerator, ea[1].denominator)
    if kind == "logat":
        return sp.log(q)
    if kind == "sqrtat":
        return sp.sqrt(q)
    if kind == "expat":
        return sp.exp(q)
    return q ** sp.Rational(ea[2].numerator, ea[2].denominator)  # powat


def search_tangent(left: list[Term], g_terms: list[Term],
                   const: Fraction, domain: tuple[float, float], xmin: float,
                   cfg: ProbeConfig, left_tex: str, right_tex: str) -> dict | None:
    """枚举有理切点，取凹侧切线，数值验证 left - line > 0。

    实测行为: 候选的切线就是凹侧自身（仿射右侧切线 ≡ 自身，e^x>=x+1
    的零 gap 因此全灭 → "没找到漂亮的有理切点直线"）；站端只验证
    left-line 一侧的 gap（proof 无 right_gap 字段）。
    """
    lexpr = sum_expr(left)
    for x0 in cfg.candidates(xmin, domain, cfg):
        m = slope_at(g_terms, x0)
        if m is None:
            continue
        num, items = intercept_parts(g_terms, m, x0)
        b = sp.Rational((num + const).numerator, (num + const).denominator)
        for c, ea in items:
            b += sp.Rational(c.numerator, c.denominator) * evalatom_sym(ea)
        line_expr = sp.Rational(m.numerator, m.denominator) * X + b
        gx, gmin = numeric_min(lexpr - line_expr, domain, cfg)
        if gmin > cfg.gap_tol:
            line_text = line_render(m, num, items, const, False)
            line_tex = line_render(m, num, items, const, True)
            return {
                "formula_latex": f"{left_tex} > {line_tex} \\ge {right_tex}",
                "left_gap_min": gmin,
                "left_gap_min_x": gx,
                "line_latex": line_tex,
                "line_text": line_text,
                "tangent_at": str(x0),
                "tangent_at_latex": fmt_frac(x0, True),
            }
    return None


# ============================ provided_line ============================

def check_line(text: str, left: list[Term], right: list[Term],
               domain: tuple[float, float], cfg: ProbeConfig) -> dict:
    """隐藏字段 line: 同文法解析后归约成浮点 (m,b)，数值检验两侧 gap。"""
    m, b = Fraction(0), Fraction(0)
    for c, a in parse_sum(parse_expression(text)):
        if a[0] == "x":
            m += c
        elif a[0] == "const":
            b += c
        else:
            raise ValueError(LINE_NOT_AFFINE)
    mf, bf = float(m), float(b)
    line_expr = sp.Rational(m.numerator, m.denominator) * X + sp.Rational(
        b.numerator, b.denominator)
    gx, gmin = numeric_min(sum_expr(left) - line_expr, domain, cfg)
    rx, rmin = numeric_min(line_expr - sum_expr(right), domain, cfg)
    return {
        "b": bf,
        "left_gap_min": gmin,
        "left_gap_min_x": gx,
        "m": mf,
        "ok": gmin > -cfg.gap_tol and rmin > -cfg.gap_tol,  # PROBE: 判据未实测
        "right_gap_min": rmin,
        "right_gap_min_x": rx,
    }


# ============================ 主流程 ============================

@dataclass
class ProbeConfig:
    """站端数值机制的接缝参数；默认值是待探测钉死的假设（见模块 docstring）。"""

    domain_lo: float = 1e-8       # 实测: 含 log/sqrt 时从 x≈1e-08 起扫
    domain_hi: float = 10.0       # PROBE: 扫描上界未钉死
    domain_pos_only: bool = True  # PROBE: 无受限原子时是否仍 x>1e-8
    domain_neg_lo: float = -10.0  # PROBE: 全体实数时的左端
    grid_n: int = 400             # PROBE: 网格数与排布未钉死
    root_xtol: float = 2e-12      # brentq 的 xtol；影响最小值点末位 ulp（PROBE）
    min_tol: float = -1e-9        # fmin < min_tol -> failed（容忍边界 ~0 极小值）
    gap_tol: float = 1e-9         # 切线 gap > gap_tol 才采纳（零 gap 被拒实测）
    curv_tol: float = 1e-7        # 混合符号 f'' 数值分类的判零阈值
    candidate_den: int = 200      # PROBE: 切点候选渐近分数的分母上限
    candidates: Callable[..., Iterator[Fraction]] = convergent_candidates


CFG = ProbeConfig()


def prove(inequality: str, line: str | None = None, cfg: ProbeConfig = CFG) -> dict:
    """POST /convex/prove 的 result 字典；ValueError 由路由层映射 400。"""
    terms = parse_inequality(inequality)
    left, right = normalize(terms)
    domain = domain_of(terms, cfg)
    curvature = {
        "difference": classify(terms, domain, cfg),
        "left": classify(left, domain, cfg),
        "right": classify(right, domain, cfg),
    }
    normalized = {
        "difference": render_terms(terms, False),
        "difference_latex": render_terms(terms, True),
        "left": render_terms(left, False),
        "left_latex": render_terms(left, True),
        "right": render_terms(right, False),
        "right_latex": render_terms(right, True),
    }
    result = {
        "curvature": curvature,
        "minimum": None,
        "normalized": normalized,
        "ok": False,
        "proof": None,
        "provided_line": None,
        "reason": "",
        "status": "",
    }
    if line:
        result["provided_line"] = check_line(line, left, right, domain, cfg)

    if (curvature["left"] not in ("convex", "affine")
            or curvature["right"] not in ("concave", "affine")):
        result["status"] = "inconclusive"
        result["reason"] = REASON_INCONCLUSIVE
        return result

    xmin, fmin = numeric_min(sum_expr(terms), domain, cfg)
    result["minimum"] = {
        "value": fmin,
        "value_text": f"{fmin:.12g}",
        "x": xmin,
        "x_text": f"{xmin:.12g}",
    }
    if fmin < cfg.min_tol:
        result["status"] = "failed"
        result["reason"] = REASON_FAILED
        return result

    result["status"] = "proved"
    result["ok"] = True
    g_terms = [(c, a) for c, a in right if a[0] != "const"]
    const = sum((c for c, a in right if a[0] == "const"), Fraction(0))
    # left 仿射 / 右侧无常数项部分 -> 无中间直线可生成（PROBE: 后一子情形
    # 站端是否共用同一条 reason 未实测）
    if curvature["left"] == "affine" or not g_terms:
        result["reason"] = REASON_NO_LINE
        return result
    proof = search_tangent(left, g_terms, const, domain, xmin, cfg,
                           normalized["left_latex"], normalized["right_latex"])
    if proof is None:
        result["reason"] = REASON_SEARCH_MISS
    else:
        result["proof"] = proof
    return result
