"""组合不等式拆解：把"常数 ⋚ 有理数"的复合不等式拆成若干基础型子证明。

界值分配机制 vendor 自 bench/archive/decompose_model.py（对全部实测记录验证过）：
每个加法/乘积项分到一条记录界（chain of running-min/max 有理记录），
乘积项再按 FACTOR_RANK 把界拆到各因子，倒数因子取底数记录。
"""

import math
import re
from fractions import Fraction

import sympy as sp

from . import render, solve
from .engine import EqualClaim, NoSolution, WrongDirection

# ------------------------------------------------------------------ sympy 表示
# 常数一律用 Symbol（e/pi/gamma/golden/catalan/gauss/varpi/zeta3），canon 排序
# 与站端 normalized_latex 完全一致；函数也全部用 undefined Function——未定义
# 应用函数按「类创建序」排 canon，ln,sin,cos,tan,arctan,sinh,tanh,exp 的创建
# 顺序恰好逐字节复刻站端排序（如 -sin-cos-tan、e+sin(1) 的项序）。
E = sp.Symbol("e")
PI = sp.Symbol("pi")
GAMMA = sp.Symbol("gamma")
GOLDEN = sp.Symbol("golden")
CATALAN = sp.Symbol("catalan")
GAUSS = sp.Symbol("gauss")
VARPI = sp.Symbol("varpi")
ZETA3 = sp.Symbol("zeta3")
LN = sp.Function("ln")
SIN = sp.Function("sin")
COS = sp.Function("cos")
TAN = sp.Function("tan")
ARCTAN = sp.Function("arctan")
SINH = sp.Function("sinh")
TANH = sp.Function("tanh")
EXP = sp.Function("exp")

LOC = {
    "e": E,
    "pi": PI,
    "π": PI,
    "gamma": GAMMA,
    "phi": GOLDEN,
    "golden": GOLDEN,
    "catalan": CATALAN,
    "gauss": GAUSS,
    "varpi": VARPI,
    "zeta3": ZETA3,
    "ln": LN,
    "arctan": ARCTAN,
    "atan": ARCTAN,
    "exp": EXP,
    "sin": SIN,
    "cos": COS,
    "tan": TAN,
    "sinh": SINH,
    "tanh": TANH,
}

BASE = {
    PI: ("pi", Fraction(1)),
    E: ("e", Fraction(1)),
    GAMMA: ("gamma", Fraction(1)),
    GOLDEN: ("golden", Fraction(1)),
    CATALAN: ("catalan", Fraction(1)),
    ZETA3: ("zeta3", Fraction(1)),
    VARPI: ("varpi", Fraction(1)),
    GAUSS: ("gauss", Fraction(1)),
}

# 乘积项里只允许非函数型因子（ln(2)*e、sin(1)*pi 等报"基础常数"错）
FUNC_KINDS = {"ln_q", "sin_q", "cos_q", "tan_q", "arctan_q", "sinh_q", "tanh_q"}
# 系数型：prove 的 power 槽位填项系数；参数型填原子参数
COEF_KINDS = {"pi", "e", "e_pi", "gamma", "golden", "catalan", "zeta3", "varpi", "gauss"}
BASE_OF = {"pi_n": "pi", "e_q": "e"}

ERR_ATOM = "当前乘积证明只支持基础常数的乘积"
ERR_POW_SPECIAL = "当前特殊常数暂不支持乘方组合"
ERR_POW_INT = "当前乘积证明只支持整数次幂；e^pi 可作为盖尔方德常数使用，但其他常数作指数暂不支持"
ERR_NUMERIC = "经数值检验，该不等式不成立；请检查不等号方向或输入内容"
ERR_NO_SPLIT = "暂未找到可由基础注意力积分证明的加法乘积拆解"
ERR_RHS_FORM = "当前乘积组合只支持“乘积 与 有理数”比较"

BAD_CHARS = re.compile(r"[^0-9A-Za-z_()+\-*/^.\sπ]")

EXPONENT_LIMIT = {"pi": 30, "e": 30}


def parse_atom(expr):
    """sympy 因子 → (kind, arg)；不支持的形状抛对应站端错误串。"""
    if expr in BASE:
        return BASE[expr]
    if expr.is_Pow:
        b, e = expr.base, expr.exp
        if b == E and e == PI:
            return "e_pi", Fraction(1)
        if b in (PI, E):
            if e.is_Integer:
                return ("pi_n" if b == PI else "e_q"), Fraction(e.p, e.q)
            raise ValueError(ERR_POW_INT)
        if b in BASE:
            raise ValueError(ERR_POW_SPECIAL)
        raise ValueError(ERR_ATOM)
    if expr.is_Function:
        name = expr.func.__name__
        a = expr.args[0] if expr.args else None
        if name == "exp":
            if a == PI:
                return "e_pi", Fraction(1)
            if a is not None and a.is_Rational:
                return "e_q", Fraction(a.p, a.q)
            raise ValueError(ERR_ATOM)
        if a is not None and a.is_Rational:
            arg = Fraction(a.p, a.q)
            table = {
                "ln": "ln_q",
                "log": "ln_q",
                "sin": "sin_q",
                "cos": "cos_q",
                "tan": "tan_q",
                "arctan": "arctan_q",
                "atan": "arctan_q",
                "sinh": "sinh_q",
                "tanh": "tanh_q",
            }
            if name in table:
                return table[name], arg
        raise ValueError(ERR_ATOM)
    raise ValueError(ERR_ATOM)


def merge_atoms(atoms):
    """pi*pi -> pi_n(2)、e*e^2 -> e_q(3)：同族幂合并（多数情况 sympify 已合并）。"""
    fam = {"pi": "pi", "pi_n": "pi", "e": "e", "e_q": "e"}
    out = []
    for kind, arg in atoms:
        if out and fam.get(out[-1][0]) == fam.get(kind) == "pi" and kind != "e_pi":
            out.append(("pi_n", out.pop()[1] + arg))
        elif out and fam.get(out[-1][0]) == fam.get(kind) == "e":
            out.append(("e_q", out.pop()[1] + arg))
        else:
            out.append((kind, arg))
    out = [("pi", a) if k == "pi_n" and a == 1 else (k, a) for k, a in out]
    out = [("e", a) if k == "e_q" and a == 1 else (k, a) for k, a in out]
    return out


# ------------------------------------------------------------------ 数值与记录链

VARPI_V = float((sp.gamma(sp.Rational(1, 4)) ** 2 / (2 * sp.sqrt(2 * sp.pi))).evalf(40))
GAUSS_V = float((sp.gamma(sp.Rational(1, 4)) ** 2 / (2 * sp.sqrt(2 * sp.pi**3))).evalf(40))


def atom_value(kind, arg):
    """原子的 float 数值——只用于选记录界；界本身是 Fraction。"""
    f = float(arg)
    return {
        "pi": lambda: math.pi * f,
        "e": lambda: math.e * f,
        "pi_n": lambda: math.pi**f,
        "e_q": lambda: math.exp(f),
        "e_pi": lambda: math.e**math.pi * f,
        "ln_q": lambda: math.log(f),
        "sin_q": lambda: math.sin(f),
        "cos_q": lambda: math.cos(f),
        "tan_q": lambda: math.tan(f),
        "arctan_q": lambda: math.atan(f),
        "sinh_q": lambda: math.sinh(f),
        "tanh_q": lambda: math.tanh(f),
        "gamma": lambda: f * 0.5772156649015329,
        "golden": lambda: f * (1 + math.sqrt(5)) / 2,
        "catalan": lambda: 0.915965594177219 * f,
        "zeta3": lambda: 1.2020569031595942 * f,
        "varpi": lambda: VARPI_V * f,
        "gauss": lambda: GAUSS_V * f,
    }[kind]()


def term_true(coef, atoms):
    """项真值 coef·Πatom_value（float）。"""
    v = float(coef)
    for kind, arg in atoms:
        v *= atom_value(kind, arg)
    return v


class Chain:
    """某常数值的有理记录界：upper = ceil(nk)/k 递减记录，lower = floor 递增记录。"""

    def __init__(self, value, k_min, k_max=3000):
        """遍历 k∈[k_min,k_max] 收 running-min(upper)/running-max(lower) 记录。"""
        up, lo = [], []
        best_u, best_l = math.inf, -math.inf
        for k in range(k_min, k_max + 1):
            n = math.floor(value * k)
            if Fraction(n, k) > Fraction(str(value)):
                n -= 1
            u = Fraction(n + 1, k)
            if u < best_u:
                up.append(u)
                best_u = u
            lo_ = Fraction(n, k)
            if best_l < lo_ < value:
                lo.append(lo_)
                best_l = lo_
        self.upper = up
        self.lower = lo

    def largest_upper(self, cap):
        """不超过 cap 的最大上界记录。"""
        ok = [b for b in self.upper if b <= cap]
        return ok[0] if ok else None

    def largest_upper_strict(self, cap):
        """严格小于 cap 的最大上界记录。"""
        ok = [b for b in self.upper if b < cap]
        return ok[0] if ok else None

    def smallest_lower(self, floor):
        """严格大于 floor 的最小下界记录。"""
        ok = [b for b in self.lower if b > floor]
        return ok[0] if ok else None

    def smallest_lower_geq(self, floor):
        """不小于 floor 的最小下界记录。"""
        ok = [b for b in self.lower if b >= floor]
        return ok[0] if ok else None


CHAINS = {}


def chain(kind, arg, side):
    """记录链带惰性 k_min：站端记录表只收贴近真值的界——upper 须 ≤ 6v/5、
    lower 须 ≥ 5v/6，等价于跳过开头的粗记录（ln5 上界从 5/3 起、cos 下界从
    1/2 起均由此而来）。"""
    key = (kind, arg, side)
    if key not in CHAINS:
        v = atom_value(kind, arg)
        kmin = 1
        if side == "up":
            while math.ceil(v * kmin) / kmin > v * 1.2:
                kmin += 1
        else:
            while math.floor(v * kmin) < 1:
                kmin += 1
        CHAINS[key] = Chain(v, kmin)
    return CHAINS[key]


def val_chain(value, side):
    """任意项值的记录链（非单原子项走 k_min=1 全链）。"""
    key = ("v", round(value, 15), side)
    if key not in CHAINS:
        CHAINS[key] = Chain(value, 1)
    return CHAINS[key]


FACTOR_RANK = {
    "e": 0,
    "e_q": 0,
    "e_pi": 0,
    "pi": 1,
    "golden": 2,
    "gamma": 4,
    "ln_q": 5,
    "pi_n": 3,
    "zeta3": 9,
    "catalan": 8,
    "sin_q": 6,
    "cos_q": 6,
    "tan_q": 7,
    "arctan_q": 7,
    "sinh_q": 7,
    "tanh_q": 7,
    "varpi": 10,
    "gauss": 10,
}

# steps 序：乘积项在前（项间 canon 序），单原子项按此秩降序、同型按 arg 升序
STEP_RANK = {
    "golden": 100,
    "gamma": 90,
    "ln_q": 80,
    "sin_q": 70,
    "cos_q": 69,
    "tan_q": 60,
    "arctan_q": 59,
    "sinh_q": 58,
    "tanh_q": 57,
    "catalan": 56,
    "zeta3": 55,
    "varpi": 54,
    "gauss": 53,
    "e_q": 52,
    "pi_n": 51,
    "e_pi": 50,
    "pi": 10,
    "e": 0,
}


def is_reciprocal(atoms):
    """atoms 含负幂因子 → 倒数形 a·b^{-1}。"""
    return any(a < 0 for _, a in atoms)


def pick_processed(terms, ti, share, comp):
    """非 resid 项的（带符号）项界；None 表示分不到。"""
    coef, atoms = terms[ti]
    v = term_true(coef, atoms)
    if len(atoms) == 1:
        kind, arg = atoms[0]
        vv = abs(v)
        if abs(coef) == 1:
            ch_up = chain(kind, arg, "up")
            ch_lo = chain(kind, arg, "lo")
        else:
            ch_up = val_chain(vv, "up")
            ch_lo = val_chain(vv, "lo")
        if comp == "<":
            if coef > 0:
                return ch_up.largest_upper(share)
            b = ch_lo.smallest_lower_geq(-share)
            return -b if b is not None else None
        if coef > 0:
            return ch_lo.smallest_lower(share)
        b = ch_up.largest_upper(-share)
        return -b if b is not None else None
    pch_up = val_chain(abs(v), "up")
    pch_lo = val_chain(abs(v), "lo")
    if comp == "<":
        if coef > 0:
            return pch_up.largest_upper(share) or share
        b = pch_lo.smallest_lower_geq(-share)
        return -b if b is not None else -share
    if coef > 0:
        return pch_lo.smallest_lower(share) or share
    b = pch_up.largest_upper(-share)
    return -b if b is not None else -share


def resid_bound(terms, ti, resid, comp, has_product):
    """resid 项的（带符号）项界。"""
    coef, atoms = terms[ti]
    if comp == "<":
        if coef > 0 and not has_product:
            if len(atoms) == 1 and coef == 1:
                b = chain(atoms[0][0], atoms[0][1], "up").largest_upper(resid)
            else:
                b = val_chain(abs(term_true(coef, atoms)), "up").largest_upper(resid)
            if b is not None:
                return b
        return resid
    if coef >= 0:
        return resid
    if len(atoms) == 1 and coef == -1:
        b = chain(atoms[0][0], atoms[0][1], "up").largest_upper(-resid)
    else:
        b = val_chain(abs(term_true(coef, atoms)), "up").largest_upper(-resid)
    if b is not None:
        return -b
    return resid


def factor_split(terms, ti, bound_signed, comp):
    """乘积项界拆到因子：最低 FACTOR_RANK 因子拿记录，尾因子 verbatim。"""
    coef, atoms = terms[ti]
    n = len(atoms)
    trues = [atom_value(k, a) for k, a in atoms]
    order = sorted(range(n), key=lambda j: FACTOR_RANK.get(atoms[j][0], 50))
    out = {}
    for pos, j in enumerate(order):
        kind, arg = atoms[j]
        prod_others = Fraction(1)
        for k in range(n):
            if k == j:
                continue
            prod_others *= out[k] if k in out else Fraction(str(trues[k]))
        cap = bound_signed / (abs(coef) * prod_others)
        if pos < n - 1:
            if comp == "<":
                b = chain(kind, arg, "up").largest_upper(cap)
            else:
                b = chain(kind, arg, "lo").smallest_lower(cap)
            if b is None:
                b = cap
        else:
            b = cap
        out[j] = b
    return out


def reciprocal_split(terms, ti, comp, R, first_j):
    """a·b^{-1} vs R：两因子都拿记录，按显示序处理。"""
    coef, atoms = terms[ti]
    num_j = next(j for j, (k, a) in enumerate(atoms) if a > 0)
    den_j = next(j for j, (k, a) in enumerate(atoms) if a < 0)
    nk, na = atoms[num_j]
    dk, da = atoms[den_j]
    den_true = Fraction(str(atom_value(dk, -da)))
    num_true = Fraction(str(atom_value(nk, na)))
    out = {}
    if comp == ">":
        if first_j == num_j:
            nb = chain(nk, na, "lo").smallest_lower(R * den_true / abs(coef))
            db = chain(dk, -da, "up").largest_upper_strict(nb * abs(coef) / R) if nb else None
        else:
            db = chain(dk, -da, "up").largest_upper(num_true * abs(coef) / R)
            nb = chain(nk, na, "lo").smallest_lower(R * db / abs(coef)) if db else None
    else:
        if first_j == num_j:
            nb = chain(nk, na, "up").largest_upper(R * den_true / abs(coef))
            db = chain(dk, -da, "lo").smallest_lower(nb * abs(coef) / R) if nb else None
        else:
            db = chain(dk, -da, "lo").smallest_lower(num_true * abs(coef) / R)
            nb = chain(nk, na, "up").largest_upper(R * db / abs(coef)) if db else None
    if nb is None or db is None:
        return None
    out[num_j], out[den_j] = nb, db
    return out


def choose_resid(terms, order, comp):
    """resid 在 order 里的位置：'<' 有乘积→首个非乘积项；否则末个负系数项；'>'→首项。"""
    has_product = any(len(terms[t][1]) > 1 for t in order)
    if comp == "<":
        if has_product:
            for pos, t in enumerate(order):
                if len(terms[t][1]) == 1:
                    return pos
            return 0
        negs = [pos for pos, t in enumerate(order) if terms[t][0] < 0]
        if negs:
            return negs[-1]
        return 0
    return 0


# ------------------------------------------------------------------ 渲染

CONST_LABEL = {
    "pi": "\\pi",
    "e": "e",
    "gamma": "\\gamma",
    "golden": "\\phi",
    "catalan": "C",
    "gauss": "G",
    "varpi": "\\varpi",
    "zeta3": "\\zeta(3)",
    "e_pi": "e^{\\pi}",
}
FUNC_NAME = {
    "ln_q": "ln",
    "sin_q": "sin",
    "cos_q": "cos",
    "tan_q": "tan",
    "arctan_q": "arctan",
    "sinh_q": "sinh",
    "tanh_q": "tanh",
}


def frac(x):
    r"""站端 \dfrac 有理数渲染（分母 1 裸排；负号进分子）。"""
    x = Fraction(x)
    if x.denominator == 1:
        return str(x.numerator)
    return f"\\dfrac{{{x.numerator}}}{{{x.denominator}}}"


def atom_label(kind, arg):
    """step label / 单项渲染用的原子字形。"""
    if kind in CONST_LABEL:
        return CONST_LABEL[kind]
    if kind == "pi_n":
        return "\\pi" if arg == 1 else f"\\pi^{{{arg}}}"
    if kind == "e_q":
        return "e" if arg == 1 else f"e^{{{arg}}}"
    return f"\\{FUNC_NAME[kind]}\\left({arg}\\right)"


def factor_label(kind, arg):
    """乘积因子的字形：负幂因子 -> \\dfrac{1}{底数字形}。"""
    if arg < 0:
        return f"\\dfrac{{1}}{{{atom_label(BASE_OF.get(kind, kind), -arg)}}}"
    return atom_label(kind, arg)


def term_piece(coef, atoms):
    """decomposition_latex 左侧/单项 normalized 的项渲染。"""
    if len(atoms) == 1:
        kind, arg = atoms[0]
        if coef == 1:
            return atom_label(kind, arg)
        if coef == -1:
            return "-" + atom_label(kind, arg)
        return frac(coef) + atom_label(kind, arg)
    if coef == 1:
        head = ""
    elif coef == -1:
        head = "-"
    else:
        head = frac(coef) + "\\cdot "
    return head + "\\cdot ".join(factor_label(k, a) for k, a in atoms)


def join_signed(pieces):
    """用 + 连接各段；段首自带 - 时省略 +。"""
    out = pieces[0]
    for p in pieces[1:]:
        out += p if p.startswith("-") else "+" + p
    return out


def flip(comp):
    """不等号镜像。"""
    return "<" if comp == ">" else ">"


# ------------------------------------------------------------------ 主流程


def single_atom(expr):
    """expr 恰为单个基础原子（系数 1）→ (kind,arg)，否则 None。"""
    try:
        coef = Fraction(1)
        factors = []
        for f in sp.Mul.make_args(expr):
            if f.is_Rational:
                coef *= Fraction(f.p, f.q)
            else:
                factors.append(f)
        if coef != 1 or len(factors) != 1:
            return None
        atoms = merge_atoms([parse_atom(f) for f in factors])
        if len(atoms) != 1:
            return None
        return atoms[0]
    except (ValueError, TypeError):
        return None


def make_step(kind, arg, coef_disp, comp, bound, is_factor=False):
    """跑 prove + render 生成一个 step dict；失败包成站端 400 文案。

    is_factor：乘积/倒数拆出的因子证明——power 取原子参数、bound 原样下发
    （站端对 e^2/pi 里的 e² 直接证 e²-7，不除以幂次）。"""
    bound = Fraction(bound)
    if kind in COEF_KINDS or is_factor:
        power = str(coef_disp if kind in COEF_KINDS else arg)
        bound_eff = bound
    else:
        power, bound_eff = str(arg), bound / abs(coef_disp)
    lab = atom_label(BASE_OF.get(kind, kind), abs(arg)) if arg < 0 else atom_label(kind, arg)
    bl = frac(bound)
    try:
        res = solve.prove(kind, power, comp, str(bound_eff))
    except WrongDirection:
        inner = "要证明的式子不等号方向反了"
    except NoSolution:
        limit = EXPONENT_LIMIT.get(kind, 10)
        inner = f"在指数不超过{limit}的范围内未找到{comp}方向的解"
    except EqualClaim as exc:
        inner = str(exc)
    except ValueError as exc:
        inner = str(exc)
    except Exception as exc:
        inner = str(exc)
    else:
        equation = render.render_equation(res["parameters"], kind, power, comp, str(bound_eff))
        return {
            "bound": str(bound),
            "bound_latex": bl,
            "coefficient": str(coef_disp),
            "comparison": comp,
            "equation": equation,
            "label": lab,
            "type": kind,
        }
    raise ValueError(f"基础证明「{lab} {comp} {bl}」生成失败：{inner}")


def decompose_inequality(problem):
    """站端 /decompose_inequality 主流程；所有 400 走 ValueError。

    problem 字段原样回显的是去空白后的文本（站端同 /calculate 一样删
    " " 和 "\\n"）；sympify 失败在站端是 500 而非 400（SympifyError 恰是
    ValueError 子类，故在此换成 RuntimeError 让路由走 catch-all）。"""
    problem = problem.replace(" ", "").replace("\n", "")
    txt = problem.replace("°", "*pi/180").replace("^", "**")
    if txt.count(">") + txt.count("<") != 1:
        raise ValueError("请输入一个只包含一个 > 或 < 的不等式")
    comp = ">" if ">" in txt else "<"
    lt_s, rt_s = txt.split(comp)
    if not lt_s.strip() or not rt_s.strip():
        raise ValueError("表达式不能为空")
    if BAD_CHARS.search(lt_s) or BAD_CHARS.search(rt_s):
        raise ValueError("表达式包含暂不支持的字符")

    try:
        lhs = sp.expand(sp.sympify(lt_s, locals=LOC))
        rhs = sp.sympify(rt_s, locals=LOC)
    except Exception as exc:
        raise RuntimeError(f"sympify failed: {exc}") from exc

    # ---- 拆 lhs 成 (coef, atoms) 项；rational 项与乘积内函数因子按站端规则处理
    raw_terms = []  # (coef, factors)
    for term in sp.Add.make_args(lhs):
        coef = Fraction(1)
        factors = []
        for f in sp.Mul.make_args(term):
            if f.is_Rational:
                coef *= Fraction(f.p, f.q)
            else:
                factors.append(f)
        raw_terms.append((coef, factors))
    terms = []  # (coef, atoms)
    rationals = Fraction(0)
    for coef, factors in raw_terms:
        if not factors:
            rationals += coef
            continue
        atoms = merge_atoms([parse_atom(f) for f in factors])
        if len(atoms) > 1 and any(k in FUNC_KINDS for k, _ in atoms):
            raise ValueError(ERR_ATOM)
        terms.append((coef, atoms))

    # ---- rhs 归类
    rhs_numeric = rhs.is_Rational or rhs.is_Float
    ratom = None if rhs_numeric else single_atom(rhs)
    if not rhs_numeric and ratom is None:
        raise ValueError(ERR_RHS_FORM)

    # ---- 常数 ⋚ 常数：链式两证明
    if ratom is not None:
        if not (
            len(raw_terms) == 1 and len(terms) == 1 and len(terms[0][1]) == 1 and terms[0][0] == 1
        ):
            raise ValueError(ERR_RHS_FORM)
        lk, la = terms[0][1][0]
        rk, ra = ratom
        lv, rv = atom_value(lk, la), atom_value(rk, ra)
        if not (lv > rv if comp == ">" else lv < rv):
            raise ValueError(ERR_NUMERIC)
        if comp == ">":
            B = chain(lk, la, "lo").smallest_lower(Fraction(str(rv)))
        else:
            B = chain(lk, la, "up").largest_upper(Fraction(str(rv)))
        if B is None:
            raise ValueError(ERR_NO_SPLIT)
        steps = [
            make_step(lk, la, Fraction(1), comp, B),
            make_step(rk, ra, Fraction(1), flip(comp), B),
        ]
        nl = f"{atom_label(lk, la)}{comp}{atom_label(rk, ra)}"
        dl = f"{atom_label(lk, la)}{comp}{frac(B)}{comp}{atom_label(rk, ra)}"
        return {
            "basic_count": 2,
            "decomposition_latex": dl,
            "direct_basic": False,
            "normalized_latex": nl,
            "problem": problem,
            "steps": steps,
            "success": True,
        }

    R = Fraction(rhs.p, rhs.q) if rhs.is_Rational else Fraction(str(rhs))

    has_product = any(len(atoms) > 1 for _, atoms in terms)
    direct = len(raw_terms) == 1 and len(terms) == 1 and len(terms[0][1]) == 1

    if direct:
        # 单原子项不做数值预检，直接跑 prove——假命题的错误按站端被包成
        # 「基础证明…生成失败：要证明的式子不等号方向反了」400
        coef, atoms = terms[0]
        kind, arg = atoms[0]
        done = R if comp == "<" else resid_bound(terms, 0, R, comp, False)
        step = make_step(
            kind, arg, abs(coef), comp if coef > 0 else flip(comp), done if coef > 0 else -done
        )
        nl = term_piece(coef, atoms) + comp + frac(R)
        return {
            "direct_basic": True,
            "normalized_latex": nl,
            "problem": problem,
            "steps": [step],
            "success": True,
        }

    # ---- 拆解路径
    # 有 bug 的数值预检：倒数因子按底数值计（e/pi -> e·pi）；只对多项式子启用
    buggy = 0.0
    for coef, factors in raw_terms:
        v = float(coef)
        for f in factors:
            kind, arg = parse_atom(f) if len(factors) else None
            v *= atom_value(kind, abs(arg))
        buggy += v
    if not (buggy > float(R) if comp == ">" else buggy < float(R)):
        raise ValueError(ERR_NUMERIC)

    if has_product and rationals:
        raise ValueError(ERR_NO_SPLIT)
    if not has_product:
        R -= rationals
        terms = [(c, a) for c, a in terms]
    # 单原子负幂项（1/pi 这类）在加法里拆不动
    if any(len(atoms) == 1 and atoms[0][1] < 0 for _, atoms in terms):
        raise ValueError(ERR_NO_SPLIT)
    if not terms:
        raise ValueError(ERR_NO_SPLIT)

    order = [i for i, (_, a) in enumerate(terms) if len(a) > 1] + sorted(
        (i for i, (_, a) in enumerate(terms) if len(a) == 1),
        key=lambda i: (-STEP_RANK[terms[i][1][0][0]], terms[i][1][0][1]),
    )

    n = len(terms)
    done = {}
    if n == 1:
        ti = order[0]
        if is_reciprocal(terms[ti][1]):
            done[ti] = None
        elif comp == "<":
            done[ti] = R
        else:
            done[ti] = resid_bound(terms, ti, R, comp, False)
    else:
        resid_pos = choose_resid(terms, order, comp)
        resid_ti = order[resid_pos]
        proc = order[:resid_pos] + order[resid_pos + 1 :]

        # '<' 无乘积：站端在取反后的 '>' 形上分配记录界
        # （e-pi-gamma<0 实际按 gamma+pi-e>0 分配，resid 位置仍按 '<' 规则）
        flip_alloc = comp == "<" and not has_product
        aterms = [(-c, a) for c, a in terms] if flip_alloc else terms
        aR = -R if flip_alloc else R
        acomp = ">" if flip_alloc else comp

        def contrib(tj):
            """已分项取记录界，未分项回退项真值。"""
            if tj in done:
                return done[tj]
            return Fraction(str(term_true(*aterms[tj])))

        for ti in proc:
            share = aR - sum(contrib(tj) for tj in range(n) if tj != ti)
            b = pick_processed(aterms, ti, share, acomp)
            if b is None:
                raise ValueError(ERR_NO_SPLIT)
            done[ti] = b
        share_r = aR - sum(contrib(tj) for tj in range(n) if tj != resid_ti)
        done[resid_ti] = resid_bound(aterms, resid_ti, share_r, acomp, has_product)
        if flip_alloc:
            done = {t: -b for t, b in done.items()}

    # ---- steps + decomposition_latex
    steps = []
    sign_factor = -1 if (comp == "<" and not has_product) else 1
    lhs_pieces, bound_pieces, total = [], [], Fraction(0)
    for ti in order:
        coef, atoms = terms[ti]
        lhs_pieces.append(term_piece(coef * sign_factor, atoms))
        if len(atoms) == 1:
            kind, arg = atoms[0]
            b_atom = done[ti] if coef > 0 else -done[ti]
            step = make_step(kind, arg, abs(coef), comp if coef > 0 else flip(comp), b_atom)
            steps.append(step)
            sb = b_atom if coef * sign_factor > 0 else -b_atom
            bound_pieces.append("-" + frac(b_atom) if sb < 0 else frac(b_atom))
            total += sb
        else:
            if is_reciprocal(atoms):
                fsplit = reciprocal_split(terms, ti, comp, R, 0)
            elif done[ti] is None:
                fsplit = None
            else:
                fsplit = factor_split(terms, ti, abs(done[ti]), comp)
            if fsplit is None:
                raise ValueError(ERR_NO_SPLIT)
            parts = []
            pval = abs(coef)
            for j, (kind, arg) in enumerate(atoms):
                b = fsplit[j]
                skind = BASE_OF.get(kind, kind) if arg < 0 else kind
                scomp = comp if coef > 0 else flip(comp)
                if arg < 0:
                    scomp = flip(scomp)
                steps.append(make_step(skind, abs(arg), abs(arg), scomp, b, is_factor=True))
                parts.append(frac(1 / b) if arg < 0 else frac(b))
                pval *= b if arg > 0 else 1 / b
            if abs(coef) != 1:
                parts.insert(0, frac(abs(coef)))
            piece = "\\cdot ".join(parts)
            if coef * sign_factor < 0:
                piece = "-" + piece
                total -= pval
            else:
                total += pval
            bound_pieces.append(piece)
    lhs_dl = join_signed(lhs_pieces)
    bounds_dl = join_signed(bound_pieces)
    if comp == "<" and not has_product:
        # 取反后的 '>' 形：左式 = 移项符号项，尾链目标为 -R
        R_show, dl_comp = -R, ">"
    else:
        R_show, dl_comp = R, comp
    dl = f"{lhs_dl}{dl_comp}{bounds_dl}={frac(total)}"
    if total != R_show:
        dl += f"{dl_comp}{frac(R_show)}"

    # ---- normalized_latex
    if has_product:
        if len(terms) == 1:
            nl = term_piece(terms[0][0], terms[0][1]) + comp + frac(R)
        else:
            nl = sp.latex(lhs) + comp + frac(R)
    else:
        moved = (
            lhs - sp.sympify(rt_s, locals=LOC)
            if comp == ">"
            else sp.sympify(rt_s, locals=LOC) - lhs
        )
        nl = sp.latex(moved) + ">0"

    return {
        "basic_count": len(steps),
        "decomposition_latex": dl,
        "normalized_latex": nl,
        "problem": problem,
        "steps": steps,
        "success": True,
    }
