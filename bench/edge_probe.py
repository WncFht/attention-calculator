r"""边缘探针：对 zhuyidao.net 发起 golden 网格未覆盖的请求，并本地重放比对。

每条探针先打站端（限速 ≥1.05s），再用 app.test_client() 重放同一请求，
追加写入 bench/data/edge-probes.jsonl：

    {"id", "endpoint", "method", "kind", "request", "note",
     "http_status", "raw", "ours_status", "ours_body",
     "match", "byte_match"}

match = 状态码一致且响应体语义一致（JSON 按解析后比较，HTML/纯文本按原文）；
byte_match = 响应体逐字节一致（站端 JSON 键排序 + \uXXXX 转义与本方不同，
属已知系统性差异，单独记录不计入 match）。

用法：.venv/bin/python bench/edge_probe.py [--start N] [--only PATTERN]
断点续跑：已写入的 id 自动跳过。
"""

import argparse
import contextlib
import json
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

BASE = "https://zhuyidao.net"
OUT = Path(__file__).parent / "data" / "edge-probes.jsonl"
MIN_INTERVAL = 1.05
TIMEOUT = 100

SESSION = requests.Session()
SESSION.headers["User-Agent"] = "attention-calculator-bench/0.1 (edge probes)"

# /get_integral_image 的合法基线参数（pi 1 < 22/7 的站端解）
IMG = {
    "m": "3",
    "n": "3",
    "a_val": "47/120",
    "b_val": "-13/120",
    "c_val": "0",
    "u_val": "120",
    "au_val": "47",
    "bu_val": "-13",
    "cu_val": "0",
    "type": "pi",
    "coef": "1",
    "comparison": "<",
    "rational": "22/7",
}

PROBES: list[dict] = []


def calc(pid, note, **form):
    """POST /calculate 表单探针。字段值为 None 表示该字段不发送。"""
    fields = {k: v for k, v in form.items() if v is not None}
    PROBES.append(
        {
            "id": pid,
            "endpoint": "/calculate",
            "method": "POST",
            "kind": "form",
            "request": fields,
            "note": note,
        }
    )


def img(pid, note, **query):
    """GET /get_integral_image 探针；缺省字段用 IMG 基线补齐。"""
    fields = dict(IMG)
    fields.update({k: v for k, v in query.items()})
    fields = {k: v for k, v in fields.items() if v is not None}
    PROBES.append(
        {
            "id": pid,
            "endpoint": "/get_integral_image",
            "method": "GET",
            "kind": "query",
            "request": fields,
            "note": note,
        }
    )


def raw(pid, endpoint, method, kind, request, note):
    PROBES.append(
        {
            "id": pid,
            "endpoint": endpoint,
            "method": method,
            "kind": kind,
            "request": request,
            "note": note,
        }
    )


# ---------------------------------------------------------------- 1. 二者相等
# golden 只钉了三个 Niven 点；下面这些类型的常数在特定 power 下恰好取有理数值，
# 或界与常数在 float64 下相等——探测站端"相等"判定的覆盖面与精度。

calc(
    "equal:cos_deg60-lt",
    "cos60°=1/2 精确（Niven 点，golden 未含 cos_q_degree 60）",
    type="cos_q_degree",
    power="60",
    comparison="<",
    rational="1/2",
)
calc(
    "equal:cos_deg60-gt",
    "同上反方向",
    type="cos_q_degree",
    power="60",
    comparison=">",
    rational="1/2",
)
calc(
    "equal:pi_n0-lt",
    "pi^0=1 vs 界1：pi_n power=0 是否触发相等",
    type="pi_n",
    power="0",
    comparison="<",
    rational="1",
)
calc("equal:pi_n0-gt", "同上反方向", type="pi_n", power="0", comparison=">", rational="1")
calc(
    "equal:pi_n0-lt2",
    "pi^0=1<2 为真：若无相等判定会进入搜索",
    type="pi_n",
    power="0",
    comparison="<",
    rational="2",
)
calc(
    "equal:e_q0-lt",
    "e^0=1 vs 界1（e_q 0 已知 500，验证是否先到相等）",
    type="e_q",
    power="0",
    comparison="<",
    rational="1",
)
calc("equal:e_q0-gt", "同上反方向", type="e_q", power="0", comparison=">", rational="1")
calc(
    "equal:cosh0-lt",
    "cosh0=1 vs 界1（本方 1/q 崩溃→500；站端可能先判相等）",
    type="cosh_q",
    power="0",
    comparison="<",
    rational="1",
)
calc("equal:cosh0-gt", "同上反方向", type="cosh_q", power="0", comparison=">", rational="1")
calc("equal:cosh0-lt2", "cosh0=1<2 为真", type="cosh_q", power="0", comparison="<", rational="2")
calc("equal:sinh0-gt0", "sinh0=0 vs 界0", type="sinh_q", power="0", comparison=">", rational="0")
calc("equal:sinh0-lt0", "0<0 为假", type="sinh_q", power="0", comparison="<", rational="0")
calc("equal:tanh0-gt0", "tanh0=0 vs 0", type="tanh_q", power="0", comparison=">", rational="0")
calc(
    "equal:coth0-lt2",
    "coth0 极点：500 还是其它",
    type="coth_q",
    power="0",
    comparison="<",
    rational="2",
)
calc(
    "equal:atan0-gt0",
    "arctan0=0 vs 0（已知 atan0→500）",
    type="arctan_q",
    power="0",
    comparison=">",
    rational="0",
)
calc(
    "equal:atan0-lt1", "arctan0=0<1 为真", type="arctan_q", power="0", comparison="<", rational="1"
)
calc(
    "equal:acot0-lt2",
    "arccot0=pi/2<2 为真，q=0 退化核",
    type="arccot_q",
    power="0",
    comparison="<",
    rational="2",
)
calc(
    "equal:acot0-gt1",
    "arccot0=pi/2>1 为真",
    type="arccot_q",
    power="0",
    comparison=">",
    rational="1",
)
calc(
    "equal:pi0-lt0",
    "0·pi=0 vs 界0：零目标恒等式",
    type="pi",
    power="0",
    comparison="<",
    rational="0",
)
calc("equal:pi0-gt0", "0>0 为假但两侧皆 0", type="pi", power="0", comparison=">", rational="0")
calc("equal:e0-lt0", "0·e vs 0", type="e", power="0", comparison="<", rational="0")
calc("equal:e0-gt0", "同上反方向", type="e", power="0", comparison=">", rational="0")
calc(
    "equal:gamma0-lt1",
    "0·gamma<1 为真（本方 bound/abs(0) 除零→500）",
    type="gamma",
    power="0",
    comparison="<",
    rational="1",
)
calc("equal:gamma0-lt0", "0·gamma vs 0", type="gamma", power="0", comparison="<", rational="0")
calc("equal:golden0-gt0", "0·phi vs 0", type="golden", power="0", comparison=">", rational="0")
calc("equal:zeta0-lt1", "0·zeta3<1", type="zeta3", power="0", comparison="<", rational="1")
calc("equal:catalan0-lt1", "0·C<1", type="catalan", power="0", comparison="<", rational="1")
calc("equal:e_pi0-lt1", "0·e^pi<1", type="e_pi", power="0", comparison="<", rational="1")
calc("equal:varpi0-lt1", "0·varpi<1", type="varpi", power="0", comparison="<", rational="1")
calc("equal:gauss0-lt1", "0·G<1", type="gauss", power="0", comparison="<", rational="1")
calc(
    "equal:ln1-gt0",
    "ln1=0 vs 界0：值域拒绝与相等谁先",
    type="ln_q",
    power="1",
    comparison=">",
    rational="0",
)
calc(
    "equal:sinpi_half",
    "sin_pi_q 1/2 值域外（<1/2 严格）且 sin=1 vs 界1",
    type="sin_pi_q",
    power="1/2",
    comparison="<",
    rational="1",
)
calc(
    "equal:sinpi_half_unred",
    "同上但 power 写 500/1000（未约分边界）",
    type="sin_pi_q",
    power="500/1000",
    comparison="<",
    rational="1",
)
calc(
    "equal:sindeg90",
    "sin_q_degree 90 值域外且 sin90°=1 vs 界1",
    type="sin_q_degree",
    power="90",
    comparison="<",
    rational="1",
)
calc(
    "equal:cosdeg0",
    "cos_q_degree 0 值域外且 cos0°=1 vs 界1",
    type="cos_q_degree",
    power="0",
    comparison="<",
    rational="1",
)

# float64 与常数恰好相等的界：站端相等判定若在 float64 上做则报"二者相等"，
# 若在精确有理数上做则按真假方向走搜索。
calc(
    "equal:pi-float-lt",
    "界 = float(pi) 精确有理值：测相等判定的精度",
    type="pi",
    power="1",
    comparison="<",
    rational="884279719003555/281474976710656",
)
calc(
    "equal:pi-float-gt",
    "同上反方向",
    type="pi",
    power="1",
    comparison=">",
    rational="884279719003555/281474976710656",
)
calc(
    "equal:e-float-lt",
    "界 = float(e)",
    type="e",
    power="1",
    comparison="<",
    rational="6121026514868073/2251799813685248",
)
calc(
    "equal:e-float-gt",
    "同上反方向",
    type="e",
    power="1",
    comparison=">",
    rational="6121026514868073/2251799813685248",
)
calc(
    "equal:pi_n52-float",
    "界 = float(pi^{5/2})：升幂路径的浮点判定",
    type="pi_n",
    power="5/2",
    comparison=">",
    rational="4923959516357971/281474976710656",
)
calc(
    "equal:sin1-float",
    "界 = float(sin1)",
    type="sin_q",
    power="1",
    comparison=">",
    rational="3789648413623927/4503599627370496",
)

# ---------------------------------------------------------------- 2. 指数上限
# 10^16 上限内可表达的最紧界（误差 ~1e-31），远超已观测的最大 (m,n)，
# 用来逼出"在指数不超过N的范围内未找到"文案并确认各类型的 N。

calc(
    "cap:pi-gt-tight",
    "pi 真方向 4.9e-32 紧界：耗尽搜索应报 cap（30?）",
    type="pi",
    power="1",
    comparison=">",
    rational="6134899525417045/1952799169684491",
)
calc(
    "cap:pi-lt-tight",
    "pi 真方向 2.3e-31 紧界（< 方向）",
    type="pi",
    power="1",
    comparison="<",
    rational="5706674932067741/1816491048114374",
)
calc(
    "cap:e-gt-tight",
    "e 真方向 6.5e-32：e 的 cap 是 10 还是 30",
    type="e",
    power="1",
    comparison=">",
    rational="2124008553358849/781379079653017",
)
calc(
    "cap:e-lt-tight",
    "e < 方向紧界",
    type="e",
    power="1",
    comparison="<",
    rational="9581113603440437/3524694718233772",
)
calc(
    "cap:sin_q-tight",
    "sin1 真方向 9.6e-33：sin_q cap（已知 ≥7）",
    type="sin_q",
    power="1",
    comparison=">",
    rational="8199121568317184/9743795943467975",
)
calc(
    "cap:cos_q-tight",
    "cos1 真方向 3.2e-32：cos_q cap",
    type="cos_q",
    power="1",
    comparison=">",
    rational="293104830616638/542483027433479",
)
calc(
    "cap:e_q-tight",
    "e^1 真方向紧界：e_q cap（≥9）",
    type="e_q",
    power="1",
    comparison=">",
    rational="2124008553358849/781379079653017",
)
calc(
    "cap:e_q3-tight",
    "e^3 紧界",
    type="e_q",
    power="3",
    comparison=">",
    rational="9060589063489840/451100167157089",
)
calc(
    "cap:ln_q-tight",
    "ln2 紧界：ln_q cap/结构性失败",
    type="ln_q",
    power="2",
    comparison=">",
    rational="1554903831458736/2243252046704767",
)
calc(
    "cap:zeta3-tight",
    "zeta3 紧界：复验 cap 10 文案",
    type="zeta3",
    power="1",
    comparison=">",
    rational="5517909911604193/4590389936699688",
)
calc(
    "cap:atan3-tight",
    "arctan3 紧界：复验 cap 10",
    type="arctan_q",
    power="3",
    comparison=">",
    rational="125014145208443/100087721339793",
)
calc(
    "cap:ln_qsq-tight",
    "ln^2 2 紧界 1.1e-32",
    type="ln_q_square",
    power="2",
    comparison=">",
    rational="4646620020445232/9671330777074349",
)

# ---------------------------------------------------------------- 3. power=0 全型扫描
for t in [
    "pi",
    "e",
    "pi_n",
    "e_q",
    "ln_q",
    "ln_q_square",
    "sin_q",
    "cos_q",
    "tan_q",
    "cot_q",
    "sin_q_degree",
    "cos_q_degree",
    "sin_pi_q",
    "cos_pi_q",
    "arctan_q",
    "arccot_q",
    "sinh_q",
    "cosh_q",
    "tanh_q",
    "coth_q",
    "artanh_q",
    "arcoth_q",
    "gamma",
    "golden",
    "catalan",
    "zeta3",
    "e_pi",
    "varpi",
    "gauss",
]:
    calc(
        f"zero:{t}", f"power=0 全型扫描：{t} 0 < 1", type=t, power="0", comparison="<", rational="1"
    )
# 变体：真方向的 0<1 对退化核更有意义，加 '>' 反方向对照
for t in ["e", "golden", "catalan", "zeta3", "e_pi", "varpi", "gauss", "gamma"]:
    calc(f"zero-gt:{t}", f"{t} 0 > 1（假方向）", type=t, power="0", comparison=">", rational="1")
calc(
    "zero:power-0/5",
    "power=0/5 未约分零（sin_q 值域）",
    type="sin_q",
    power="0/5",
    comparison="<",
    rational="1/2",
)
calc("zero:power-00", "power=00 前导零", type="pi", power="00", comparison="<", rational="4")

# ---------------------------------------------------------------- 3b. 退化输入
calc(
    "deg:bound0-pi-gt",
    "界=0 真方向（golden 有 pi>0 成功样例）",
    type="pi",
    power="1",
    comparison=">",
    rational="0",
)
calc("deg:bound0-e-gt", "e>0 真方向界0", type="e", power="1", comparison=">", rational="0")
calc("deg:bound0-sin-gt", "sin1>0 界0", type="sin_q", power="1", comparison=">", rational="0")
calc("deg:bound0-ln-gt", "ln2>0 界0", type="ln_q", power="2", comparison=">", rational="0")
calc("deg:bound0-gamma-gt", "gamma>0 界0", type="gamma", power="1", comparison=">", rational="0")
calc(
    "deg:bound0-acot-gt",
    "arccot1=pi/4>0 界0",
    type="arccot_q",
    power="1",
    comparison=">",
    rational="0",
)
calc("deg:bound0-cosh-gt", "cosh1>0 界0", type="cosh_q", power="1", comparison=">", rational="0")

calc(
    "deg:cap-num",
    "rational 分子恰 10^16",
    type="pi",
    power="1",
    comparison="<",
    rational="10000000000000000",
)
calc(
    "deg:cap-num-1",
    "rational 分子 10^16-1（可通过格式，pi<它真）",
    type="pi",
    power="1",
    comparison="<",
    rational="9999999999999999",
)
calc(
    "deg:cap-den",
    "rational 分母恰 10^16",
    type="pi",
    power="1",
    comparison="<",
    rational="1/10000000000000000",
)
calc(
    "deg:cap-den-1",
    "分母 10^16-1（界≈0，pi<它假）",
    type="pi",
    power="1",
    comparison="<",
    rational="1/9999999999999999",
)
calc(
    "deg:cap-both",
    "分子分母都 10^16-1",
    type="pi",
    power="1",
    comparison="<",
    rational="9999999999999999/9999999999999999",
)
calc(
    "deg:cap-pow-num",
    "power 分子恰 10^16",
    type="pi",
    power="10000000000000000",
    comparison="<",
    rational="9999999999999999",
)
calc(
    "deg:cap-pow-num-1",
    "power 10^16-1：LHS=10^16·pi > 一切可表达界",
    type="pi",
    power="9999999999999999",
    comparison="<",
    rational="9999999999999999",
)
calc(
    "deg:cap-pow-den",
    "power 分母恰 10^16",
    type="pi",
    power="1/10000000000000000",
    comparison="<",
    rational="1",
)

calc("deg:pow-plus", "power 带前导 +", type="pi", power="+1", comparison="<", rational="4")
calc("deg:rat-plus", "rational 带前导 +", type="pi", power="1", comparison="<", rational="+4")
calc("deg:pow-space-l", "power 前导空格", type="pi", power=" 1", comparison="<", rational="4")
calc("deg:pow-space-t", "power 尾随空格", type="pi", power="1 ", comparison="<", rational="4")
calc("deg:rat-space-l", "rational 前导空格", type="pi", power="1", comparison="<", rational=" 4")
calc("deg:rat-space-t", "rational 尾随空格", type="pi", power="1", comparison="<", rational="4 ")
calc(
    "deg:rat-space-mid",
    "rational 中间空格 22 /7",
    type="pi",
    power="1",
    comparison="<",
    rational="22 /7",
)
calc("deg:rat-space-mid2", "rational 22/ 7", type="pi", power="1", comparison="<", rational="22/ 7")
calc("deg:pow-decimal", "power 小数 1.5", type="e_q", power="1.5", comparison="<", rational="5")
calc("deg:pow-decimal2", "power .5", type="e_q", power=".5", comparison="<", rational="5")
calc("deg:pow-decimal3", "power 1.0", type="e_q", power="1.0", comparison="<", rational="3")
calc(
    "deg:pow-lead0",
    "power=007 前导零（7pi 系数）",
    type="pi",
    power="007",
    comparison="<",
    rational="23",
)
calc(
    "deg:rat-lead0", "rational=022/7 前导零", type="pi", power="1", comparison="<", rational="022/7"
)
calc(
    "deg:pow-unreduced",
    "power=2/4 未约分（e_q q=1/2）",
    type="e_q",
    power="2/4",
    comparison="<",
    rational="2",
)
calc("deg:type-unicode", "type=π unicode", type="π", power="1", comparison="<", rational="4")
calc("deg:pow-unicode", "power=π", type="pi", power="π", comparison="<", rational="4")
calc("deg:rat-unicode", "rational=π", type="pi", power="1", comparison="<", rational="π")
calc(
    "deg:pow-fullwidth",
    "power=１２ 全角数字",
    type="pi",
    power="１２",
    comparison="<",
    rational="38",
)
calc(
    "deg:rat-fullwidth",
    "rational=２２/７ 全角数字",
    type="pi",
    power="1",
    comparison="<",
    rational="２２/７",
)
calc("deg:pow-negzero", "power=-0", type="pi", power="-0", comparison="<", rational="1")
calc("deg:rat-negzero", "rational=-0", type="pi", power="1", comparison="<", rational="-0")
calc("deg:pow-1/1", "power=1/1 未约分 1", type="e_q", power="1/1", comparison="<", rational="3")
calc(
    "deg:bigpower-e_q",
    "e_q power=100：e^100≈2.7e43 远超可表达界",
    type="e_q",
    power="100",
    comparison=">",
    rational="9999999999999999",
)
calc(
    "deg:bigpower-e_q-lt",
    "同上 < 方向（假）",
    type="e_q",
    power="100",
    comparison="<",
    rational="9999999999999999",
)
calc(
    "deg:bigpower-ln",
    "ln_q power=10^6：ln(1e6)=13.8<14",
    type="ln_q",
    power="1000000",
    comparison="<",
    rational="14",
)
calc(
    "deg:bigpower-pi_n",
    "pi_n power=20 >1（pi^20 真，核 ln^19）",
    type="pi_n",
    power="20",
    comparison=">",
    rational="1",
)
calc(
    "deg:pi_n-11",
    "pi_n power=11（超出 BETA/ETA 表 k≤10）",
    type="pi_n",
    power="11",
    comparison=">",
    rational="1",
)
calc(
    "deg:pi_n-1/2", "pi_n 1/2：pi^1 vs 界^2", type="pi_n", power="1/2", comparison="<", rational="2"
)
calc(
    "deg:domain-edge-sin",
    "sin_q q=355/113>pi：值域外紧贴 pi",
    type="sin_q",
    power="355/113",
    comparison="<",
    rational="1",
)
calc(
    "deg:domain-edge-sin2",
    "sin_q q=312689/99532<pi：值域内但 sin<0",
    type="sin_q",
    power="312689/99532",
    comparison=">",
    rational="0",
)
calc(
    "deg:domain-edge-tan",
    "tan_q q=355/226 与 pi/2 比较",
    type="tan_q",
    power="355/226",
    comparison="<",
    rational="100",
)
calc(
    "deg:domain-edge-cosdeg",
    "cos_q_degree 90 值域外",
    type="cos_q_degree",
    power="90",
    comparison="<",
    rational="1",
)
calc(
    "deg:domain-edge-sindeg89",
    "sin_q_degree 89 值域内贴边",
    type="sin_q_degree",
    power="89",
    comparison="<",
    rational="1",
)
calc(
    "deg:domain-edge-sindeg-frac",
    "sin_q_degree 89/2=44.5 分数度",
    type="sin_q_degree",
    power="89/2",
    comparison="<",
    rational="1",
)

# ---------------------------------------------------------------- 3c. 字段/方法
calc("fld:comp-empty", "comparison 字段缺失", type="pi", power="1", comparison=None, rational="4")
calc(
    "fld:comp-emptystr", "comparison='' 显式空串", type="pi", power="1", comparison="", rational="4"
)
calc("fld:power-emptystr", "power='' 显式空串", type="pi", power="", comparison="<", rational="4")
calc(
    "fld:rational-emptystr",
    "rational='' 显式空串",
    type="pi",
    power="1",
    comparison="<",
    rational="",
)
calc("fld:comp-ge", "comparison=≥ unicode", type="pi", power="1", comparison="≥", rational="4")
calc(
    "fld:comp-space",
    "comparison=' < ' 带空格",
    type="pi",
    power="1",
    comparison=" < ",
    rational="4",
)
calc("fld:comp-dbl", "comparison='<<'", type="pi", power="1", comparison="<<", rational="4")
calc("fld:type-PI", "type=PI 大写", type="PI", power="1", comparison="<", rational="4")
calc("fld:type-space-l", "type=' pi' 前导空格", type=" pi", power="1", comparison="<", rational="4")
calc("fld:type-space-t", "type='pi ' 尾随空格", type="pi ", power="1", comparison="<", rational="4")
calc("fld:type-empty", "type='' 显式空串", type="", power="1", comparison="<", rational="4")
calc(
    "fld:type-missing",
    "type 缺失：pi 缺省回退（e 参数当 pi 解）",
    power="1",
    comparison="<",
    rational="22/7",
)
calc(
    "fld:empty-post",
    "完全空表单",
)
calc("fld:only-type", "只有 type", type="pi")
calc("fld:no-rational", "缺 rational", type="pi", power="1", comparison="<")
calc("fld:no-power", "缺 power", type="pi", comparison="<", rational="4")
calc(
    "fld:extra-field",
    "多余字段应被忽略",
    type="pi",
    power="1",
    comparison="<",
    rational="22/7",
    foo="bar",
    extra="1",
)

raw(
    "mth:get-calculate",
    "/calculate",
    "GET",
    "query",
    {"type": "pi", "power": "1", "comparison": "<", "rational": "22/7"},
    "GET /calculate：站端 405 还是其它",
)
raw(
    "mth:get-calculate-slash",
    "/calculate/",
    "GET",
    "query",
    {"type": "pi", "power": "1", "comparison": "<", "rational": "22/7"},
    "GET /calculate/ 尾斜杠",
)
raw(
    "mth:post-calculate-query",
    "/calculate",
    "POST",
    "query",
    {"type": "e", "power": "1", "comparison": "<", "rational": "3"},
    "POST 但参数在 query string：站端读 args 还是只读 form",
)
raw(
    "mth:post-calculate-json",
    "/calculate",
    "POST",
    "json",
    {"type": "e", "power": "1", "comparison": "<", "rational": "3"},
    "POST JSON body：form 解析不到字段时的回退",
)
raw(
    "mth:post-calculate-slash",
    "/calculate/",
    "POST",
    "form",
    {"type": "pi", "power": "1", "comparison": "<", "rational": "22/7"},
    "POST /calculate/ 尾斜杠",
)
raw(
    "mth:post-image",
    "/get_integral_image",
    "POST",
    "form",
    IMG,
    "POST /get_integral_image：笔记记 500，本方路由只挂 GET→405",
)
raw(
    "mth:get-decompose",
    "/decompose_inequality",
    "GET",
    "query",
    {"problem": "pi<22/7"},
    "GET /decompose_inequality",
)
raw(
    "mth:post-calculate-dup",
    "/calculate",
    "POST",
    "form-dup",
    [("type", "pi"), ("type", "e"), ("power", "1"), ("comparison", "<"), ("rational", "22/7")],
    "重复 type 键：Flask 取第一个",
)

# ---------------------------------------------------------------- 4. image 边
img("img:baseline", "基线 200 锚点")
img("img:empty", "空 query", **{k: None for k in IMG})
img("img:m-alpha", "m=abc", m="abc")
img("img:m-neg", "m=-1 负指数渲染", m="-1")
img("img:m-float", "m=3.5", m="3.5")
img("img:u-zero", "u_val=0 除零", u_val="0")
img("img:u-neg", "u_val=-120", u_val="-120")
img(
    "img:gamma-deg",
    "gamma2 退化参数组（u_val=0,m=n=0,a=1/4,c=1）站端是否渲染",
    m="0",
    n="0",
    a_val="1/4",
    b_val="0",
    c_val="1",
    u_val="0",
    au_val="0",
    bu_val="0",
    cu_val="0",
    type="gamma",
    coef="2",
    comparison=">",
    rational="1/2",
)
img("img:type-foo", "type=foo 未知", type="foo")
img("img:type-missing", "type 缺失", type=None)
img("img:coef-missing", "coef 缺失（站端有默认值?）", coef=None)
img("img:rational-missing", "rational 缺失", rational=None)
img(
    "img:bogus-params",
    "合法格式但非真实解的参数组",
    m="99",
    n="99",
    a_val="1/2",
    b_val="3/4",
    c_val="5/6",
    u_val="7",
    au_val="1",
    bu_val="3",
    cu_val="5",
    rational="4",
)
img("img:comp-eq", "comparison='='：渲染还是报错", comparison="=")
img("img:comp-missing", "comparison 缺失", comparison=None)
img("img:rational-neg", "rational=-3", rational="-3")
img("img:rational-dec", "rational=3.5", rational="3.5")
img("img:coef-alpha", "coef=abc", coef="abc")
img("img:coef-unreduced", "coef=2/4 原文回显", coef="2/4")
img("img:coef-lead0", "coef=01 前导零回显", coef="01")
img("img:rational-lead0", "rational=022/007 回显", rational="022/007")
img("img:rational-latex", "rational=\\frac{22}{7} LaTeX 形（站端接受）", rational="\\frac{22}{7}")
img("img:rational-spaces", "rational=' 22/7 ' 带空格", rational=" 22/7 ")
img("img:a-zero-den", "a_val=1/0", a_val="1/0")
img("img:a-alpha", "a_val=abc", a_val="abc")
img("img:missing-cval", "缺 c_val 字段", c_val=None)
img("img:extra-param", "多余 query 字段", unused="zzz")

# ---------------------------------------------------------------- 5. unified_form
for t, pw, cp, rt in [
    ("pi", "3/2", "<", "5"),
    ("pi", "7", ">", "22"),
    ("e", "5/3", ">", "9/2"),
    ("e", "2", "<", "6"),
    ("pi_n", "7", "<", "3021"),
    ("e_q", "7/2", "<", "34"),
    ("ln_q", "20", "<", "3"),
    ("ln_q_square", "10", "<", "6"),
    ("sin_q", "5/4", ">", "9/10"),
    ("cos_q", "5/4", "<", "1/3"),
    ("tan_q", "5/4", ">", "3"),
    ("cot_q", "5/4", "<", "1/3"),
    ("sinh_q", "3", ">", "10"),
    ("cosh_q", "3", "<", "11"),
    ("tanh_q", "2", ">", "24/25"),
    ("coth_q", "1/2", ">", "2"),
    ("artanh_q", "1/5", "<", "1/4"),
    ("arcoth_q", "4", "<", "3/10"),
    ("arccot_q", "3", ">", "3/10"),
    ("golden", "3", "<", "5"),
    ("catalan", "2", "<", "2"),
    ("zeta3", "2", ">", "2"),
    ("e_pi", "3", ">", "69"),
    ("varpi", "2", ">", "5"),
    ("gauss", "2", ">", "3/2"),
    ("gamma", "3", "<", "2"),
    ("sin_q_degree", "45", ">", "7/10"),
    ("cos_pi_q", "2/5", "<", "2/5"),
]:
    calc(
        f"uf:{t}-{pw}-{cp}",
        f"unified_form 猎捕：{t} {pw} {cp} {rt}",
        type=t,
        power=pw,
        comparison=cp,
        rational=rt,
    )

# ---------------------------------------------------------------- 6. 页面
raw("page:root", "/", "GET", "none", {}, "首页 HTML")
raw("page:en", "/en", "GET", "none", {}, "英文页 HTML：与 / 的差异")
raw("page:en-slash", "/en/", "GET", "none", {}, "/en/ 尾斜杠")
raw("page:404", "/nonexistent-page", "GET", "none", {}, "404 页面形态")
raw("page:robots", "/robots.txt", "GET", "none", {}, "robots.txt")


# ---------------------------------------------------------------- 运行器


def site_request(p):
    """按探针描述发站端请求，返回 (status, raw_text)。"""
    if p["method"] == "GET":
        resp = SESSION.get(BASE + p["endpoint"], params=p["request"] or None, timeout=TIMEOUT)
    elif p["kind"] in ("form", "form-dup"):
        resp = SESSION.post(BASE + p["endpoint"], data=p["request"], timeout=TIMEOUT)
    elif p["kind"] == "query":
        resp = SESSION.post(BASE + p["endpoint"], params=p["request"], timeout=TIMEOUT)
    elif p["kind"] == "json":
        resp = SESSION.post(BASE + p["endpoint"], json=p["request"], timeout=TIMEOUT)
    else:
        raise ValueError(p["kind"])
    return resp.status_code, resp.text


def ours_request(client, p):
    """同一请求打到本方 test_client。"""
    if p["method"] == "GET":
        resp = client.get(p["endpoint"], query_string=p["request"] or None)
    elif p["kind"] in ("form", "form-dup"):
        resp = client.post(p["endpoint"], data=p["request"])
    elif p["kind"] == "query":
        resp = client.post(p["endpoint"], query_string=p["request"])
    elif p["kind"] == "json":
        resp = client.post(p["endpoint"], json=p["request"])
    else:
        raise ValueError(p["kind"])
    return resp.status_code, resp.get_data(as_text=True)


def bodies_equal(site_raw, ours_raw):
    """语义比较：两边都能解析成 JSON 就按 dict 比，否则按原文比。"""
    try:
        s, o = json.loads(site_raw), json.loads(ours_raw)
    except (json.JSONDecodeError, TypeError):
        return site_raw.strip() == ours_raw.strip()
    return s == o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()

    from attention_calculator.server import app

    client = app.test_client()

    done = set()
    if OUT.exists():
        for line in OUT.read_text().splitlines():
            if line.strip():
                with contextlib.suppress(json.JSONDecodeError):
                    done.add(json.loads(line)["id"])
    todo = [p for p in PROBES if p["id"] not in done]
    if args.only:
        todo = [p for p in todo if args.only in p["id"]]
    todo = todo[args.start :]
    print(f"{len(done)} done, {len(todo)} to probe", flush=True)

    last = 0.0
    for i, p in enumerate(todo, 1):
        wait = MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        t0 = time.monotonic()
        try:
            status, raw_body = site_request(p)
        except (OSError, requests.RequestException) as exc:
            status, raw_body = -1, f"REQUEST_FAILED: {exc}"
        elapsed = time.monotonic() - t0
        last = time.monotonic()

        try:
            ostatus, obody = ours_request(client, p)
        except Exception as exc:  # 本地崩溃也记录
            ostatus, obody = -2, f"LOCAL_CRASH: {type(exc).__name__}: {exc}"

        rec = dict(p)
        rec.update(
            {
                "http_status": status,
                "raw": raw_body,
                "elapsed_ms": round(elapsed * 1000, 1),
                "ours_status": ostatus,
                "ours_body": obody,
                "byte_match": raw_body == obody,
                "match": status == ostatus and bodies_equal(raw_body, obody),
            }
        )
        with OUT.open("a") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        flag = "" if rec["match"] else "  <<< MISMATCH"
        print(
            f"[{i}/{len(todo)}] {p['id']}: site={status} ours={ostatus} "
            f"({rec['elapsed_ms']:.0f}ms){flag}",
            flush=True,
        )


if __name__ == "__main__":
    main()
