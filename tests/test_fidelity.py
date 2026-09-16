"""Float64 方向预检的保真测试（solve.prove 层）。

站端机制（探测钉死，详见 docs/fidelity-notes.md）：进核前用 float64 求值
命题常数 c，把有理界与 c 做 **精确** 比较（Fraction vs float 的 Python 语义，
bound 侧不舍入）——'>' 命题 bound>c、'<' 命题 bound<c 即报"方向反了"。
float 等值但严格偏离 c 的界照样判反；bound==c（或 bound<C 真命题）越过预检
进扫描后，两方向的非正解耗尽都报"未找到解"，由 float64 真假兜底决定是否
改报"方向反了"。唯一例外是 zeta3 '>'，站端用 float(bound)>c。
下述有理界全部取自 bench/data/{edge,fidelity,probes}.jsonl 的实测记录。
"""

import math
from fractions import Fraction

import pytest

from attention_calculator import solve
from attention_calculator.engine import NoSolution, WrongDirection

F = Fraction


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        # 真命题、界比 c 高但 float 等值（c 低于真值）→ 站端仍报方向反了
        ("pi", "1", ">", "6134899525417045/1952799169684491"),  # π−4.9e-32
        ("e", "1", ">", "2124008553358849/781379079653017"),  # e−6.5e-32
        ("e_q", "1", ">", "2124008553358849/781379079653017"),  # 同界
        ("sin_q", "1", ">", "8199121568317184/9743795943467975"),  # sin1−9.6e-33
        ("ln_q", "2", ">", "1554903831458736/2243252046704767"),  # ln2−1.3e-32
        # ln_q_square 的 c 是 log(2)**2 二次舍入，界在其上 3.7e-17
        ("ln_q_square", "2", ">", "4646620020445232/9671330777074349"),
        # '<' 方向：float 等值但严格低于 c 的界（catalan 的 c 比真值高 1.1e-17）
        ("catalan", "1", "<", "915965594177219/1000000000000000"),
        # 各型 '<' bound<c 判反（探测记录原样）
        ("zeta3", "1", "<", "169174469410865/140737488355328"),  # Cf-1ulp
        ("e", "1", "<", "765128314358509/281474976710656"),  # Cf-1ulp
        ("pi", "1", "<", "7074237752028439/2251799813685248"),
        ("sin_q", "1", "<", "7579296827247853/9007199254740992"),
        ("e_q", "3", "<", "5653576037675269/281474976710656"),
        ("ln_q", "2", "<", "3121657384082679/4503599627370496"),
        ("arctan_q", "3", "<", "5625202075141471/4503599627370496"),
        ("pi_n", "5/2", "<", "2461979758178985/140737488355328"),
        ("e_pi", "1", "<", "6513525919879993/281474976710656"),
        ("varpi", "1", "<", "5904348712226991/2251799813685248"),
        ("gamma", "1", "<", "649887063340739/1125899906842624"),
        ("golden", "1", "<", "7286977268806823/4503599627370496"),
        ("gamma", "3", "<", "1"),
        # 系数型 power=0：c=0，'>0' 触发预检（0>0 本身不触发，见下）
        ("pi", "0", ">", "1"),
    ],
)
def test_direction_precheck_wrong_direction(kind, power, comp, bound):
    """bound 精确越过 float64 常数 → WrongDirection（无论 float 是否等值）。"""
    with pytest.raises(WrongDirection):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "kind,power,comp,bound",
    [
        # 界恰等于 c：'<' 越过预检进扫描，非正解按站端语义报未找到
        ("pi", "1", "<", "884279719003555/281474976710656"),  # == fl(π)
        ("e", "1", "<", "6121026514868073/2251799813685248"),  # == fl(e)
        ("catalan", "1", "<", "8250284617241437/9007199254740992"),  # == Cf
        ("zeta3", "1", "<", "5413583021147681/4503599627370496"),  # == Cf
        ("zeta3", "1", "<", "2706791510573841/2251799813685248"),  # Cf+1ulp
        # '>' 界低于 c（c 向上舍入的真命题紧界）：越过预检，搜索耗尽
        ("cos_q", "1", ">", "293104830616638/542483027433479"),  # cos1−3.2e-32
        ("arctan_q", "3", ">", "125014145208443/100087721339793"),  # atan3−7e-31
        # zeta3 '>' 是 float 判定：bound_f==Cf 的紧界与 Cf+0ulp 都放行耗尽
        ("zeta3", "1", ">", "5413583021147681/4503599627370496"),  # == Cf
        ("zeta3", "1", ">", "461424925/383862797"),  # Cf+4.3e-17
        # '>' 界恰等于 c（bound==Cf 的 dyadic）：预检不触发，扫描非正耗尽
        # 后 float 兜底 diff==0 → 站端"未找到>方向的解"而非"方向反了"
        ("cos_q", "1", ">", "1216652631687587/2251799813685248"),  # == fl(cos1)
        ("golden", "1", ">", "910872158600853/562949953421312"),  # == fl(φ)
        ("sin_q", "1/2", ">", "539785169252447/1125899906842624"),  # == fl(sin½)
        # trig_pi 无外层预检：bound<C 真命题、bound_f<Cf——defer 扫描只见
        # 伪非正计划（BIAS18 (1,8)），耗尽后按 float 兜底报"未找到"
        ("sin_pi_q", "1/5", ">", "587785252292473/1000000000000000"),
    ],
)
def test_direction_precheck_no_solution(kind, power, comp, bound):
    """未触发预检且搜索无符号定号解 → NoSolution。"""
    with pytest.raises(NoSolution):
        solve.prove(kind, power, comp, bound)


@pytest.mark.parametrize(
    "bound",
    [
        "2706791510573841/2251799813685248",  # Cf+1ulp
        "5413583021147683/4503599627370496",  # Cf+2ulp
    ],
)
def test_zeta3_gt_float_compare(bound):
    """zeta3 '>' 例外：bound_f 超过 Cf 才报方向反了（站端 float 判定）。"""
    with pytest.raises(WrongDirection):
        solve.prove("zeta3", "1", ">", bound)


def test_precheck_respects_domain_order():
    """域校验先于方向：ln_q q<=1、trig 值域外由核内 check_input 报错。"""
    with pytest.raises(ValueError):
        solve.prove("ln_q", "1", ">", "1")
    with pytest.raises(ValueError):
        solve.prove("ln_q_square", "0", "<", "1")
    with pytest.raises(ValueError):
        solve.prove("sin_q", "0", "<", "1")
    with pytest.raises(ValueError):
        solve.prove("sin_q", "355/113", ">", "0")
    with pytest.raises(ValueError):
        solve.prove("tan_q", "355/226", "<", "100")


def test_precheck_skips_negative_bound():
    """负界在站端被格式校验挡掉；核内 check_input 负责报错，不进预检。"""
    with pytest.raises(ValueError):
        solve.prove("ln_q", "3", "<", "-1/2")


def test_precheck_power_zero_zero():
    """power=0 且界=0：'>' 的 0>0 不触发预检，核给出退化 ∫0dx 证明（站端 200）。"""
    out = solve.prove("pi", "0", ">", "0")
    assert out["parameters"]["a_val"] == "0"
    assert out["parameters"]["b_val"] == "0"


def test_c_values():
    """钉死的各型 float64 常数：naive 表达式，不一定是正确舍入。"""
    assert solve.direction_f("ln_q_square", F(2)) == math.log(2) ** 2
    assert solve.direction_f("e_pi", F(1)) != math.exp(math.pi)  # 正确舍入
    assert solve.direction_f("e_q", F(3)) == math.exp(3.0)
    assert solve.direction_f("pi_n", F(5, 2)) == math.pi**2.5
    assert solve.direction_f("ln_q", F(1)) is None  # 域外
    assert solve.direction_f("ln_q_square", F(5)) is None  # 核内特判
    assert solve.direction_f("sin_q", F(355, 113)) is None  # 值域外
    assert solve.direction_f("artanh_q", F(1, 2)) is None  # 未探测，不预检
