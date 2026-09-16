"""tests for the mode=exact correctness judge (bench/judge_correct.py, W2).

Small deterministic corpus through judge_case: verdicts against certified
ground truth, zero-tolerance invariant (every emitted proof passes
exact_check), and the exact-vs-site divergence clusters.
"""

import importlib.util
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent.parent / "bench"


def load_judge():
    spec = importlib.util.spec_from_file_location("judge_correct", BENCH / "judge_correct.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


judge = load_judge()


def rec(kind, power, comp, rational):
    return judge.judge_case(
        {"type": kind, "power": power, "comparison": comp, "rational": rational}
    )


CLEAN = [
    # (type, power, comp, rational, expected verdict)
    ("e", "1", ">", "8/3", "proved"),
    ("e", "1", "<", "8/3", "rejected-false"),
    ("pi", "1", ">", "3", "proved"),
    ("pi", "1", "<", "22/7", "proved"),
    ("pi", "1", ">", "22/7", "rejected-false"),
    ("sin_pi_q", "1/6", "<", "1/2", "equal-claim"),  # Niven 点 sin(π/6)=1/2
    ("sin_pi_q", "1/6", ">", "1/2", "equal-claim"),
    ("ln_q", "1/2", "<", "0", "rejected"),  # 域外 q<=1
    ("pi", "1", "<", "22/7/2", "rejected"),  # 畸形输入
    # 站端奇异矩系统崩溃簇：exact 修复后正常证出
    ("ln_q_square", "5", "<", "13/5", "proved"),
    # 首轮扫出并已由上游修复的三个缺陷（原 KNOWN_GAPS，修复后收紧）：
    ("sinh_q", "3/5", ">", "311853313088795/489832024562742", "proved"),  # float64 预检门控
    ("e_q", "1/2", ">", "0", "proved"),  # exp_family target 零键
    ("pi_n", "0", ">", "0", "rejected"),  # 指数 0 超出矩表 [1,10]，干净域拒
]


@pytest.mark.parametrize(
    ("kind", "power", "comp", "rational", "want"), CLEAN, ids=[c[0] for c in CLEAN]
)
def test_verdicts(kind, power, comp, rational, want):
    assert rec(kind, power, comp, rational)["verdict"] == want


def test_every_emitted_proof_verifies():
    """零容忍不变量：小型混合语料上不允许任何 BUG:/FLAG:。"""
    corpus = [
        {"type": k, "power": p, "comparison": c, "rational": r}
        for k, p, c, r, _ in CLEAN
        if not r.endswith("/2/2")  # 跳过畸形输入（本身 rejected，无证明可验）
    ]
    for case in corpus:
        out = rec(case["type"], case["power"], case["comparison"], case["rational"])
        assert not out["verdict"].startswith(("BUG:", "FLAG:")), out
        if out["exact"]["outcome"] == "emitted":
            assert out["exact"]["identity_ok"] and out["exact"]["nonneg"]


def test_direction_oracle_crosscheck():
    """certified_cmp 与独立真值实现不得在同一条目上分歧。"""
    out = rec("pi", "1", "<", "22/7")
    assert "truth_mismatch" not in out


def test_site_divergence_clusters():
    """已知站端失真簇的归因：ln²q 奇异系统 500 → exact 证出。"""
    out = rec("ln_q_square", "5", "<", "13/5")
    assert out["site"]["outcome"] == "internal-error"
    assert out["divergence"] == "ln2-singular-500"
    # float64 方向预检在紧界上误判：exact 证出而站端误报方向反了
    out = rec("arctan_q", "1/5", "<", "197395559849880759/1000000000000000000")
    assert out["verdict"] == "proved"
    assert out["site"]["outcome"] == "wrong-direction"
    assert out["divergence"] == "float64-direction"


def test_adversarial_corpus_deterministic():
    gen = judge.generate_cases
    a, b = gen(), gen()
    assert a == b and len(a) > 1000
    assert all(set(c) == {"type", "power", "comparison", "rational"} for c in a)
