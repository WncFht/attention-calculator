# Benchmark：golden 数据集与评测

benchmark 是本项目一切结论的依据。**任何实现决策（核族参数、搜索顺序、失败边界）都要落到 golden 数据上的数字**。

## 数据格式 `bench/data/golden.jsonl`

每行一个 case：

```json
{
    "type": "pi",
    "power": "1",
    "comparison": "<",
    "rational": "22/7",
    "success": true,
    "parameters": {
        "m": 3,
        "n": 3,
        "a_val": "47/120",
        "b_val": "-13/120",
        "c_val": "0",
        "au_val": 47,
        "bu_val": -13,
        "cu_val": 0,
        "u_val": 120
    },
    "equations": { "solution": "a = 47/120, b = -13/120" },
    "equation": "\\dfrac{22}{7} - \\pi = \\int_0^1 ... > 0",
    "elapsed_ms": 12.3
}
```

失败 case 记 `"success": false, "error": "<服务端原文>"`。

`combo.jsonl`：`/decompose_inequality` 的完整返回 + 原始 problem 字符串。

## 覆盖设计

- 29 个 type 全覆盖；带参数 q 的类型（ln_q、sin_q、…）取多个 q 值。
- 每个 (type, q)：围绕常数真值取若干有理界——连分数收敛子 + 不同间距的分数网格 + 反向（预期"方向反了"错误）。
- 目标规模 ≥ 1500 个 case；另加 ≥ 50 组合不等式。

## 评测指标（parity）

1. **成功率**：同一 case 我们 solver 是否给出证明（按服务端 success 对齐）。
2. **证明有效**：用 `bench/verify.py` 对我们的输出做 (a) 数值验证恒等式成立（mpmath 高精度积分）(b) 被积函数不变号。
3. **一致性**（次要）：参数/积分式与 golden 完全一致的比例——不要求完全一致，但一致率高说明搜索序复现对了。
