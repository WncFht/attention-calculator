# Benchmark：数据集、判官与评测口径

benchmark 是本项目一切结论的依据。**任何实现决策（核族参数、搜索顺序、失败边界）都要落到数据上的数字**。

`bench/data/` 全部 gitignored（站端采集语料，体积大，不入库——rsync/手动同步）； `bench/probes/*.json` 入库——是 tests/ 里内联断言的引用源；`bench/out/` 是评测输出。

## 数据集 ↔ 判官对照表

| 数据文件                            | 条数 | 采集/生成                       | 判官                                 | 口径                                                 |
| ----------------------------------- | ---- | ------------------------------- | ------------------------------------ | ---------------------------------------------------- |
| `golden.jsonl`                      | 3454 | `harvest.py`（`cases.py` 例表） | `parity.py`                          | /calculate + /get_integral_image 状态 + 响应体逐字节 |
| `golden.jsonl`                      | 同上 | 同上                            | `verify.py`                          | 恒等式数学真值（mpmath 50dps），站点 bug 复现计入假  |
| `combo.jsonl`                       | 87   | `probe_decompose.py`            | `parity_decompose.py`                | /decompose_inequality 完整响应逐字节                 |
| `decompose.jsonl`                   | 182  | 同上                            | 同上                                 | 逐题 JSON 级（transport_error 只需 400+error 键）    |
| `edge-probes.jsonl`                 | 376  | `edge_probe{,2,3,3b}.py`        | `replay_edge.py`                     | 边缘探针离线重放（不经网络）                         |
| `fuzz-probes.jsonl`                 | 922  | `fuzz.py --seed 1`              | `replay_fuzz.py`                     | 随机差分探针离线重放                                 |
| `fidelity-probes.jsonl`             | 228  | `fidelity_probe{,2,3,4}.py`     | `replay_capture.py`                  | float 保真边界探针重放                               |
| `probes.jsonl`                      | 167  | `probe_api.py` 等散采           | `replay_capture.py`                  | 160 吻合 + 7 站端瞬态 500 按 SKIP-TRANSIENT 跳过     |
| `health-probes.jsonl`               | 420  | `probe_health.py`               | `parity_health.py`                   | /health/calculate 逐字节                             |
| `convex-probes.jsonl`               | 333  | `probe_convex.py`               | `parity_convex.py`                   | /convex/prove 逐字节                                 |
| `convex-golden.jsonl`               | 323  | convex 采集                     | —（孤儿）                            | convex-behavior.md 建模用子集，无判官消费            |
| `site-*.html`                       | 8 页 | curl 抓取                       | tests/test_convex.py + test_pages.py | 页面逐字节 parity（缺席即 skip）                     |
| `decompose/`、`*.log`、`summary.md` | —    | 采集副产物                      | —                                    | 溯源档案                                             |

全部判官当前**零分歧**。一键全量：`bench/run_all.sh`。

## golden.jsonl 格式

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

## 评测指标（parity）

1. **成功率**：同一 case 我们 solver 是否给出证明（按服务端 success 对齐）。
2. **证明有效**：用 `bench/verify.py` 对我们的输出做 (a) 数值验证恒等式成立（mpmath 高精度积分）(b) 被积函数不变号。
3. **一致性**（次要）：参数/积分式与 golden 完全一致的比例——不要求完全一致，但一致率高说明搜索序复现对了。

## 目录组织

- 顶层 `*.py` = 活跃判官（`parity*.py`、`replay_*.py`、`verify.py`、`report.py`）
    - 语料生产者（`harvest.py`、`cases.py`、`fuzz.py`、`probe_*.py`、`edge_probe*.py`、 `fidelity_probe*.py`）+ `cf_bounds.py`（CF 界生成器）。
- `archive/` = 已收敛的一次性分析脚本（sim_\*/trace_\* 假设甄别、rule_search、analyze_decomp、`convex_model.py`/`decompose_model.py` 离线模型）——结论已钉进 src/ 与 tests/，保留作溯源，不再运行。
