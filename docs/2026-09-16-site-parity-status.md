# 站点复刻里程碑状态报告（v1.0.0-site-parity，2026-09-16）

本报告冻结「字节级复刻 zhuyidao.net」阶段的终点状态，作为项目目标切换（→ 数学正确性）的基线。对应 git tag `v1.0.0-site-parity`（`f64c230`，master HEAD）。详细机制见 HANDOFF.md、kernel-spec.md、verify-report.md；本文件只做状态归档与失真面盘点。

## 交付物与判官基线

复刻对象是站点行为而非数学正确性：响应体逐字节一致是最高判据，站端 bug 原样保留。终态判官全绿：

| 判官 | 口径 | 数字 | 入口 |
|---|---|---|---|
| golden parity | 响应 status+body 逐字节；eq_match 为 /get_integral_image 逐字节 | 3454/3454；1588/1588 | `bench/parity.py` |
| decompose parity | combo 字节级 + decompose JSON 级 | 87/87 + 182/182 | `bench/parity_decompose.py` |
| health parity | 姊妹应用逐字节 | 420/420 | `bench/parity_health.py` |
| convex parity | 姊妹应用逐字节 | 333/333 | `bench/parity_convex.py` |
| edge 重放 | 边缘探针离线重放 | 376/376 | `bench/replay_edge.py` |
| fuzz 重放 | 随机探针离线重放 | 922/922 | `bench/replay_fuzz.py` |
| capture 重放 | fidelity 228 + probes 160（7 条站端瞬态 500 skip） | 全绿 | `bench/replay_capture.py` |
| 页面字节级 | 8 条路由 vs 抓取页面 | 全绿 | `tests/test_pages.py`、`test_convex.py` |
| pytest | 单测 + 上述判官包装 | 619 绿 + 21 skip（语料缺席） | `pytest tests/` |
| verify | 恒等式数学真值（50dps + acb + 精确矩三层） | 1438 valid / **146 false-identity** / 1866 方向自洽 | `bench/verify.py` |

全量评测 `bench/run_all.sh` 一条命令跑完。CI 含 format（ruff/shfmt/markdownlint/autocorrect/prettier/taplo/actionlint/gitleaks）与 test（uv sync + pytest）两个 job。

## 架构分层

`server.py`（Flask 路由 + 站点校验序 + JSON 序列化协议）→ `solve.py`（方向预检 + 类型分发 + 失败文案映射 `failure_text`）→ `engine.py`（(m,n) 搜索序、ℚ 上高斯消元、非负检查、WrongDirection/NoSolution 原语）→ `kernels/` 8 族 29 型（矩闭式/递推 + prove/render_equation）。`render.py` 逐字 LaTeX 排版（`\cdot`/空格怪癖逐型对齐）；`decompose.py` 组合不等式拆解；`convex.py` 站端 JS 逐字移植；`health.py` 纯表单算术；`integrand.py` 高精度常数重建（verify 专用，不入证明路径）。证明路径全 `Fraction` 精确算术，float 只做数值校验。

## 失真面盘点：站点 bug 清单（数学正确性阶段的修复清单）

复刻阶段刻意保留的全部数学错误，按机制分四簇。verify 在 3454 条 golden 上共核出 **146 条假恒等式**（恒等式等号不成立；不等式本身仍为真），全归属下述簇，无 UNKNOWN：

| 簇 | 条数 | 位置 | 机制 | 数学后果 |
|---|---|---|---|---|
| trig-bias | 30 | `kernels/trig_pi.py` `BIAS18_NUM` | (m,n)=(1,8) 档用含偏置 δ(α)·(C−1) 的存储矩公式解 (a,b)，33 点有理插值恢复 | 印刷恒等式两边差恰为 a·δ(α)·(1−C) |
| gauss-window | 3 | `kernels/beta_family.py` `transposed_lt_proof` | gauss '<' 在 (G, ~0.855] 窗口内走 m=5 转置错误方程组 | 发出的 (a,b) 不满足真矩方程，G⁻¹ 系数 ~1026 对声称 ~0.83 |
| 系数缩放 | ~113 | `solve.py` power≠1 路径 | power 乘进被积函数、LHS 声称 power·C−bound，真值核验 claimed = power×actual | 例：`gamma 2>1` 印 `2γ−1=∫2·kernel`，真积分为 γ−1/2 |
| ln_q_square 崩溃 | — | `kernels/log_family.py` q∈{5,7} 分流 | 站端求解必崩：真→500、假→404 方向反了 | 该输入域无证明产出 |

其余刻意复现但不产生假恒等式的行为：渲染端逐字回显不约分（`3140/1000` 照印）、`<` 方向扫描中 WrongDirection 映射为「未找到解」而非「方向反了」、varpi/gauss 兜底模板（cu_val=1/2，恒等式本身成立）、trig_pi 非正候选 defer 跳过、float64 方向兜底（zeta3/gamma 精度持平边界按不 falsify 处理）。

## 覆盖边界与已知缺口

- 29 型全部覆盖；文章目录另有 ~16 种休眠类型（arcsin/arccos、ln π、Γ(1/3)、Γ(1/4)、ψ'、erf、ζ 一般、Glaisher、Li₂、Si/Cin、Ein/Ei、arsinh/arcosh、(ln q)ⁿ 等），站端同样 400 拒收，复刻侧未实现——kernel-spec.md 已载其中多型的核与矩公式。
- 方向判定依赖 float64 预检/兜底：极端接近边界的命题（|C−r| 小于 float 分辨力）方向判定无认证保障，golden 内未遇到反例但理论上存在。
- verify 4 条 harness-error：varpi/gauss power=None 记录的重建缺口（verify 侧，非复现侧）。
- decompose 选界策略由 benchmark 反推钉死，机制上不能保证拆出的子界都可证。
- bench/data 为 gitignored 语料（rsync 同步），CI 上判官测试 skip。

## 后续

目标切换为数学正确性：产出**真**恒等式证明而非站点同款输出，并大规模扩展类型与命题覆盖。路线与工作流划分见 `docs/2026-09-16-math-correctness-plan.md`。
