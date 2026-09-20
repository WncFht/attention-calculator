# attention-calculator

_English: a reproduction and extension of the "attention calculator" at [zhuyidao.net](https://zhuyidao.net/) — given `constant ⋚ rational`, it automatically constructs a definite-integral identity proof. In `mode=exact` the identity is exact over ℚ and ships with a machine-checkable certificate._

复现并扩展 [zhuyidao.net](https://zhuyidao.net/)「注意力计算器」：输入"常数 ⋚ 有理数"，自动构造定积分恒等式证明。

项目分两阶段：

- **site 模式（已冻结，tag `v1.0.0-site-parity`）**：站端行为逐字节复刻——含站端 bug（错误的转置解、伪恒等式、崩溃不对称性）原样重现。范围 = 主站 29 种常数类型 + 姊妹应用 `/convex`、`/health`。判官全绿，基线与失真面盘点见 `docs/2026-09-16-site-parity-status.md`。
- **exact 模式（当前主线）**：数学正确性路径——恒等式在 ℚ 上精确成立、定号有精确证书、方向判定经递增精度认证；每条证明附机器可检 `certificate`；并扩展站端没有的 exact-only 类型（当前 26 型）。方案见 `docs/2026-09-16-math-correctness-plan.md`，当前状态 `docs/2026-09-16-exact-status.md`。

## 原理一句话

每种目标常数配一族被积函数 `B_{m,n}(x)·P(x)·K(x)`，使所有矩落在含目标常数的低维 ℚ-线性空间里；对 P 的 2~4 个系数在 ℚ 上解线性方程组（消掉寄生常数、常数项配齐目标界），再沿 (m,n) 网格（m+n 升序、|m−n| 小者优先）找到第一个使 P 在积分域上不变号的解。详见 `docs/kernel-spec.md`。

## 结构

- `src/attention_calculator/` — `engine.py`（搜索序/精确消元/定号判定）、`moment.py`（ℚ 字典矩向量：常数符号 ↦ Fraction 系数）、`kernels/`（20 族实现：站端 29 型 + exact-only 26 型）、`solve.py`（方向预检 + 分发 + `certified_cmp` 认证方向 + `prove_exact`）、`exact_check/`（19 族 ℚ 字典相等复核器，W0；gamma14/34/12 三型走 `certificate.py` 的证明 DAG 递归复核）、`certificate.py`（机器可检证书，W5）、`pade.py`（ln/arctan 的 Padé 插值第二证法，W4）、`agm.py`（gauss/varpi 的 AGM 区间第二证明器，W7）、`euler_gamma.py`（γ 的 Euler–Maclaurin 第二证明器，W7）、`decompose.py` + `decompose_exact.py`（组合式拆解：站点版 / 可证构造版）、`render.py`（LaTeX 逐字回显）、`server.py` + `templates/` + `static/`（Flask 站点）、`integrand.py`（常数高精值/被积函数重建）、`convex.py`、`health.py`（姊妹应用）
- `bench/` — 判官套件与站端采集语料：`parity*.py`（golden/decompose/health/convex 字节级）、`replay_*.py`（edge/fuzz/capture 离线重放）、`verify.py`（恒等式 50dps 数值审计，已被 exact_check 取代为正确性裁决）、`judge_correct.py` + `cases_correct.py`（exact 模式正确性判官 + 对抗语料，W2）、`run_all.sh`（一键全量）；`data/*.jsonl` 不入库（gitignored）
- `docs/` — `kernel-spec.md`（29 型数学规格 + exact-only 型清单，唯一事实源）、`api-spec.md`、`*-notes.md`（各域机制规格）、`2026-09-16-*.md`（exact 阶段方案/证书规格/W3 调研/里程碑状态）、`HANDOFF.md`（site 阶段完成记录）、`ROADMAP.md`（site 阶段过程记录）
- `tests/` — pytest（1724 过 + 21 跳：21 条全为 test_trig.py 参数拆分的有意互跳；另有语料缺席守卫，本机语料齐全不触发）
- `tools/`、`scripts/` — 离线复核（`tools/verify_cert.py`）、jsonl 统计（`tools/jlstat.py`）、markdown 整段化（`tools/unwrap_md.py`）、dev server 生命周期（`scripts/dev-server.sh`）、入库文件 lint 合集（`scripts/lint-tracked.sh`）等维护脚本

## 运行

```bash
uv sync                                # 或 pip install -e .
python -m attention_calculator.server  # waitress 起在 8080（等价于安装后的 attention-calculator 命令）
```

- `/`（及 `/en`、`/attention/*`）——与站端逐字节一致的页面；其前端引用 jsdelivr/hertzen 的 CDN，离线/不可达时公式不排版
- `/demo` —— 同页面、CDN 依赖换成本地 `static/vendor/` 副本；首次先跑 `scripts/fetch-vendor-assets.sh`（组件与许可证清单见 `scripts/THIRD-PARTY.md`）
- API：`POST /calculate`（表单加 `mode=exact` 走正确性路径，响应附 `certificate`）、`GET /get_integral_image`、`POST /decompose_inequality`（同认 `mode=exact`，走可证界分配）、`POST /convex/prove`、`POST /health/calculate`
- 证书离线复核：`.venv/bin/python tools/verify_cert.py cert.json`（或 stdin）

## 评测

```bash
.venv/bin/python bench/parity.py bench/data/golden.jsonl   # golden 3454/3454 字节级
.venv/bin/python bench/parity_decompose.py               # combo 87 + decompose 182
.venv/bin/python bench/parity_health.py                  # 420/420
.venv/bin/python bench/parity_convex.py                  # 333/333
.venv/bin/python bench/replay_edge.py                    # 376/376
.venv/bin/python bench/replay_fuzz.py                    # 922/922
.venv/bin/python bench/replay_capture.py                 # fidelity+probes 全绿
.venv/bin/python bench/verify.py bench/data/golden.jsonl # 恒等式数值审计（site 模式）
.venv/bin/python bench/judge_correct.py bench/data/golden.jsonl --adversarial  # exact 正确性判官
.venv/bin/python -m pytest tests/                        # 1724 过 + 21 跳
bench/run_all.sh                                         # 一键全量
```

站点判官当前零分歧；exact 判官口径见 `bench/README.md` 与 `docs/2026-09-16-math-correctness-plan.md`。
