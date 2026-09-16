# attention-calculator

复现 [zhuyidao.net](https://zhuyidao.net/)「注意力计算器」：输入"常数 ⋚ 有理数"，自动构造定积分恒等式证明。复刻目标是**站点行为本身**——站端的 bug（错误的转置解、伪恒等式、崩溃不对称性）一律原样复现。范围 = 主站 29 种常数类型 + 姊妹应用 `/convex`（凹凸不等式）、`/health`（健康计算器）。

## 原理一句话

每种目标常数配一族被积函数 `B_{m,n}(x)·P(x)·K(x)`，使所有矩落在含目标常数的低维 ℚ-线性空间里；对 P 的 2~3 个系数在 ℚ 上解线性方程组（消掉寄生常数、常数项配齐目标界），再沿 (m,n) 网格（m+n 升序、|m−n| 小者优先）找到第一个使 P 在积分域上不变号的解。详见 `docs/kernel-spec.md`。

## 结构

- `src/attention_calculator/` — `engine.py`（搜索序/精确消元/定号判定）、`kernels/`（8 族实现 29 型）、`solve.py`（方向预检 + 分发）、`render.py`（LaTeX 逐字回显）、`decompose.py`（复合式拆解 + 记录链界值分配）、`server.py` + `templates/` + `static/`（Flask 站点）、`integrand.py`（常数高精值/被积函数重建）、`convex.py`、`health.py`（姊妹应用）
- `bench/` — 判官套件与站端采集语料：`parity*.py`（golden/decompose/health/convex 字节级）、`replay_*.py`（edge/fuzz/capture 离线重放）、`verify.py`（恒等式真值核验）、`run_all.sh`（一键全量）；`data/*.jsonl` 不入库（gitignored，rsync 同步）
- `docs/` — `kernel-spec.md`（29 型数学规格，唯一事实源）、`api-spec.md`、`decompose-notes.md`、`sibling-apps.md`、`fidelity-notes.md`、`verify-report.md`、`HANDOFF.md`（完成记录）等
- `tests/` — pytest（607 过 + 21 跳）
- `tools/`、`scripts/` — 预计算与维护脚本

## 运行

```bash
uv sync                                # 或 pip install -e .
python -m attention_calculator.server  # waitress 起在 8080
```

- `/`（及 `/en`、`/attention/*`）——与站端逐字节一致的页面；其前端引用 jsdelivr/hertzen 的 CDN，离线/不可达时公式不排版
- `/demo` —— 同页面、CDN 依赖换成本地 `static/vendor/` 副本；首次先跑 `scripts/fetch-vendor-assets.sh`
- API：`POST /calculate`、`GET /get_integral_image`、`POST /decompose_inequality`、`POST /convex/prove`、`POST /health/calculate`

## 评测

```bash
.venv/bin/python bench/parity.py bench/data/golden.jsonl   # golden 3454/3454 字节级
.venv/bin/python bench/parity_decompose.py               # combo 87 + decompose 182
.venv/bin/python bench/parity_health.py                  # 420/420
.venv/bin/python bench/parity_convex.py                  # 333/333
.venv/bin/python bench/replay_edge.py                    # 376/376
.venv/bin/python bench/replay_fuzz.py                    # 922/922
.venv/bin/python bench/replay_capture.py                 # fidelity+probes 全绿
.venv/bin/python bench/verify.py bench/data/golden.jsonl # 恒等式真值核验
.venv/bin/python -m pytest tests/                        # 607 过 + 21 跳
bench/run_all.sh                                         # 同步 devbox 后全量跑
```

全部判官当前零分歧；评测口径与语料清单见 `bench/README.md`。
