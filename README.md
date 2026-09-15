# attention-calculator

复现 [zhuyidao.net](https://zhuyidao.net/)「注意力计算器」：输入"常数 ⋚ 有理数"，自动构造定积分恒等式证明。

## 原理一句话

每种目标常数配一族被积函数 `B_{m,n}(x)·P(x)·K(x)`，使所有矩落在含目标常数的低维 ℚ-线性空间里；对 P 的 2~3 个系数在 ℚ 上解线性方程组（消掉寄生常数、常数项配齐目标界），再沿 (m,n) 网格（m+n 升序、|m−n| 小者优先）找到第一个使 P 在积分域上不变号的解。详见 `docs/kernel-spec.md`。

## 结构

- `src/attention_calculator/` — `moment.py`（精确矩向量）、`engine.py`（解方程/搜索序/符号判定）、`kernels/`（29 型的核族实现）、`solve.py`（调度）、`decompose.py`（组合拆解）、`server.py` + `templates/`（Flask 站点）、`integrand.py`（参数→被积函数，验证用）
- `bench/` — golden 数据集（`data/golden.jsonl`，采自线上 API）、`verify.py`（恒等式数值+符号验证）、`parity.py`（与线上一致性评测）、`harvest.py`
- `docs/` — kernel-spec（数学规格）、api-spec（协议）、behavior-notes、decompose-notes、verify-notes
- `tools/` — 预计算/维护脚本

## 运行

```bash
uv sync          # 或 pip install -e .
python -m attention_calculator.server   # waitress 起在 8080
```

## 评测

```bash
.venv/bin/python bench/parity.py bench/data/golden.jsonl
.venv/bin/python bench/verify.py bench/data/golden.jsonl
```
