# 复现路线图与分工

目标产物：与 zhuyidao.net 行为一致的完整实现（solver + Web + 组合拆解），全部结论由 `bench/data/golden.jsonl` 上的 parity/verify 数字背书。

## 阶段

1. **benchmark 先行**（进行中）：`bench/harvest.py` 在 devbox 采集 29 型 × 多界 golden 数据集 + combo 集；`bench/verify.py` 数值+符号双重验证每条恒等式；`bench/parity.py` 跑本方 solver 对 golden 的一致率。
2. **核族实现**（进行中）：8 个 family 模块按 kernel-spec.md 实现 `prove`/`render_equation`，矩用闭式/递推（文章印的展开式做 ground truth 测试）。
3. **装配**：solve.py 调度 → server.py 端点 → 前端复刻 → decompose.py。
4. **评测收口**：全量 parity + verify；不一致的 case 逐个归因（搜索序/上限/模板差异）迭代。

## Agent 分工（agentId 尾号 → 任务 → 产出文件）

- harvest → bench/{cases,harvest}.py, bench/data/*
- verify → integrand.py, bench/verify.py, docs/verify-notes.md
- quadlog → kernels/quadlog.py（pi, pi_n, arctan/arccot, catalan, zeta3）
- exp → kernels/{exp_family,hyperbolic}.py（e, e_q, e_pi, sinh/cosh/tanh/coth）
- log → kernels/log_family.py（ln_q, ln_q_square, artanh, arcoth）+ docs/log-notes.md
- trig → kernels/{trig_q,trig_pi}.py（sin/cos/tan/cot 弧度 + π倍数/角度）
- special → kernels/{beta_family,gamma}.py（phi, varpi, gauss, gamma）
- decompose → decompose.py + docs/decompose-notes.md
- webapp → server.py, render.py, templates/, static/
- behavior → docs/behavior-notes.md（线上行为规格）, bench/probe_api.py

## 已确认的线上行为（我亲手探测）

- 响应仅参数化系数 {m,n,a,b,c,au,bu,cu,u}；方程由 /get_integral_image 二次渲染。
- pi 用偶次幂族 x^{2m}(1-x²)^n·(a+bx²)/(1+x²)（非文章 type1 的三系数全幂族）。
- ln_q 分母升幂 (1+(q-1)x)^s，s 随情形变（ln2→2、ln3→3）。
- 错误："要证明的式子不等号方向反了"、"在指数不超过10的范围内未找到<方向的解"。
- e^π 走 [0,π] 上 sin 基 e^x 核；ϖ/G 走 4k+r 阶的 1/√(1-x⁴) 核；γ 是两段式复合积分。

## 风险

- 各类型 (m,n) 参数编码不一致（pi 的 m 对应 x^{2m}，varpi 的 m 对应 x^{4m+3}）——靠探测钉死。
- 组合拆解的界分配策略未公开，只能由响应数据反推。
- γ/ln A 等复合型依赖 ln q 子界，族间有接口依赖（log_family.ln_bound_proof）。
- 服务端原作者用 Mathematica 预计算+查表；我们用闭式矩现算，速度预计仍毫秒级，parity 数字会说话。
