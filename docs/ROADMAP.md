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

- **端到端验证**：pi<22/7 现场返回 m=3,n=3,a=47/120,b=-13/120 —— 与引擎原型（mn_order+solve_moment+poly_nonneg + J_k 递推）逐位一致，架构与搜索序已对齐。
- **power 语义**：e 型的 power 是**系数**（power=2 → 2e，渲染 `19/3 - 2e`）；e_q 的 power 是**指数**（power=2 → e²，核 e^{2x}）。
- **API 字段陷阱**：/calculate 的字段是 `type`；误发 `constant=e` 会被**静默忽略并缺省回退到 type=pi**（响应 `"type":"pi"`）。曾因此误判"e 型恒等式为假"——实为 π 型参数被 e 渲染器渲染；e 型实测完全正常（双向支持：`e<3`→(0,1)P=x、`e>8/3`→(1,1)P=x/3）。**parity/探测脚本必须发 `type=` 字段**。
- **搜索序 tie-break**：同 m+n、同 |m−n| 时 **m 小者优先**（pi>8/3→(0,1)；全部非对称 pi golden 均 m<n；400+ 条参数逐位一致）。
- **指数上限按单指数计**：m≤cap 且 n≤cap（pi 实测用到 (15,18)，m+n=33>30）。
- **ln 类分母幂**：`s=max(m,n,1)`（ln_q 61/61 + ln_q_square 31/31 数值拟合全中；s 不回传）。
- **错误码分层**：格式错误 400；值域/方向/无解 404；image 端点失败统一 500。错误响应体只有 `{"error"}` 无 `success` 键。
- **校验顺序（成对探针钉死）**：type（缺省→pi 回退；显式空串→400）→ comparison → **右侧整组**（格式 `^\d+(/\d+)?$` → 分母 0 → 分子或分母 ≥10^16）→ **左侧整组**（同三项，文案"左侧系数…"）→ 类型值域 404 → 求解。
- **/get_integral_image 只接受 GET+query**：对它 POST 表单一律 500（grid:*:img 探针曾误判为类型覆盖缺失；pi_n 等全类型 GET 均正常返回方程）。
- **方向假优先于无解（float64 判定）**：命题数值为假时报"方向反了"——即使搜索耗尽都没搜出恒≤0 的 P（实测 arctan_q 3：`>5/4`→方向反了，`<5/4`→未找到解；arctan3≈1.2490）。判定在 **float64** 精度进行：对 zeta3/gamma 发送 float 与常数相等但方向为假的界，站端放行进入搜索最终报"未找到解"（float 差严格 <0 才算反）。实现：solve.prove 捕获 NoSolution 后做 `float(C)−float(r)` 严格判号，假则重抛 WrongDirection；恒等式精确成立保证真命题永远搜不出恒≤0 的 P，此路径只触达假命题。
- **unified_form 恒 {}**：76 个成功响应全部如此——死字段。
- **服务端确定性**：同请求重复返回完全相同参数。
- **恒等式可信**：verify 50dps 全量 574/574 valid + 抽样精验——站点输出数学上全部成立，parity 参数可直接信任。
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
