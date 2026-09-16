# si_q / cin_q 落地笔记（2026-09-16）

W3 实现笔记。调研底稿 `docs/2026-09-16-w3-research-sicin.md` 的全部结构断言经实现复核；本文件记录实现选型、与调研稿的偏差、以及 kernel-spec 草稿（由 leader 并入 `docs/kernel-spec.md`）。

## 文件与类型

- `src/attention_calculator/kernels/sicin.py` — 一族两型：`si_q`（Si(q)=∫₀^q sin t/t dt）与 `cin_q`（Cin(q)=∫₀^q (1−cos t)/t dt = γ+ln q−Ci(q)）。同一台矩机器，按目标常数选核（与 varpi/gauss 按方向选核同构）。
- `src/attention_calculator/exact_check/sicin.py` — `check(kind, power, comp, bound, params)`。
- `tests/test_exact_sicin.py` — 63 通过 + 2 skip（注册后转绿）。

## 域与折叠

- `q ∈ ℚ \ {0}`；q=0 拒收（T/U 递推带 1/q，且 Si(0)=Cin(0)=0 是有理数，不属常数型）。
- 奇偶折叠：Si 奇、Cin 偶。统一写成 σ（si 取 sign(q)，cin 取 +1）：命题 `s·(C(q)−r)>0` 在 q̄=|q| 上解目标 `{C_q: sσ, 1: −sr}`、核方向取 kdir=sσ。σ·C(|q|)=C(q) 使认证恒等式即字面命题——负 q 的 si '<' 自动落到 '>' 阶梯核（Si(q)<r ⟺ Si(q̄)>−r），不需要特判。
- 渲染端同样取 |q| 进核（`render_equation` 里 `sp.sin(|q|x)`），LHS 常数项按原值回显。

## 参数编码

在 ln_q_cube 的九槽扩 `d_val`/`du_val` 之上再加 **`t_val`**（阶梯级别，非负整数）：`{m, n, a_val…d_val, au_val…du_val, u_val, t_val, unified_form}`。`u_val` 为四个系数分母 lcm。solution 串沿用 `a = X, b = Y, c= Z, d= W`。checker 对 `t_val < 0` 直接 ValueError——t<0 会把多项式尾丢光、退成裸核 sin(qx)/x，q>π 时定号性失效，属证书可构造的造假面，必须在复核侧拦截（详见下「安全边界」）。

## 符号与矩

- si 侧 span：`{si_q, sin_q, cos_q, 1}`；cin 侧：`{cin_q, sin_q, cos_q, 1}`。目标向量 `{C_q: kdir, "1": −s·bound}`（滤零）。
- 矩表：`t_moments` = `[{si_q:1}] + trig_q.sin_moments(kmax−1, q)`（T_k=S_{k−1} 平移复用）；`u_moments` 由 T 表生成 U_l=1/l−C_l、C_l=sin q/q−(l−1)/q·T_{l−1}（l=1 时 (l−1) 因子吞掉 T_0 的 si_q 符号）。
- 组合矩 `basis_moment(kind, kdir, t, table, q, m, n, j) = kdir·(special − trunc)`：special 为 Σ_{l,i}C(m,l)(−1)^iC(n,i)·table[l+i+j]，trunc 为纯有理 β(m,n,s)=∫x^s(1+x)^m(1−x)^n 加权和。`M = kdir·(special−trunc)` 一条式子统一四个方向——'>' 是「函数−下截断」、'<' 是「上截断−函数」，只差整体符号。
- β 上界：si 侧用到 s≤j+2R（R=2t+1），cin 侧 s≤j+2R−1（R=2t+2），r 循环界见 `trunc_moment`。
- 表深 kmax = 2·LIMIT+3（m+n+j ≤ 10+10+3）。

## 搜索

`mn_order(10)` 外层 (m,n)、内层 t∈{0,1,2}——t 是廉价轴（同一 (m,n) 下 4×4 求解成本与 t 无关），t-inner 让紧界在浅 (m,n) 处用高 t 命中，避免把 (m,n) 预算全烧在 t=0 的弱边上。命中即 emit；耗尽 NoSolution。

**结构观察：WrongDirection 在本族不可达。** C_q 行恒为 ±δ_{j0}（'>' 取 +、'<' 取 −），目标分量 kdir=±1 令解恒满足 a₀=+1>0——解出的 P 不可能恒非正，假命题只能 NoSolution 耗尽。exact 路径方向由 certified_cmp 前置判定，语义一致。

## 三次非负判据

复用 `kernels.ln_pow.cubic_nonneg`（QQ(√D) 临界点精确判号，qsign 原语），未写 Sturm——sympy.polys 的 Sturm 在 QQ 上虽精确但引入新机制，ln_pow 的 A∓B₀√D 临界值闭式已被 ln_q_cube 测试覆盖。**候选：`cubic_nonneg`/`qsign` 宜提升到 engine.py**（当前两处 import 自 kernels.ln_pow：kernels/sicin.py 与 exact_check/sicin.py；引擎化后 ln_pow 与 sicin 共引）。

## 安全边界（checker 侧）

- `t_val` 必须 ≥0——负值会静默退到裸核（见上「参数编码」）。
- 核定号性是定理不是数值检查：交错截断界对 u>0 严格成立（1−cos u≥0 逐次积分），t≥0、q̄>0 时 K 在 (0,1] 严格正、唯一零点 x=0。checker 不重复数值扫核符号，与 gauss_erf 的 e^{·}>0 同等待遇。
- `bool(integrand)` 非恒零守卫已并入 nonneg。
- m,n<0 的伪造参数自然死亡（空基→恒等式失败）；q=0 的伪造证书在 1/q 上 ZeroDivisionError→verify_cert False。

## 覆盖实测（实现侧复核）

- q=1 两型 6 位 floor/ceil 全过：t=0 在 (0,0) 覆盖 ~1e-4 级界，1e-6 界由 t=1/2 在同一 (0,0) 命中——(0,0) 基极强（四维行列式非退化，调研稿已证 1/q⁴）。
- q>π：si '>' 7/2 在 (0,0,t=2) 命中 183/100；cin '>' 10 命中 28/10（m=1,n=6,t=2）、'<' 命中 3（m=0,n=10,t=2）。cin '<' 3 需要 t≥1——t=0 边 3.104 不够锐，与调研表一致。
- si 20 '>' 真命题（Si(20)=1.548>1）预算内 NoSolution——深振荡区是本族诚实短板，'>' 需 t~q/2π 量级才压得住尾项，'<' 侧大 q 宜用 t=0+高 (m,n)（si 20 '<' 9 命中）。
- 与调研稿的唯一行为差：Si(7/2)<19/10 的样例证明在 (1,2,t=0)，本实现因 t-inner 搜索序在更早的 (0,0,t=1) 命中——同为有效证明，非复现错误。

## kernel-spec.md 草稿（「exact-only 型清单」追加）

```markdown
- `si_q` / `cin_q`（kernels/sicin.py）：`Si(q) ⋚ bound`、`Cin(q) ⋚ bound`，q≠0 由 power 槽携带；Si 奇（σ=sign q）、Cin 偶（σ=+1），负 q 折进目标符号与核方向 kdir=s·σ（矩机只见 |q|）。四维 span{C_q, sin_q, cos_q, 1}——C_q 只在基项 l=i=j=0 出现（F(0) 钉死），故 (i) P 必须三次 a+bx+cx²+dx³，参数扩 `d_val`/`du_val` 并新增 `t_val`；(ii) '<' 方向裸核下原理不可证，必须用互补核。核为 Taylor 余量阶梯 t=0,1,2：si '>' `(sin qx−T_{4t+3}(qx))/x`、si '<' `(T_{4t+1}(qx)−sin qx)/x`、cin '>' `((1−cos qx)−S_{2t+2}(qx))/x`、cin '<' `(S_{2t+1}(qx)−(1−cos qx))/x`，皆 (0,1] 严格正（交错截断引理）。基 (1+x)^m(1−x)^n——x^m(1−x)^n 在 m≥1 使 C_q 行全零恒奇异。矩：T_0=Si(q)、T_k=S_{k−1}（trig_q 平移）；U_0=Cin(q)、U_l=1/l−C_l、C_l=sin q/q−(l−1)/q·T_{l−1}；阶梯尾项只添纯有理 β(m,n,s)。P(0)=+1 被目标方程强制 → WrongDirection 不可达，假命题统一 NoSolution。非负判据复用 ln_pow.cubic_nonneg（QQ(√D) 精确三次规则）。LIMIT=10、T_LEVELS={0,1,2}。覆盖短板：q 深振荡区 '>' 需 t~q/2π，q=20 级预算内诚实无解。推导与覆盖表见 docs/2026-09-16-w3-research-sicin.md。
```

## 注册所需改动（leader 侧，本 session 未动）

- `kernels/__init__.py`：`EXACT_TYPES += ["si_q", "cin_q"]`。
- `solve.py`：`FAMILY += {"si_q": "sicin", "cin_q": "sicin"}`。
- `integrand.py` `constant_mpf`：`"si_q": lambda: mp.si(q)`、`"cin_q": lambda: mp.euler + mp.log(abs(q)) - mp.ci(abs(q))`（Cin 偶函数，certified_cmp 对负 q 也走 |q|）。
- 可选：`certificate._SYMBOL_TEX` 加 `"si_q": "\\mathrm{Si}", "cin_q": "\\mathrm{Cin}"`（缺省回退 \mathrm{si_q} 也能印）。
- 无需 EXPONENT_LIMIT 条目（默认 10 即 LIMIT）。
