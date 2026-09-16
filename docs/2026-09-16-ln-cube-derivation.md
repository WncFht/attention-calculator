# ln_q_cube 推导笔记（2026-09-16）

新类型 `ln_q_cube`（exact 模式专属，站点无此型）：证 `(ln q)³ ⋚ bound`，q 为有理数且 q > 1，power 槽携带 q。

## 核与矩空间

沿用 ln 族构造：基 `x^m(1-x)^n`、域 [0,1]、分母幂 `s = max(m, n, 1)`，核取 `ln²(1+cx)/(1+cx)^s`（c = q−1）。与 ln_q（ln⁰ 核出 ln¹）、ln_q_square（ln¹ 核出 ln²）同构——核的对数幂次恒比目标低一，目标常数由 `u⁻¹` 项产出。

代入 u = 1+cx：`x^k = c^{-k} Σ_r C(k,r)(−1)^{k−r} u^r`、`dx = du/c`，矩化为 `∫₁^q u^p ln²u du` 的和（p = r−s）。闭式：

- p = −1：`J_{-1} = ln³q / 3`——目标常数由此进入矩空间；
- p ≠ −1：对 `u^p ln²u` 逐次 IBP 得 `J_p = q^{p+1}(ln²q/(p+1) − 2 ln q/(p+1)² + 2/(p+1)³) − 2/(p+1)³`。

故矩空间恰为 4 维 span{ln³q, ln²q, ln q, 1}（s = max(m,n,1) 下 s−1 ≤ m+n−1 恒小于最高矩次数 m+n+3，p=−1 项总在、ln³ 分量总在），需要三次 P = a+bx+cx²+dx³ 共 4 个未知量。备选方案核查：ln³ 核在同样 s 规则下 `u^{-1}ln³u → ln⁴q/4`，矩空间变 5 维、P 升四次，更差；ln² 核取 s = 0 或 s ≥ m+n+5 则 p=−1 不出现、丢目标常数；3 系数 P 对 4 符号超定，generic 无解——deg-3 是构造内最小维数。

矩公式经 200dps 数值核验（`mpmath.quad` 对照，c ∈ {1, 1/2, 2, 1/3, 4, 3}、s ≤ 5、k ≤ 7 全部吻合至 ~1e-200）。

## 三次非负判据（精确，QQ(√D)）

`P(x) = a+bx+cx²+dx³` 在 [0,1] 的最小值只可能落在端点或 P′ = 3dx²+2cx+b 的实根处。d = 0 时回落 engine.poly_nonneg 的二次判据；d ≠ 0 时：

1. 端点：a < 0 或 a+b+c+d < 0 即否；
2. D = c²−3bd < 0：P′ 无实根、P 单调，端点非负即足；
3. D ≥ 0：每个根 `x*± = (−c ± √D)/(3d)`，先判 0 < x* < 1，再对落在内部的根验 P(x*) ≥ 0。

关键化简：用 3dx*² = −2cx*−b 把 P 在临界点降为一次式，得 `P(x*±) = A ∓ B₀√D`，其中 `A = a + (2c³−9bcd)/(27d²)`、`B₀ = 2D/(27d²) ≥ 0`。符号判定全部落在 QQ(√D)：元素写作 α+β√D，异号情形归约为 α² 与 β²D 的有理比较；D 为平方有理数（x* 本身有理）与 D = 0 都被同一比较自然覆盖。两个判断（根在 (0,1) 内、P(x*) 符号）共用同一个 qsign 原语。该判据是充要的，非充分条件；B₀√D 公式经 80dps 下 3000 组随机系数核验，判据整体与数值最小值在 50000 组随机三次上零分歧。

## 参数编码与注册

响应沿用九个站点槽位语义，第四个系数扩 `d_val`/`du_val` 两字段：`u_val = lcm(a,b,c,d 分母)`，`*u_val` 为整分子，`*_val` 为约分文本；`solution` 串按族内怪癖续排 `a = X, b = Y, c= Z, d= W`。注册侧（merge 时由 leader 接线）：`kernels.TYPES += "ln_q_cube"`、`solve.FAMILY["ln_q_cube"] = "ln_pow"`、`integrand.constant_mpf["ln_q_cube"] = mp.log(q)**3`、exact_check 按 FAMILY 自动落到 `exact_check/ln_pow.py`；域校验 q > 1 复用 ln 族文案。

## 实测覆盖（m,n ≤ 10 预算内）

q ∈ {2, 3, 3/2, 5, 10, 4/3, 5/4} 双方向均可证；ln³2 的 12 位小数紧界（gap ~1e-12）在 (7,10)/(9,10) 命中，15 位紧界诚实 NoSolution；假命题在扫描中由恒 ≤0 的解认证为 WrongDirection。
