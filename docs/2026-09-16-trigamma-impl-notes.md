# psi1_q 实现笔记（2026-09-16）

W3 exact-only 型 `psi1_q`：`ψ′(q) ⋚ r`，q ∈ ℚ_{>0}，power 槽载 q（系数恒 1）。依据 `docs/2026-09-16-w3-research-trigamma.md` 的望远镜核对设计落地：`K_> = x^{q−1}(x−1−lnx)/(1−x)`、`K_< = x^{q+N−2}(1−x+x·lnx)/(1−x)`，N = 最小非负整数使 q+N>1（即 q>1 → 0，0<q≤1 → 1）。本文档记录实现侧与该调研的差异与编码约定；kernel-spec 合并草案在末节。

## 与调研的三处实现差异

1. **常数尾不需要**。调研按"a,b ≥ 0"推导出 `>` 侧窗口 [1/q, ψ′)、`<` 侧 (ψ′, U_0(0)]，松散端要常数尾兜底。实际站端非负判据是 `a≥0 且 a+b≥0`（engine.poly_nonneg），而解出的 P 恒有 `a+b=1`（ψ′ 行系数 ±1 配 target ±1）——判据退化为 `a≥0` 单边条件：`>` 在首个 `r ≤ L_{m+1}` 的 m 命中（覆盖全部 r<ψ′，含 r<1/q 与负界，b<0 合法），`<` 在首个 `r ≥ U_{m+1}` 命中（覆盖全部 r>ψ′，含 r>U_0(0)）。松散侧统一在 m=0 出证（P 形如大 a 配负 b），两核对全部真命题完备，无尾编码字段。
2. **N 不进参数**。N 由 q 决定（`shift_n`），且 `U_u(N+1) = U_{u+1}(N)`——更大的 N 只是重标 m、可达 U 集更窄，canonical N 严格最优；参数只有标准九槽，checker 从 q 重算 N。
3. **一维 m 扫描**。n≥1 时 (1−x)^n 消掉核极点、ψ′ 行归零（Σ(−1)^t C(n,t)=δ_{n0}），二维搜索结构退化；`solve` 直接以 `(m, 0)` 计划喂 `engine.search`，发射参数 n 恒 0。假命题两侧都既解不出非负也解不出非正 P（a+b=1>0 堵死 nonpos）——扫描耗尽报 NoSolution（核内诚实回答；走 solve.prove(exact=True) 时 certified_cmp 早已判反）。

## 编码与符号

- 矩符号键即类型名 `"psi1_q"`；矩公式 `>`：`{psi1_q: +1, "1": −L_u}`、`<`：`{psi1_q: −1, "1": U_u}`，`L_u = T_u(q)+1/(q+u)`、`U_u = T_{N+u}(q)+1/(q+N+u−1)`、`T_k = Σ_{i<k}(q+i)^{−2}`（`tail_sums` 增量构造，`basis_moment` 仍是通用二项式组合）。
- 参数即 emit 标准九槽（P=a+bx 两系数，c_val/cu_val=0）；无核选择子字段——comp 即方向核选择（varpi/gauss 先例），无 N 字段。
- LIMIT=256（beta 族 EXACT_LT_LIMIT 同档）：一维扫描单档是 2×2 有理解，满扫 ~10ms；q~1 紧界深度 ~1/(2m²)，覆盖到 δ~1e-5。需要 EXPONENT_LIMIT["psi1_q"]=256 让报错文案对齐。
- 域：`q>0` 全域。`q≤0` 报 ValueError（0,−1,−2,… 是 ψ′ 极点；负非整数是调研标注的可选扩展，首版不做）。certified_cmp 对 q≤0 求值 mp.psi 要么极点异常（→None→核内 ValueError）要么负非整数有限值（符号判定后仍核内 ValueError）。
- 收敛：`ψ′−L_m ≈ 1/(2(q+m)²)`、`U_m−ψ′` 同阶二次——1e-3 界 m~20+，1e-4 界 m~70~140，1e-6 界 m~700（超预算诚实 NoSolution）。

## kernel-spec.md 草案节（合并用，建议加在 exact-only 型清单）

- `psi1_q`（kernels/trigamma.py）：`ψ′(q) ⋚ bound`，q>0 有理数由 power 槽携带，系数恒 1。方向核 `>`：`x^{q−1}(x−1−lnx)/(1−x)`；`<`：`x^{q+N−2}(1−x+xlnx)/(1−x)`，N=0（q>1）或 1（0<q≤1），两核定号均为 ln t≤t−1 的改写。矩 `∫x^uK_> = ψ′(q)−L_u`、`∫x^uK_< = U_u−ψ′(q)`（L_u=T_u+1/(q+u)、U_u=T_{N+u}+1/(q+N+u−1)、T_k=Σ_{i<k}(q+i)^{−2}）落 span{ψ′(q),1}；(1−x)^n 消 ψ′ 行 → 一维 m 扫描，P=a+bx 恒 a+b=1，站端判据下全真命题可证无常数尾。EXPONENT_LIMIT=256。
- constant_mpf：`"psi1_q": lambda: mp.psi(1, q)`（power 即参数 q，**不乘**系数——与 sin_q 同，勿照 gamma/varpi 乘 q）。
- direction_f：无需 float 预检项（exact-only 型不进站端方向预检路径）。
- 渲染：程序化 render_equation，LHS `\psi_1\left(q\right)`；被积函数印合并幂 `x^{m+q−1}`/`x^{m+q+N−2}`。

## 测试与自检

`tests/test_exact_trigamma.py` 65 过 2 跳（跳的是待注册的全管线用例）。要点：tight floor/ceil 双向（9 个 q × 1e-3、6 组 1e-4）、文档已验恒等式逐条复现（ψ′(1)>8/5 落 m=2 P=(2+3x)/5、ψ′(1)<7/4 落 m=0 P=x、ψ′(1/3)<21/2 落 m=0 a=1/9 b=8/9 与调研一致）、松散侧无常数尾覆盖（r<1/q、r>U_0、负界、bound=0）、假命题核内 NoSolution、q≤0 ValueError、n≥1 伪造参数恒等式诚实失败、零被积函数守卫。

数值自检：矩公式 `∫x^uK_> = ψ′(q+u)−1/(q+u)`、`∫x^uK_< = 1/(q+N+u−1)−ψ′(q+N+u)` 在 6 个 q × u=0..3 以 mp.quad 复核（端点奇异档 u=0 为 quad 噪声，u≥1 全干净；调研已 100dps 全量核验）。
