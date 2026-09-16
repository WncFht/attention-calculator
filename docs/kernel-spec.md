# 注意力计算器复现 — 数学规格

来源：作者知乎专栏文章（本地归档 `obsidian/output/zhihu-mathematical/量化调酒师的数学文章/`）+ 对 zhuyidao.net 线上 API 的实测。

## 总框架

要证 `C ⋚ r`（C 为目标常数、r 为有理数），构造恒等式

```
±(C - r) = ∫_a^b f(x) dx,  f 在 [a,b] 上恒不变号且零点有限
```

统一被积函数形式：

```
f(x) = B_{m,n}(x) · P(x) · K(x)
```

- `B_{m,n}` 是族基函数（如 `x^m(1-x)^n`、`x^{2m}(1-x^2)^n`、`sin^m x (1-sin x)^n`），在域上恒非负；
- `P` 是带 2~3 个待定系数的小多项式 `a+bx` 或 `a+bx+cx^2`（偶核用 `a+bx^2`；lemniscate 两型用 `a+bx^4`；[0,π] 与 [0,π/2] 域用 `a+b·sin x`）；
- `K` 是核函数，使所有"矩" `∫ B·x^k·K dx` 落在 `span_Q{1, C, D_1, ...}` 有限维空间里。

求解：矩对 (a,b,c) 线性 → 在 Q 上解小线性方程组（engine.gauss_solve，要求唯一解）

- 目标矩 = `s·(coef·C − r)`，s=+1 为 `>`、−1 为 `<`：C 的系数 = s·coef，常数项 = −s·r，寄生常数 D_j 系数 = 0；
- coef 是请求的 power 系数——pi/catalan/zeta3/e/e_pi/golden/varpi/gauss 直接进目标；ln/trig/hyperbolic 各型为 1；artanh/arcoth 为 1/2（回报的 a,b 自带半因子）；pi_n 的 p/q 形把常数项换成 bound^q、coef 固定为 1。

然后检查 `P` 在域上不变号；变号则换下一组 (m,n)。

## 搜索顺序（作者亲述 + 实测修正）

1. `m+n` 升序；
2. 同 `m+n` 时 `|m-n|` 小的优先（即先取平衡的）；
3. 同 `|m-n|` 的镜像对按 **m+n 奇偶**定向：和为奇数 m 小者先、和为偶数 m 大者先——sum=2 的序是 (1,1),(2,0),(0,2)，sum=3 是 (1,2),(2,1),(0,3),(3,0)，等价于沿 n = ⌈s/2⌉, ⌈s/2⌉−1, ⌈s/2⌉+1, ⌈s/2⌉−2, … 走（engine.mn_order）。实测：pi>8/3 取 (0,1)（奇数和）；ln²(3/2)<17/100 取 (2,0) 而非同样可行的 (0,2)（偶数和，决定性证据）。旧版"镜像对恒 m 小者优先"只覆盖了奇数和情形——以实测为准更正，全部 golden 参数逐位吻合。
4. e、π 两类型指数上限 30；其余类型上限 10。**上限按单指数计**——m≤30 且 n≤30，非 m+n≤30（实测 pi 用到 (15,18)，m+n=33）。耗尽报错文案随类型上限插值："在指数不超过30的范围内…"（e/pi）/"…不超过10…"（其余）。

方向反判定分三层：

- **float64 预检**（仅核会崩的类型需要）：e_q（math.exp）、sinh/cosh/tanh/coth（math.*，coth 内部另有 1/tanh(0)）、arctan_q/arccot_q（除以 q）、pi_n（表外指数 KeyError）、gamma（coef·γ）、ln_q_square q∈{5,7}（站点自身崩溃）、trig_pi 四型（mpmath 50dps）。严格为假的断言在此直接 方向反了；真/等值继续——等值多由输入检查先报 `二者相等`（Niven 点）或随核内崩溃变 500。zeta3 '>' 与负界不进预检（zeta3 阈值在 float c 上方 1 ulp 外，留给下方兜底覆盖）。
- **扫描中**：`engine.search` 命中"被积函数恒≤0"的候选抛 WrongDirection 作**原语**，`solve.prove` 拦截后按方向分流：'<' 且预检已映射（c 非 None）→ 改报"未找到解"（站端 '<' 的方向判定只在预检发生，扫描里的非正 P 直接耗尽）；'>' 不判反，落入下方统一兜底（bound==c 的 cos/golden/sin_q 与 bound<C 的 sin_pi_q 实测均报"未找到>方向的解"）；未映射的型（artanh/arcoth/trig_pi）核内自带方向判定，WD 原样上报。**trig_pi 例外**：该四型恒≤0 解跳过继续找（engine.search defer），搜完仍无非负解才抛。
- **耗尽兜底**（solve.prove）：搜索耗尽（含 '>' 途中非正被吞的 WD）后统一用 float64 比较常数与界，命题为假 → 报"方向反了"（实测 arctan 3>5/4）；为真 → 维持"未找到解"。float 差恰为 0 的边界（zeta3/gamma 精度持平的假命题）按不 falsify 处理，与站点一致。详见 fidelity-notes.md。

## 非负性判据（作者亲述，只需判断 P）

- 一次 `a+bx`：`a≥0 且 a+b≥0`。
- 二次 `a+bx+cx^2`：`c≤0` 时同上（端点非负即可）；`c>0` 时若对称轴 `-b/2c ∈ (0,1)` 需 `4ac-b²≥0`，否则端点非负。
- 偶核 `a+bx^2`：等价于 `a+b·t` 在 `t∈[0,1]` 非负 → `a≥0 且 a+b≥0`。
- `a+b·sin x`（域 [0,π] 或 [0,π/2]，sin x∈[0,1]）：同一次型判据。

## 类型表（网站 29 种 → 核与矩空间）

记 `M(...)` 为"矩 ∈ span{...}"。`s` 为分母幂（ln 类可升幂加速收敛，见技巧6）。

**ln 类分母幂实测规则**：`s = max(m, n, 1)`（log 族全部参数记录逐一数值拟合全中，见 log-notes.md；s 不在响应参数里，只体现在渲染方程与矩构造中）。

| type | 域 | 核 K(x) | 基 B_{m,n} | P | 矩空间 |
|---|---|---|---|---|---|
| pi | [0,1] | 1/(1+x²) | x^{2m}(1-x²)^n | a+bx² | {π, 1} |
| e | [0,1] | e^x | x^m(1-x)^n | a+bx | {e, 1} |
| pi_n (π^k) | [0,1] | ln^{k-1}(1/x)/(1+x²) | x^{2m+δ}(1-x²)^n，δ=1 当 k 偶 / 0 当 k 奇 | a+bx² | {π^k, 1} |
| e_q (e^q) | [0,1] | e^{qx} | x^m(1-x)^n | a+bx | {e^q, 1} |
| ln_q (q>1) | [0,1] | 1/(1+(q-1)x)^s | x^m(1-x)^n | a+bx | {ln q, 1} |
| ln_q_square ((ln q)²) | [0,1] | ln(1+(q-1)x)/(1+(q-1)x)^s | x^m(1-x)^n | a+bx+cx² | {ln²q, ln q, 1} |
| sin_q | [0,1] | sin(qx) | x^m(1-x)^n | a+bx+cx² | {sin q, cos q, 1} |
| cos_q | [0,1] | sin(qx) | x^m(1-x)^n | a+bx+cx² | {sin q, cos q, 1} |
| tan_q | [0,1] | sin(qx)，结果整体除 cos q | x^m(1-x)^n | a+bx+cx² | 解 sin q − p·cos q |
| cot_q | 同 tan | 同 tan（先取倒数） | | | |
| sin_q_degree | [0,π/2] | sin((1-2q)x)，q=度数/180 | sin^m x(1-sin x)^n | a+b sin x | {sin(πq), 1} |
| cos_q_degree | [0,π/2] | sin(2qx) | sin^m x(1-sin x)^n | a+b sin x | {cos(πq), 1} |
| sin_pi_q | [0,π/2] | sin((1-2q)x)，0<q<1/2 | sin^m x(1-sin x)^n | a+b sin x | {sin(πq), 1} |
| cos_pi_q | [0,π/2] | sin(2qx) | sin^m x(1-sin x)^n | a+b sin x | {cos(πq), 1} |
| arctan_q | [0,1] | 1/(1+(qx)²) | x^{2m}(1-x²)^n | a+bx² | {arctan q, 1} |
| arccot_q | — | 归约为 arctan(1/q) | | | |
| sinh_q | [0,1] | sinh(qx) | x^m(1-x)^n | a+bx+cx² | {sinh q, cosh q, 1} |
| cosh_q | [0,1] | sinh(qx) | 同上 | 同上 | 同上 |
| tanh_q | [0,1] | sinh(qx)，整体除 cosh q | 同上 | 同上 | 解 sinh q − p·cosh q |
| coth_q | 同 tanh（先倒数） | | | | |
| artanh_q | — | 归约 ln((1+q)/(1-q))/2 | | | |
| arcoth_q | — | 归约 ln((q+1)/(q-1))/2 | | | |
| gamma | [0,1] | 两段式：x^n·((2-x)/2 − 1/(1-x) − 1/ln x)（上界）或 x^n·(1/ln x + 1/(1-x) − 1/2)（下界），再叠加一个 ln(n+1) 的子证明积分 | — | — | {γ, ln(n+1), 1} |
| golden (φ) | [0,1] | √(4+x) | x^m(1-x)^n | a+bx | {√5, 1}（φ=(1+√5)/2） |
| catalan | [0,1] | ln(1/x)/(1+x²) | x^{2m}(1-x²)^n | a+bx² | {C, 1} |
| zeta3 | [0,1] | ln²x/(1+x²) | x^{2m+1}(1-x²)^n（奇次幂） | a+bx² | {ζ(3), 1} |
| e_pi | [0,π] | e^x | sin^m x(1-sin x)^n | a+b sin x | {e^π, 1} |
| varpi (ϖ) | [0,1] | 上界：x^{4m+3}(1-x)/√(1-x⁴)；下界：x^{4m+1}(1-x)/(π√(1-x⁴)) | 固定形 | a+bx⁴ | {ϖ,1} 或 {ϖ⁻¹,1} |
| gauss (G) | [0,1] | 下界：x^{4m}(1-x)/(π√(1-x⁴))；上界：x^{4m+2}(1-x)/√(1-x⁴) | 固定形 | a+bx⁴ | {G,1} 或 {G⁻¹,1} |

## gamma（两段式复合核：主核 + ln 子证明）

恒等式 = 主核积分 + ln(N+1) 的子证明积分，两段各分至少一半缺口：

- **N=0** 当且仅当剩余已是非负常数：'<' 时 r ≥ s₀ = 3/4（常数尾 r−3/4）、'>' 时 r ≤ r₀ = 1/2。参数形 m=n=au=bu=cu=u=0，a_val = 该常数文本，c_val 是方向标志（'<'→"2"、'>'→"1"，非多项式系数）。
- 否则取**最小 N ≥ 1** 使子积分 ≥ 主积分：'<' 判 `s_N − ln(N+1) ≤ (r+γ)/2`，'>' 判 `r_N − ln(N+1) ≥ (r+γ)/2`，全部 **float64** 比较；s_N = H_N + 1/(N+1) − 1/(2(N+2))，r_N = H_N + 1/(2(N+1))。实现侧 N 上限取 600（保证 H_N 分母可算；站点真实上限未探明，golden 内最大 N=17）。
- 子证明委托 ln_q 核（q = N+1）：'<' 需 `ln(N+1) > s_N − r`，'>' 需 `ln(N+1) < r_N − r`。回报的 m, n, au_val, bu_val, u_val 全是子证明参数；N 塞在 cu_val，c_val="0"。
- 边界：power=0 先过 float64 预检（假命题 方向反了，真/等值进 `bound/|power|` 除零 → 500）；power<0 交换方向后证 `|power|·γ ⋚ bound`。

渲染怪癖（/get_integral_image）：

- 等号空格随方向：'<' 排 `=\int_0^1`（无空格），'>' 排 `= \int_0^1`。
- 主项 `coef·x^N(kernel)`：coef 为整数时级联前写（`2 x^{11}`）、coef=1 消失；coef 为真分数时**走 sympy 折叠** `\frac{x^{N}}{d}`（站点 1/2 实测，非 `\dfrac{1}{2}x^{N}`）；N=0 时只剩 coef 文本。
- 子分数 `\frac{t·x^m(1-x)^n(au+bu·x)}{u'(Nx+1)^{s'}}`，s' = max(m,n,1)；**t/u' = coef/u 约分后的分子/分母**（coef=2、u=936 → 站点印 1/468；coef=3、u=3 → 4/1 印在分子）。s'=1 时分母展开 `u'N·x + u'`。
- u=0 两退化：au=bu=0 → 方括号内 `主项 + rat_tex(a_val·coef)`（常数尾为 0 时连方括号都省，只剩主项）；au/bu ≠ 0 → sympy zoo×分子（站端显示 ~∞ 倍因子）。
- /get_integral_image 对 gamma 的 u_val 下限放宽到 0（其余类型 ≥1）。

## varpi / gauss（lemniscate 核，含兜底模板）

- n 是方向标志（'>'→0、'<'→1）而非指数；指数在 x^{4m+res} 的 res（varpi '>'=1、'<'=3；gauss '>'=0、'<'=2）。'>' 方向 m 搜 0..10；**'<' 方向的正确解只搜 m ≤ LT_M_LIMIT**（gauss 4、varpi 10）。
- '>' 耗尽 → 兜底模板 `∫ t·poly·√(1-x⁴)/π dx + b`（poly = (1−x) 对 gauss、x(1−x) 对 varpi；矩 {G:1/3,1:−1/8} 与 {ϖ⁻¹:−1/5,1:1/8}，t 由目标常数系数定，要求 t≥0 且 b>0）。参数形 m=n=0、a_val="0"、b_val=b、cu_val="1"、au_val/u_val = 约分后的 t。
- '<' 耗尽 → 先同形模板（cu_val="2"；t·x^res(1−x)√(1-x⁴) + b，res=2 gauss/3 varpi；矩 {G⁻¹:1/5,1:−1/6} 与 {ϖ:−1/21,1:1/6}；要求 t≥0 且 **b>0 严格**——b=0 边界继续落到转置解）。
- 再不行 → **转置坏解**：m = LT_M_LIMIT+1 单发，把两条基矩按 [符号, 1] 行当方程组解 rhs=target——恒等式本身不成立，站点只靠 a+bx⁴ 非负检查拦着。前提：命题数值为真（lhs_mpf>0）且解非负。gauss 转置解恒正 → (G, ~0.855] 窗口全落这里；varpi 恒 a<0 → 永不触发。
- 渲染：cu_val 非零即兜底模板——`(au/u)·x^res(1-x)·√(1-x⁴)[/π] + b_val`，'>' 带 /π 分母、'<' 不带；gauss 的乘积括进 `\left(\right)`、varpi 裸排；`\int_0^1` 后空格数按 kind×comp 实测固定（gauss '>' 双空格、'<' 单空格；varpi 相反）。常规形里基底的 (1−x) 因子排在分母分数内：`\dfrac{(1-x)}{\pi\sqrt{1-x^4}}`（'>'）/ `\dfrac{(1-x)}{\sqrt{1-x^4}}`（'<'）；e=0 的 P 块回落 sympy 直排（会把 (au+bu·x⁴)/u 拆成两个分式）。

## trig_pi 四型（sin_pi_q / cos_pi_q / sin_q_degree / cos_q_degree）的站点实测行为

站点对这四型不是现算积分，而是代入按 (m,n) 预存的闭式（作者专栏文章载明 Mathematica 预计算管线）。674 条 golden 全部吻合以下模型：

1. **数值方向预检**。先比较 r 与常数 $C=\cos(\pi\alpha/2)$（sin 型即 $\sin(\pi q)$、cos 型即 $\cos(\pi q)$）的数值大小；方向错直接 `要证明的式子不等号方向反了`，根本不进 (m,n) 扫描。全量 golden 里 "方向正确但无证明"的记录为零——方向对必出解。等值情形（Niven 三点： sin π/6、cos π/3、sin 30°=1/2）仍由输入检查报 `二者相等`。
2. **扫描跳过非正解**。恒≤0 的 (m,n) 候选不报错、继续往后找—— `sin_pi_q 1/10 > 1292/4181` 跳过 (1,8) 的非正解后才在 (5,9) 命中； 搜完若无恒非负解，见非正则报方向反、否则报未找到。
3. **(m,n)=(1,8) 的 j=0 存储公式是坏的**。站点该档用的矩是 $M_0' = M_0 + \delta(\alpha)\cdot(C-1)$，其中

   $$\delta(\alpha) = \frac{(\alpha-1)(\alpha-2)\,Q_{16}(\alpha)}{256\,\alpha\prod_{k=1}^{9}(\alpha^2-k^2)}$$

   $Q_{16}$ 为 16 次整系数多项式（系数见 `kernels/trig_pi.py` 的 `BIAS18_NUM`）。分母恰是真 $M_0(1,8)$ 的自然分母——即存储公式 = 真公式 + $N/D\cdot(C-1)$，像是预计算结果在 a-系数多项式上叠了 $+N/D$、常数项叠 $-N/D$ 的转录/化简错误。该式由 33 个 $\alpha$ 上有理插值恢复（19 点拟合 + 14 点留出全部逐位命中）。
4. **后果**：j=0 解出的 (a,b) 与真值不同。界 $r$ 落在真 (1,8) 解不定号、 坏 (1,8) 解非负的窗口时，站点在 (1,8) 报出坏系数——即 golden 里 30 条 "假恒等式"（方程两边实际不等，差 $a\,\delta(\alpha)(1-C)$）。 若坏解恒≤0 则被跳过（见第 2 点）。故"站点 (m,n,a,b,c,u) 对应什么真实 被积函数"的答案是：**没有**——参数对应的是损坏的存储公式， 渲染出的积分等式本身不成立。

渲染侧：被积函数作为一个 sympy 乘积整体 latex，因子顺序为 数值系数 → `(1-sin x)^n` → `(a+b sin x)` → `sin(αx)` → `sin^m x`； ` \cdot ` 只插在前一字符为数字或 `}`、后因子为 `\left(` 紧跟数字处 （负号开头的 `\left(- …` 不加；裸 `\left(1-\sin x\right)` 作前因子时 其后仍是空格——它以 `)` 结尾而非 `}`/数字）。

实测样例（线上返回）：

- `pi < 22/7` → `22/7 - π = ∫₀¹ x⁶(1-x²)³(47-13x²)/(120(1+x²)) dx > 0`（m=3,n=3，P=47-13x²，u=120）
- `e > 8/3` → `e - 8/3 = ∫₀¹ x²(1-x)eˣ/3 dx > 0`（m=1,n=1，P=x/3）
- `ln 3 < 11/10` → `∫₀¹ x²(1-x)³(40x+6)/(45(2x+1)³) dx`（s=3）
- `e^π > 23` → `∫₀^π (1-sin x)⁵(547950-422240 sin x)eˣ sin²x/336633 dx`
- `ϖ < 8/3` → `∫₀¹ x²⁷(3319635x⁴+2432199)(1-x)/(353600√(1-x⁴)) dx`
- `G > 4/5` → `∫₀¹ x⁸(44x⁴+6)(1-x)/(5π√(1-x⁴)) dx` —— 即文章示例
- `sin(π/5) > 1/2` → `∫₀^{π/2} (1-sin x)(637 sin x+125) sin(3x/5)/750 dx`（1-2q=3/5）
- `γ < 3/5` → `∫₀¹ [x⁵((2-x)/2-1/(1-x)-1/ln x) + x²(1-x)⁴(18460x+2811)/(168(5x+1)⁴)] dx`（前半给 γ+ln6 项、后半是 ln6 上界子证明）

## mode=exact（数学正确性路径，2026-09-16 起）

`solve.prove(..., exact=True)` / `POST /calculate mode=exact` 走正确性路径，与站点行为解耦（site 路径冻结，见 `docs/2026-09-16-math-correctness-plan.md`）：

- **方向判定**：`solve.certified_cmp` 递增精度（80→2400 dps）数值+护栏带认证 `sign(C−r)`；假命题预检即 `方向反了`、等值 `二者相等`、常数不可实值求值（域外输入）→ None 交给核内域校验。扫描中命中非正 P（真矩下即对偶不等式认证）→ 直接 WrongDirection 上报，不做 '<'→未找到 的站点映射。耗尽 → NoSolution（真命题预算内证不出的诚实回答）。
- **核内真系统**：`module.prove(..., exact=True)` 要求核用真矩求解——trig_pi 走 `basis_moment`（非 site_bias）、beta '<' 不走转置档、ln_q_square q∈{5,7} 正常求解、power≠1 令发射参数满足印刷 LHS（∫印刷被积函数 == `power·C − bound` 的符号向量）。
- **发射自检**：`exact_check.verify(kind, power, comp, bound, params)` 对每条发射证明做 ℚ 字典相等复核（`combine(coeffs, true_basis) == target` 且 `poly_nonneg`）；失败是自家 bug → InternalError。exact_check 与模式无关：对 site 模式参数跑它即可复核站点输出的真伪（verify.py 的 ℚ 精确版）。
- **exact-only 类型**：`kernels.EXACT_TYPES` 注册无站端对应的新类型（当前：zeta5、zeta7）。/calculate 仅在 mode=exact 下接受（site 模式仍按站端口径 400）；solve.prove 对不带 exact 的调用报 ValueError；certified_cmp 常数走 `integrand.constant_mpf` 表。/get_integral_image 仍只认站端 29 型——exact 类型的渲染经 kernel 的 render_equation 程序化调用。

## 参数语义（a_val/b_val/c_val/au_val/bu_val/cu_val/u_val）

常规模板：`u = lcm(a, b, c 的分母)`，`au_val = a·u` 等整分子，`a_val` 为约分文本；`m, n` 为基指数（响应里 int），`au/bu/cu/u_val` 为字符串。逐型复用：

- trig_pi 四型：`c_val = α`（角度参数 1−2q 或 2q，非多项式系数）；degree 型回报的 type 归一为 sin_pi_q / cos_pi_q。
- artanh_q / arcoth_q：a, b 自带 1/2 因子（u 不变）；`c_val = q~` 供渲染端还原分母。
- gamma：见上节——N>0 时 m,n,au,bu,u 全是 ln_q 子证明参数，`cu_val = N`、`c_val = "0"`；N=0 时 a_val = 剩余常数、c_val 为方向标志。
- varpi/gauss 兜底：`cu_val ∈ {1,2}` 是模板标志（非 c·u），`a_val = "0"`，b_val = 附加常数项；solution 串恒为 `a = 0, b = 0`（站点原样）。
- `unified_form` 恒 {}。solution 串：`a = X, b = Y`；三系数族追加 `, c= Z`（`c=` 后无空格，站点原样）。

## 渲染层怪癖（/get_integral_image，逐字节实测）

- bound 与 coef **原样回显请求数位**（不约分：`3140/1000` → `\dfrac{3140}{1000}`），`\frac/\dfrac` 宏逐字回显；需要数值处宏按**无符数字**解析（`-\dfrac{6}{4}` 与 `\dfrac{6}{-4}` 都得 3/2）。
- 系数位只消去规范形 "1"：`coef=1/1` 印 `1\pi`、`01` 印 `01\pi`、`2/4` 印 `\dfrac{2}{4}\pi`。
- **零分子塌缩**：quadlog 族的 P=a+bx² 恒为乘积因子——a=b=0 时整个被积函数经 sympy 塌成裸 `0`（zeta3 印 `0\ln^2(x)`、catalan 印 `0\ln(1/x)`，皆实测）；golden 同样输出 `0`（实测）。
- ` \cdot ` 规则各族不同（各自对齐站端 sympy 版本）：quadlog/beta 族——左因子 `}` 结尾且右因子是 `\left(`+数字…`\right)`；exp/hyperbolic——左为 Pow 且右为数字开头的 Add；trig_q/trig_pi——左以数字或 `}` 结尾、右 `\left(`+数字或 `\frac{d}{d}` 开头。
- 空格微差：gamma '<' 排 `=\int`、`>' 排 `= \int`；pi '>' 句尾 `\mathrm{d}` 前双空格；pi_n 分指数 '>' 的 LHS 紧致排 `(π^{p/q})^q- (bound)^q`。
- ln 尾因子：catalan `\ln(1/x)`、zeta3 `\ln^2(x)`、pi_n `(\ln(1/x))^{r}`（r=0 省略；指数为负不加花括号）。
- golden 分子序：bu<0 时 P 排在 `\sqrt{x + 4}` 之前，否则在后。
- tan_q/cot_q 把 `1/cos(q)`（`1/sin(q)`）以 `\dfrac{1}{\cos(n/d)}}` 前置于积分，q 用 raw n/d 文本（宏坍缩成无符号 `d/d`）；负 q 整体前加 `- `。tanh_q/coth_q 同构，但 `\cosh(...)`/`\sinh(...)` 内嵌的是 coef **原文**（宏逐字回显，不归一成 n/d）。
- ln 族分母：s≥2 → `u'(dx+e)^s`（u'=1 省略）；s=1 → 展开 `Ax+B`（sympy 项序，如 `4 - 10 x`）。gamma 子分数同此精神（s'=1 → `u'Nx + u'`）。
- 字段校验（站点实测）：m,n ∈ [0,30]，u_val ≥ 1（gamma 放宽到 0）；a/b/c_val 走与 /calculate 相同的 'n/d' 语法；coef/rational 完全不校验、原文进 LaTeX。

## 关键数学事实（实现矩计算时用）

- `∫₀¹ x^k/(1+x²)dx`：偶 k → `±π/4 + 有理`；奇 k → `±(ln2)/2 + 有理`。递推 `J_k + J_{k-2} = 1/(k-1)`。
- `∫₀¹ x^{2k} ln^r(1/x)/(1+x²) dx` → `Q·β(r+1) + Q`（β(1)=π/4, β(2)=C, β(3)=π³/32, β(2k+1)∈Q·π^{2k+1}）。
- `∫₀¹ x^{2k+1} ln^r(1/x)/(1+x²) dx` → `Q·η(r+1) + Q`（η(s)=(1-2^{1-s})ζ(s)：η(2)=π²/12, η(3)=3ζ(3)/4, η(4)=7π⁴/720…）。
- `∫₀¹ x^k e^{qx} dx`：IBP 递推 `I_k = e^q/q − k/q·I_{k-1}` → `Q(q)e^q + Q(q)`。
- `∫₀¹ x^k/(1+cx)^s dx`：部分分式 → `Q(c)·ln(1+c) + Q(c)`。
- `∫₀^π e^x sin^j x dx`、 `∫₀^{π/2} sin^j x·sin(αx)dx`：积化和差展开逐项积分，全部落在 `Q·sin(πq)+Q` 或 `Q·e^{απ/2·…}+Q` 里（详见文章类型15/17/6的公式）。
- `∫₀¹ x^{4k+r}/√(1-x⁴) dx = (1/4)B((4k+r+1)/4, 1/2)`：J_{4k}∈Q·ϖ、J_{4k+1}∈Q·π、J_{4k+2}∈Q·G⁻¹、J_{4k+3}∈Q；基函数 x^{4m+res}(1−x) 取相邻两项，'>' 两型整体再除 π（符号映射 π→1、G⁻¹→ϖ⁻¹、ϖ→G，因 G=ϖ/π）→ 净落 span{r=1: Q+Q·ϖ⁻¹；r=3: Q+Q·ϖ；r=0: Q·G+Q；r=2: Q·G⁻¹+Q}。Γ(1/4) 换算见文章 35/36。
- `∫₀¹ x^k √(4+x) dx`：直接公式 → `Q + Q·√5`。

## Padé 插值法（ln q、arctan q 的第二套方案，进阶教程文章）

- `ln(1+x)`：`l_n = P_{n,n}` 严格下界、`u_n = P_{n+1,n}` 严格上界（Topsøe）。误差函数 `s_n = (ln-l_n)' = x^{2n}/((1+x)D_n²)`、`t_n = (u_n-ln)' = x^{2n+1}/((1+x)D̃_n²)` 恒正。
- 证 `ln(1+q) > p`：找最小 n 使 `l_n(q) ≤ p`，取 m<n 使 `l_m(q) > p`，插值 `a+b=1, a·l_n+b·l_m = p` → `ln(1+q)-p = ∫₀^q (a·s_n+b·s_m)dx > 0`。
- arctan 同理：`l_n = P_{2n,2n}` 下界、`u_n = P_{2n+1,2n+1}` 上界，`s_n,t_n` 为对应导数（恒正有理函数）。
- 网站实测的 ln 输出用的是 `(1+cx)^s` 核（非 Padé），Padé 方案作为备选/对照实现。

## 组合不等式 /decompose_inequality

输入如 `pi^2+8*pi>35`、`e*pi+phi+sin(1)<11`。实测输出：
- 归一化 `normalized_latex`；拆成各基本类型子界（`decomposition_latex`）；逐步 `steps[]` 各带完整积分式。
- 拆法：各项分配有理数界使方向一致、和/积覆盖目标。`pi^2+8pi>35` → `π²>227/23, 8π>578/23`（公分母 23）；`e*pi+phi+sin1<11` → `e<193/71, π<30317/9650, φ<8173/5050, sin1<85/101`（乘积项拆成两个因子界的积 193/71×30317/9650 < 该乘积分到的界）。
- 需要探测推断选界策略：疑似在每项的"可证最好界"（连分数/预算内最紧）附近分配。由 benchmark 数据反推。

## API 协议（实测）

```
POST /calculate  (urlencoded: type, power, comparison, rational)
→ {success, type, parameters: {m,n,a_val,b_val,c_val,au_val,bu_val,cu_val,u_val,unified_form}, equations:{solution}}
GET /get_integral_image?m&n&a_val&b_val&c_val&u_val&au_val&bu_val&cu_val&comparison&rational&coef&type
→ {equation: "LaTeX"}
POST /decompose_inequality (problem)
→ {success, problem, normalized_latex, decomposition_latex, basic_count, steps:[{type,label,coefficient,comparison,bound,bound_latex,equation}], direct_basic}
错误：{success:false, error:"要证明的式子不等号方向反了" | "在指数不超过{上限}的范围内未找到<方向的解"（上限=10，e/pi 为 30） | ...}
```
