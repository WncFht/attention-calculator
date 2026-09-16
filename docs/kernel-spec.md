# 注意力计算器复现 — 数学规格

来源：作者知乎专栏文章（本机归档 `~/src/obsidian/output/zhihu-mathematical/量化调酒师的数学文章/`；Mac 上是 `~/Desktop/obsidian/...`，别混）+ 对 zhuyidao.net 线上 API 的实测。

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
- `P` 是带 2~4 个待定系数的小多项式 `a+bx` 或 `a+bx+cx^2`（偶核用 `a+bx^2`；lemniscate 两型用 `a+bx^4`；[0,π] 与 [0,π/2] 域用 `a+b·sin x`；exact-only 的 `ln_q_cube` 升到三次 `a+bx+cx^2+dx^3` 配 4 维矩空间）；
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
- **ln_q_square 崩溃机理**（行列式级解释，exact 侧推导）：(m,n)=(1,0) 档 det `∝ −(c−4)(c+1)/(72c⁶)`、(1,1) 档 `∝ −(c−6)(c+1)³/(576c⁸)`、(2,1) 档 `∝ (c−12)`，c=q−1——q=5,7 恰好把浅档行列式清零（站端 500 的根因），q=13 可预言同型崩溃（未探测核实）；det(2,0)/det(3,0) 的根非有理。推论警示：q=5 '>' 且界 ≲0.64 时站端走的不是 500 而是会发证明——崩溃清单别过度概括。
- **发射自检**：`exact_check.verify(kind, power, comp, bound, params)` 对每条发射证明做 ℚ 字典相等复核（`combine(coeffs, true_basis) == target` 且 `poly_nonneg`）；失败是自家 bug → InternalError。exact_check 与模式无关：对 site 模式参数跑它即可复核站点输出的真伪（verify.py 的 ℚ 精确版）。
- **exact-only 类型**：`kernels.EXACT_TYPES` 注册无站端对应的新类型（当前 26 型：zeta5/7/9/11、beta4/6/8/10、ln_q_cube、ln_q_quad、arcsin_q、arsinh_q、gaussint_q、dawson_q、erfiint_q、pi_sqrt2、pi3、pi3_u、pi3_a、li2_q、psi1_q、si_q、cin_q 共 23 个矩核型 + gamma14、gamma34、gamma12 三个复合命题型，逐型规格见下「exact-only 型清单」）。/calculate 仅在 mode=exact 下接受（site 模式仍按站端口径 400）；solve.prove 对不带 exact 的调用报 ValueError；certified_cmp 常数走 `integrand.constant_mpf` 表。/get_integral_image 仍只认站端 29 型——exact 类型的渲染经 kernel 的 render_equation 程序化调用。

### exact-only 型清单

- `zeta5` / `zeta7` / `zeta9` / `zeta11`（kernels/zeta_odd.py）：`power·ζ(s)`，s ∈ {5,7,9,11}，复用 quadlog ln_moment 的 r=s−1 奇支路，η 系数乘 1−2^{1−s} 得 ζ（r=4 → 15/16 出 ζ(5)，r=10 → 1023/1024 出 ζ(11)）——实现细节见核 docstring；ζ(13)+ 同构可扩（仅加 ZETA 表行）。
- `beta4` / `beta6` / `beta8` / `beta10`（kernels/beta_even.py）：`power·β(s)`，s ∈ {4,6,8,10}，power 是系数。复用 ln_moment 的**偶支路** r=s−1：矩 (−1)^k r!·(β(r+1)−partial) 的符号列即 β(r+1) 本身，无归一化系数（与 zeta_odd 的 η→ζ 换算不同）。核 `x^{2m}(1−x²)^n(a+bx²)ln^r(1/x)/(1+x²)`；r 为奇，渲染必须印 `ln^r(1/x)`（`ln^r(x)` 差负号）。LIMIT=10：'<' 侧证到 ~5e-12（β4 (9,9)），'>' 侧 ~3e-14（β6 (9,10)）；β(4) '>' 7.7e-14 需 (11,11) 超预算。β(12)+ 同构可扩（仅加 BETA 表行），未实装。实现笔记 `docs/2026-09-16-beta-even-impl-notes.md`。
- `ln_q_cube`（kernels/ln_pow.py）：`(ln q)³ ⋚ bound`，q>1 由 power 槽携带。核 `ln²(1+cx)/(1+cx)^s`（s=max(m,n,1)，c=q−1），矩空间 4 维 span{ln³q, ln²q, ln q, 1}——ln 族"核低一次、u⁻¹ 项产目标常数"的同构升幂。P 升三次 a+bx+cx²+dx³（4 符号需 4 系数），参数在九槽基础上扩 `d_val`/`du_val`；非负判据为 QQ(√D) 精确三次规则（临界点 P(x*±)=A∓B₀√D，qsign 原语判号），推导见 docs/2026-09-16-ln-cube-derivation.md。
- `ln_q_quad`（kernels/ln_pow.py）：`(ln q)⁴ ⋚ bound`，q>1 由 power 槽携带。核 `ln³(1+cx)/(1+cx)^s`，矩空间 5 维 span{ln⁴q,…,1}——同构升幂的第二节。P 升四次（5 符号 5 系数），参数再扩 `e_val`/`eu_val`；非负判据为 Sturm 奇根计数：最低次非零系数 >0 且开区间 (0,1) 无奇数重根（剥端点因子后奇次 sqf 因子之积的 `count_roots`），充要条件。LIMIT=10 与 ln 族一致：q→1 浅档可证 15 位小数级紧界，q=2 最深 (10,10)，q=10 紧界诚实 NoSolution（预算问题，非正确性）。推导见 docs/2026-09-16-ln-quad-impl-notes.md。
- `arcsin_q`（kernels/arcsin.py）：`arcsin q ⋚ bound`，q∈(0,1) 有理数。核 `1/√(1−q²x²)`（arcsin q = ∫₀¹ q·dx/√(1−q²x²)），矩递推 `M_k = (k−1)/(kq²)·M_{k−2} − w/(kq²)`（w=√(1−q²)，M₀=arsin(q)/q、M₁=(1−w)/q²）落 span{arcsin q, w, 1}，w 为寄生常数——同 sin_q 的 cos_q 模式。退化档：`1−q²` 为有理平方（3-4-5 族 q=3/5, 4/5, 5/13…）时空间坍缩成 span{arcsin q, 1}，w 折入有理项按 2 维解（P=a+bx）——q=4/5 符号系任何预算都不可证、折叠后 (4,8) 可解。覆盖窗口不对称：q→1 需 m~1/(1−q)，q≳0.8 超预算。LIMIT=30。
- `arsinh_q`（kernels/invhyp.py）：`arsinh q ⋚ bound`，q≠0（奇函数，负 q 走 |q|）。核 `1/√(1+q²x²)`（arsinh q = ∫₀¹ q·dx/√(1+q²x²)），递推 `√(1+q²) = (k−1)·M_{k−2} + kq²·M_k`（k≥2）落 span{arsinh q, √(1+q²), 1}，P=a+bx+cx² 三系数消寄生根号。q=0 在 M₀ 的 1/q 上自然死亡。**arcosh_q 未实现**：其矩空间 {arcosh q, √(q²−1), √(2(q−1))} 无 "1" 方向，有理界进不了恒等式（详见核 docstring）。
- `gaussint_q` / `dawson_q` / `erfiint_q`（kernels/gauss_erf.py）：erf 本体不可及——每个矩的 erf 分量带 √π 无法被吸收（docs/2026-09-16-w3-research-erf.md §√π障碍），改为证三个缩放常数：`G(q)=∫₀^q e^{−t²}dt = √π·erf(q)/2`（核 `e^{−q²x²}`）、`F(q)=e^{−q²}∫₀^q e^{t²}dt` Dawson（核 `q·e^{q²(x²−1)}`，归一化使 F 本身而非 e^{q²}F/q 进 span）、`H(q)=∫₀^q e^{t²}dt = √π·erfi(q)/2`（核 `e^{q²x²}`）。三者共用一条两步 IBP 链（t=1/(2q²)，边界项落寄生符号 e^{±q²}），矩空间 span{C_q, e^{±q²}, 1}，P=a+bx+cx² 消寄生。皆奇函数：负 q 折进目标符号 σ=sign(q)（σ·C(|q|)=C(q) 直接回显原命题）；q=0 全退化拒收。LIMIT=10。
- `pi_sqrt2`（kernels/pi_sqrt2.py）：`q·S ⋚ r`，S = π√2 = Γ(1/4)Γ(3/4)（余元公式），lemniscate 格点 `(1−x⁴)^{c−1}` 族仅剩的干净单核型。'>' 核 `(1−x⁴)^{-3/4}`、基 `x^{4m+2}(1−x)(a+bx⁴)`，矩 `{S: u_{m+i}/4, 1: −v_{m+i}}`（`u_j=(3/4)_j/j!`、`v_j=(1)_j/(5/4)_j`）；'<' 核 `(1−x⁴)^{-1/4}`、基 `x^{4m+3}(1−x)(a+bx⁴)`，矩 `{1: w_{m+i}, S: −z_{m+i}}`（`w_j=(1/3)(1)_j/(7/4)_j`、`z_j=(1/16)(5/4)_j/(2)_j`）。无倒数归一、无兜底模板；`n` 为方向标志（同 beta_family）。power 为纯系数不进核；q≤0 免校验。单轴 m 搜索 LIMIT=256：渐近 `u_m/v_m→4/S`、`z_m/w_m→1/S` 保证真界可达。实现笔记 `docs/2026-09-16-pi-sqrt2-impl-notes.md`。
- `pi3` / `pi3_u` / `pi3_a`（kernels/dixon.py）：Dixon/B(1/3,1/3) 格点族，常数 π₃ = B(1/3,1/3) ≈ 5.2999、U = Γ(2/3)²/Γ(1/3) ≈ 0.6845、A = 2π/√3 ≈ 3.6276，power 是系数 q，q=0 拒收。核 `K = (1−x³)^{−k/3}`（k∈{1,2}），矩 `T_k = (1/3)B((k+1)/3, 1−k/3)` 按 k mod 3 分三类；基底 `x^{3m+res}(1−x)(a+bx³)` 配对相邻残类成 2 维 span，'>'-倒数方向乘外因子 `√3/(2π)=1/A`（varpi `/π` 的三次类比）做符号改名。六个 (kind,comp) 格与 (kernel,res,outer) 配置双射；span 集 {pi3, pi3_inv, U, A, G3, 1}。参数 `m` 为真实指数级、`n` 是方向标志、`cu_val=1` 为 '+b' 兜底；单轴 m=0..512 搜索，边比 ~1/m 收敛。渲染倒数方向印倒数 LHS（同 varpi '>'/gauss '<' 例）。实现笔记 `docs/2026-09-16-dixon-impl-notes.md`。
- `li2_q`（kernels/li2.py）：`Li₂(q) ⋚ bound`，q∈ℚ、q<1、q≠0 由 power 槽携带（q=1 时 ln 列塌缩、q>1 支点落入域内、q=0 核恒零，均拒收）。核 `K=σ(φ−L)`，φ=−ln(1−qx)/x，σ=±1 按方向；矩空间 span{Li₂(q), ln(1−q), 1} 恰三维，Li₂ 只在 A₀ 出现 → 基强制 (1−x)ⁿ、P 常数项钉死 +1。矩 `A₀=Li₂(q)`、`A_k=(q^{−k}−1)/k·ln(1−q)+q^{−k}/k·Σ_{l≤k}q^l/l`。位移族：q>0 '>' 用 taylor p_d、'<' 用 chord 与 lift；q<0 时 φ 严格凹，'>' 用 chord + alt 奇数阶、'<' 用 alt 偶数阶或 tan0。搜索轴 n≤16 × (fam,d≤64)，参数扩 `fam`/`d` 两键；EXPONENT_LIMIT=16。q→1⁻ 是已知弱区（间隙 ~q^d）；q<−1 只有线性见证。实现笔记 `docs/2026-09-16-li2-impl-notes.md`。
- `psi1_q`（kernels/trigamma.py）：`ψ′(q) ⋚ bound`，q>0 有理数由 power 槽携带（系数恒 1，constant_mpf 不乘 q）。方向核 '>'：`x^{q−1}(x−1−lnx)/(1−x)`；'<'：`x^{q+N−2}(1−x+xlnx)/(1−x)`，N=0（q>1）或 1（0<q≤1），两核定号均为 ln t≤t−1 的改写。矩 `∫x^uK_> = ψ′(q)−L_u`、`∫x^uK_< = U_u−ψ′(q)` 落 span{ψ′(q),1}；(1−x)^n 消 ψ′ 行 → 一维 m 扫描，P=a+bx 恒 a+b=1，全真命题可证无常数尾。EXPONENT_LIMIT=256。实现笔记 `docs/2026-09-16-trigamma-impl-notes.md`。
- `si_q` / `cin_q`（kernels/sicin.py）：`Si(q) ⋚ bound`、`Cin(q) ⋚ bound`，q≠0 由 power 槽携带；Si 奇（σ=sign q）、Cin 偶（σ=+1），负 q 折进目标符号与核方向。四维 span{C_q, sin_q, cos_q, 1}——C_q 只在基项 l=i=j=0 出现，故 P 必须三次（参数扩 `d_val`/`du_val` 并新增 `t_val`），且 '<' 方向裸核下原理不可证、必须用互补核。核为 Taylor 余量阶梯 t=0,1,2（交错截断引理保证 (0,1] 严格正）；基 `(1+x)^m(1−x)^n`——x^m(1−x)^n 在 m≥1 使 C_q 行全零恒奇异。P(0)=+1 被目标方程强制 → WrongDirection 不可达，假命题统一 NoSolution；非负判据复用 ln_pow.cubic_nonneg。LIMIT=10、T_LEVELS={0,1,2}；覆盖短板：q 深振荡区 '>' 需 t~q/2π。实现笔记 `docs/2026-09-16-sicin-impl-notes.md`。
- `gamma14` / `gamma34` / `gamma12`（kernels/gamma_special.py）：复合命题型而非矩核型——Γ(1/4) 被 quarter-lattice 奇偶引理锁在偶次幂（docs/2026-09-16-w3-research-gamma-quarter.md），落地走作者原机制：g²=2ϖ·√(2π) 拆成一条 varpi 界 + 一条 pi 界再开方转移；Γ(3/4)=π√2/g 多一级除法转移；Γ(1/2)=√π 单条 pi 界开方。证书是 `rule`+`witness`+`children` 的小证明 DAG（规则 sqrt_mul/sqrt_div/sqrt/pos_transfer），certificate.verify_cert 递归复核子证书并在 ℚ 上重算转移算术。

### exact 路径的预算与域差异（与站端口径对照）

- **搜索上限**：各核 `LIMIT` 默认 10（EXPONENT_LIMIT 例外表：pi/e=30 同站端，arcsin_q=30 因 q→1 需 m~1/(1−q)，pi_sqrt2=256，pi3 三型=512，li2_q=16 为 n 轴上限，psi1_q=256）；beta 族 '<' 在 exact 下不再走转置档、上限放宽为 `EXACT_LT_LIMIT=256`（beta_family.py），更深的真证明仍诚实 NoSolution。
- **pi_n exact 域**（quadlog.py）：`power.numerator < 1 or power.denominator > 64` → ValueError 域拒；分子不设上界——exact 侧新增 beta_pi/eta_pi 生成器（Euler/Bernoulli 系）支持表外指数，site 侧仍限查表。
- **arctan_q / arccot_q exact**：power=0 → ValueError 域拒（矩基例带 1/q；arctan 0=0 有理、arccot 0=π/2 不在 span）。
- **hyperbolic exact**（hyperbolic.py）：q<0 先按奇偶归约到 |q| 再解（cosh 偶、sinh/tanh/coth 奇——site 路径同型直接 WrongDirection，exact 修掉对真命题撒谎）；q=0 退化输入域拒。
- **方向判定**：`solve.certified_cmp` 用 mpmath 80→240→800→2400 dps 递增 + 护栏带 `2^(30−dps)·max(1,|c|)` 认证 sign(C−r)，返回 +1/−1/0/None；None（不可实值求值）交给核内域校验。

### exact_check 复核器契约（checker 作者必读，4 个 agent 独立踩过的坑）

`exact_check.<family>.check(kind, power, comp, bound, params) -> {"identity_ok","nonneg","integrand","target"}`。两条硬性规则：

1. **Moment dict 必须滤零值键**：`combine/add/scale` 产出的 dict 省略零系数；手搭 target 字面量若保留 `"1": Fraction(0)` 会造成字典相等误判——构造后过 `{k: v for k, v in t.items() if v}`。
2. **nonneg 必须含非恒零守卫**：`identity_ok=True` 挡不住 `0 = ∫0 dx > 0` 式空洞"证明"（站端 golden/power=0 记录里真有这样的成功记录——等式成立但命题为假）。nonneg 判定要并入 `bool(integrand)`（或等价的非零被积函数检查）。

另：trig_pi 复核时若记录的 `c_val`（α）与精确 α 不一致，不要强行重映射 `"C"` 键——让恒等式诚实失败，而不是凑出匹配。

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
- 证 `ln(1+q) > p`：`l_n` 单调递增至 ln——取同侧最弱下界 `l_1` 与首个越界下标 `m`（`l_m(q) > p`，必 m>n），插值 `a+b=1, a·l_1+b·l_m = p` → `ln(1+q)-p = ∫₀^q (a·s_1+b·s_m)dx > 0`。（本节旧文"取 m<n"方向写反，以 `docs/2026-09-16-pade-notes.md` 与实现为准；`p<l_1(q)` 时无同侧逼近可配，退化为单项证书 `resid = l_m(q)-p > 0`。）
- arctan 同理：`l_n = P_{2n,2n}` 下界、`u_n = P_{2n+1,2n+1}` 上界，`s_n,t_n` 为对应导数（恒正有理函数）。
- 网站实测的 ln 输出用的是 `(1+cx)^s` 核（非 Padé）。**落地状态（2026-09-16）**：`pade.py` 已作为 mode=exact 的在线兜底——`prove_exact` 对 ln_q/arctan_q 在 (m,n) 搜索耗尽后回落 Padé（预算 MAX_N=50），响应 `prover:"pade"`、无 parameters、证书为 Padé schema（`docs/2026-09-16-certificate-spec.md` §Padé 证书变体）。

## AGM 区间法（gauss/varpi 的第二套方案）

`agm.py`（顶层，与 pade.py 同级）：G = 1/M(1,√2) 的 AGM 迭代 `a_{n+1}=(a_n+b_n)/2, b_{n+1}=√(a_n·b_n)` 用 isqrt 有理包络（a_n↑M≤b_n↓ 三明治），全程 Fraction/int 无浮点。gauss 在 (m,n) 耗尽后直出区间证书（`prover:"agm"`，`agm_iter`/`agm_digits`/`lo`/`hi` 字段，复核器重放到逐字相等）；varpi = π·G 走 composite DAG（gamma_special 新规则 `pi_div_agm`：pi 子证 + AGM 子证，见证 (A,B)）。实测地板：gauss 双向 ~1e-2000 间隙（digits=4096 档，11 步迭代），varpi ~1e-16（pi oracle 地板所限）。细节 `docs/2026-09-16-agm-impl-notes.md`。

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
