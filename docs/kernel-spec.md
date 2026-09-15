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
- `P` 是带 2~3 个待定系数的小多项式 `a+bx` 或 `a+bx+cx^2`（偶核用 `a+bx^2`）；
- `K` 是核函数，使所有"矩" `∫ B·x^k·K dx` 落在 `span_Q{1, C, D_1, ...}` 有限维空间里。

求解：矩对 (a,b,c) 线性 → 在 Q 上解小线性方程组
- C 的系数 = ±1（± 号由方向定）
- 每个寄生常数 D_j 系数 = 0
- 常数项 = r

然后检查 `P` 在域上不变号；变号则换下一组 (m,n)。

## 搜索顺序（作者亲述）

1. `m+n` 升序；
2. 同 `m+n` 时 `|m-n|` 小的优先（即先取平衡的）；
3. e、π 两类型指数上限 30；其余类型上限 10（报错文案"在指数不超过10的范围内未找到…方向的解"）。

同方向连续失败时若算出的被积函数恒≤0，返回"要证明的式子不等号方向反了"。

## 非负性判据（作者亲述，只需判断 P）

- 一次 `a+bx`：`a≥0 且 a+b≥0`。
- 二次 `a+bx+cx^2`：`c≤0` 时同上（端点非负即可）；`c>0` 时若对称轴 `-b/2c ∈ (0,1)` 需 `4ac-b²≥0`，否则端点非负。
- 偶核 `a+bx^2`：等价于 `a+b·t` 在 `t∈[0,1]` 非负 → `a≥0 且 a+b≥0`。
- `a+b·sin x`（域 [0,π] 或 [0,π/2]，sin x∈[0,1]）：同一次型判据。

## 类型表（网站 29 种 → 核与矩空间）

记 `M(...)` 为"矩 ∈ span{...}"。`s` 为分母幂（ln 类可升幂加速收敛，见技巧6）。

**ln 类分母幂实测规则**：`s = max(m, n, 1)`（61 条 ln_q golden 逐一数值拟合全中；s 不在响应参数里，只体现在渲染方程与矩构造中）。

| type | 域 | 核 K(x) | 基 B_{m,n} | P | 矩空间 |
|---|---|---|---|---|---|
| pi | [0,1] | 1/(1+x²) | x^{2m}(1-x²)^n | a+bx² | {π, 1} |
| e | [0,1] | e^x | x^m(1-x)^n | a+bx | {e, 1} |
| pi_n (π^k) | [0,1] | ln^{k-1}(1/x)/(1+x²) | x^m(1-x²)^n，m+k 奇 | a+bx² | {π^k, 1} |
| e_q (e^q) | [0,1] | e^{qx} | x^m(1-x)^n | a+bx | {e^q, 1} |
| ln_q (q>1) | [0,1] | 1/(1+(q-1)x)^s | x^m(1-x)^n | a+bx | {ln q, 1} |
| ln_q_square ((ln q)²) | [0,1] | ln(1+(q-1)x)/(1+(q-1)x)^s | x^m(1-x)^n | a+bx+cx² | {ln²q, ln q, 1} |
| sin_q | [0,1] | sin(qx) | x^m(1-x)^n | a+bx+cx² | {sin q, cos q, 1} |
| cos_q | [0,1] | sin(qx) | x^m(1-x)^n | a+bx+cx² | {sin q, cos q, 1} |
| tan_q | [0,1] | sin(qx)，结果整体除 cos q | x^m(1-x)^n | a+bx+cx² | 解 sin q − p·cos q |
| cot_q | 同 tan | 同 tan（先取倒数） | | | |
| sin_q_degree | [0,π/2] | sin((1-2q)x)，q=度数·π/180 | sin^m x(1-sin x)^n | a+b sin x | {sin(πq), 1} |
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

实测样例（线上返回）：

- `pi < 22/7` → `22/7 - π = ∫₀¹ x⁶(1-x²)³(47-13x²)/(120(1+x²)) dx > 0`（m=3,n=3，P=47-13x²，u=120）
- `e > 8/3` → `e - 8/3 = ∫₀¹ x²(1-x)eˣ/3 dx > 0`（m=1,n=1，P=x/3）
- `ln 3 < 11/10` → `∫₀¹ x²(1-x)³(40x+6)/(45(2x+1)³) dx`（s=3）
- `e^π > 23` → `∫₀^π (1-sin x)⁵(547950-422240 sin x)eˣ sin²x/336633 dx`
- `ϖ < 8/3` → `∫₀¹ x²⁷(3319635x⁴+2432199)(1-x)/(353600√(1-x⁴)) dx`
- `G > 4/5` → `∫₀¹ x⁸(44x⁴+6)(1-x)/(5π√(1-x⁴)) dx` —— 即文章示例
- `sin(π/5) > 1/2` → `∫₀^{π/2} (1-sin x)(637 sin x+125) sin(3x/5)/750 dx`（1-2q=3/5）
- `γ < 3/5` → `∫₀¹ [x⁵((2-x)/2-1/(1-x)-1/ln x) + x²(1-x)⁴(18460x+2811)/(168(5x+1)⁴)] dx`（前半给 γ+ln6 项、后半是 ln6 上界子证明）

## 关键数学事实（实现矩计算时用）

- `∫₀¹ x^k/(1+x²)dx`：偶 k → `±π/4 + 有理`；奇 k → `±(ln2)/2 + 有理`。递推 `J_k + J_{k-2} = 1/(k-1)`。
- `∫₀¹ x^{2k} ln^r(1/x)/(1+x²) dx` → `Q·β(r+1) + Q`（β(1)=π/4, β(2)=C, β(3)=π³/32, β(2k+1)∈Q·π^{2k+1}）。
- `∫₀¹ x^{2k+1} ln^r(1/x)/(1+x²) dx` → `Q·η(r+1) + Q`（η(s)=(1-2^{1-s})ζ(s)：η(2)=π²/12, η(3)=3ζ(3)/4, η(4)=7π⁴/720…）。
- `∫₀¹ x^k e^{qx} dx`：IBP 递推 `I_k = e^q/q − k/q·I_{k-1}` → `Q(q)e^q + Q(q)`。
- `∫₀¹ x^k/(1+cx)^s dx`：部分分式 → `Q(c)·ln(1+c) + Q(c)`。
- `∫₀^π e^x sin^j x dx`、 `∫₀^{π/2} sin^j x·sin(αx)dx`：积化和差展开逐项积分，全部落在 `Q·sin(πq)+Q` 或 `Q·e^{απ/2·…}+Q` 里（详见文章类型15/17/6的公式）。
- `∫₀¹ x^{4k+r}/√(1-x⁴) dx = (1/4)B((4k+r+1)/4, 1/2)`：r=1 → `Q + Q·ϖ⁻¹`；r=3 → `Q + Q·ϖ`；r=0 → `Q·G + Q`（除以π）；r=2 → `Q·G⁻¹ + Q`。Γ(1/4) 换算见文章 35/36。
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
错误：{success:false, error:"要证明的式子不等号方向反了" | "在指数不超过10的范围内未找到<方向的解" | ...}
```
