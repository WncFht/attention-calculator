# W3-li2 实现记录：li2_q kernel + checker（2026-09-16）

实现 `src/attention_calculator/kernels/li2.py` 与 `src/attention_calculator/exact_check/li2.py`，调研依据 `docs/2026-09-16-w3-research-li2.md`。本文记录参数编码、符号名、对调研文档的偏离、以及建议并入 `docs/kernel-spec.md` 的条目草稿。

## 与调研文档的偏差

1. **q<0 凸性结论成立，但引理方向以本文为准。** 调研文档"q<0 时 φ 凹"经复核正确：x³φ″(x) = B(y)，y = −qx > 0，B(y) = y²/(1+y)² + 2y/(1+y) − 2ln(1+y) = (3y²+2y)/(1+y)² − 2ln(1+y)，B(0)=0 且 **B′(y) = −2y²/(1+y)³ < 0**（推导中易错成 +2y²/(1+y)³，符号反了会得到凸的错误结论），故 B(y) < 0、φ″ < 0 对一切 u=−q>0 成立——凹性无 |q| 范围限制。q>0 时 φ″>0 由级数逐项成立。
2. **lifted 族未实现，由更弱的替代覆盖。** 调研文档的 q<0 升级族 chord+t·x(1−x) 与 tan0+t·x² 需要证 t ≤ t*(q)（一个一维超越最值），无干净的精确符号论证，按预定方案不实现。替代：−1≤q<0 时交错级数给出 d 轴几何收敛族（φ = Σ(−1)^i u^i x^{i−1}/i，ux≤1 时逐项递减，奇数阶截断 p_d ≤ φ 证 '>'，偶数阶 p_d ≥ φ 证 '<'，d=2 即 tan0）。实测覆盖优于原 lifted 族：q=−1/2 两方向 1e-6 界在 d≤12、n=0 即解（文档原族只到 ~1e-4）。q<−1 处级数在 x>1/u 发散，只剩 chord（'>'）/tan0（'<'）线性见证，覆盖有限属已知限制。
3. **'>' q>0 的"微紧化"不需要单列**：p_d + 尾项首项 q^{d+1}x^d/(d+1) 就是 p_{d+1}，taylor 族 d 轴已含。

## 参数编码

emit 九槽之外新增两键（kernel.prove 写入 parameters，checker 与 renderer 读取）：

- `fam`（str）：位移多项式族标签，同时是 checker 的符号引理键——`taylor` / `lift` / `alt` / `chord` / `tan0`。
- `d`（int）：族参数（截断次数；chord 恒 0，tan0 恒 2）。

m 恒为 0（Li₂ 只在 A_0 的结构推论），a_val 恒 "1"（方向由核符号因子 σ=±1 吸收，P 不翻转）。

## 矩符号表

`{"li2_q", "ln_1mq", "1"}`：`li2_q` 即目标常数 Li₂(q)（与 kind 同名，随 gauss_erf 惯例）；`ln_1mq` 为 ln(1−q)。cert 渲染建议 `_SYMBOL_TEX` 加 `"li2_q": "\\mathrm{Li}_2"`、`"ln_1mq": "\\ln(1-q)"`。

## kernel-spec.md 条目草稿（请并入 exact-only 清单节）

- `li2_q`（kernels/li2.py）：`Li₂(q) ⋚ bound`，q∈ℚ、q<1、q≠0 由 power 槽携带（q=1 时 ln 列塌缩为 span{π²/6,1} 且 '<' 本质不可能，q>1 支点落入域内，q=0 核恒零，均拒收）。核 `K=σ(φ−L)`，φ=−ln(1−qx)/x，σ=+1（'>'）/−1（'<'）；矩空间 span{Li₂(q), ln(1−q), 1} 恰三维，Li₂ 只在 A₀ → 基强制 (1−x)ⁿ（m≡0）、P 常数项钉死 +1，(b,c) 由 {ln,1} 列 2×2 解出（search 层仍走 3×3 solve_moment，j=0 行给出 a=1）。矩公式 A₀=Li₂(q)、A_k=(q^{−k}−1)/k·ln(1−q)+q^{−k}/k·Σ_{l≤k}q^l/l。位移族与符号论证（checker `certified` 表即非负引理）：'>' q>0 用 taylor p_d（尾项逐项 ≥0）；'<' q>0 用 chord（凸函数弦在上方）与 lift p_d+C′x^d、C′=q^{d+1}/((d+1)(1−q))（尾项放缩）；q<0 时 φ 严格凹（x³φ″=B(−qx)，B′(y)=−2y²/(1+y)³<0），'>' 用 chord + alt 奇数阶（−1≤q），'<' 用 alt 偶数阶（−1≤q）或 tan0（任意 q<0，覆盖 q<−1）。搜索轴 n≤16 × (fam,d≤64)；参数扩 `fam`/`d` 两键。P 沿用站点二次非负判据。q→1⁻ 是已知弱区（间隙 ~q^d，q=0.99 的 1e-3 界超预算）；q<−1 只有线性见证。

## 注册需求（leader 接线）

- `kernels/EXACT_TYPES += "li2_q"`；`solve.FAMILY["li2_q"] = "li2"`。
- `integrand.constant_mpf`：`"li2_q": lambda: mp.polylog(2, q)`（power 是 q 本身，不是系数）。
- `solve.direction_f`：exact 型走 certified_cmp，无需 float 预检条目。
- `engine.EXPONENT_LIMIT["li2_q"] = 16`（failure_text 的"指数不超过 N"取 n 上限）。
- `certificate._SYMBOL_TEX`：见上方符号表。
- domain_error：check_input 已抛中文 ValueError（q=0 / q≥1 两条文案）。

## 已验证样例（测试内常驻）

- `Li₂(1/2) − 4/9 = ∫₀¹ (1 − 8x/3 + 16x²/9) φ dx`（taylor d=0, n=0；顶点贴零）
- `59/100 − Li₂(1/2) = ∫₀¹ (1 − 33x/25 + 18x²/25)(chord−φ) dx`（n=0）
- `Li₂(−1/2) + 9/20`、`−2/5 − Li₂(−1/2)`（alt d=3/d=2, n=0）
- q=−2 '>'/`<`（chord/tan0 线性见证，q<−1 路径）
- 1e-6 双向界：q∈{1/2, 2/3, 9/10, −1/2, −9/10, −1} 全解；q=9/10 需 d≈57。
