# W3-dixon 实现笔记：pi3 / pi3_u / pi3_a（Γ(1/3) 格三常数）

日期：2026-09-16。依据 `docs/2026-09-16-w3-research-gamma-third.md`（矩递推、六配置、覆盖结构均已在 100 dps 复核）落地。交付：`kernels/dixon.py`、`exact_check/dixon.py`、`tests/test_exact_dixon.py`（59 项全过），未提交注册表改动（registry 归 leader 合并）。

## 命名与符号

三个 exact-only 型，一型一常数、各覆盖双向：

- `pi3` — π₃ = B(1/3,1/3) = √3·Γ(1/3)³/(2π) ≈ 5.299916（Dixon 椭圆周期，OEIS A197374）。常数符号就叫 π₃，`pi3` 与 `varpi`/`gauss`/`zeta3` 同命名法。
- `pi3_u` — U = A/π₃ = Γ(2/3)²/Γ(1/3) ≈ 0.684463。U 无标准名（调研文档以 U 记之），`pi3_` 前缀标明 Γ(1/3) 格家族。
- `pi3_a` — A = B(1/3,2/3) = 2π/√3 ≈ 3.627599。**没有用建议名 `pi_sqrt3`**：`pi_sqrt2` 已落地（π·√2，见 kernels/pi_sqrt2.py），该命名法下 `pi_sqrt3` 会读成 π√3 ≈ 5.4414——与 A = 2π/√3 ≈ 3.6276 不是同一个常数。改用族前缀 `pi3_` + 文档符号 `a`。

Moment 符号表（与调研文档一致）：`{"pi3", "pi3_inv", "U", "A", "G3", "1"}`，其中 `G3 = π₃/A = U⁻¹` 只出现在 `pi3_u '<'` 的 span 里，印刷 LHS 时印 `U^{-1}` 不印 G3。符号间换算（π₃·U=A、π₃·√3/(2π)=G₃、A·√3/(2π)=1、U·√3/(2π)=π₃⁻¹）编译期硬编进 `OUTER_FACTOR` 改名表——机制同 beta_family 的 `DIV_PI`（站端 '/π' 的位置由 √3/(2π)=1/A 扮演）。

## 配置与参数编码

六个 (kind, comp) 对到 (kernel, res, outer) 是双射（`CONFIG` 表）：

| (kind, comp) | kernel (1−x³)^(−k/3) | res | √3/(2π) | span | 证明形 |
|---|---|---|---|---|---|
| pi3 '<' | k=2 | 2 | — | {1, pi3} | r − q·π₃ |
| pi3 '>' | k=1 | 0 | ✓ | {1, pi3_inv} | q − r·π₃⁻¹ |
| pi3_u '>' | k=1 | 1 | — | {U, 1} | q·U − r |
| pi3_u '<' | k=2 | 0 | ✓ | {G3, 1} | r·U⁻¹ − q |
| pi3_a '<' | k=1 | 2 | — | {1, A} | r − q·A |
| pi3_a '>' | k=2 | 1 | — | {A, 1} | q·A − r |

基底 `x^{3m+res}(1−x)(a+bx³)`，P 的非负判据与 a+bx⁴ 同（端点规则 a≥0 ∧ a+b≥0，x³ 单调）。参数编码沿用 beta_family/pi_sqrt2：**n 是方向标志**（'>'→0、'<'→1），res 不进参数——(kind, comp) 已唯一决定 (kernel, res, outer)，checker 与 renderer 都从 `CONFIG` 取，伪造 n 只改簿记、不改被核查的被积函数。

`+b` 兜底模板（varpi/gauss 的 cu_val 机制移植版）：`∫ t·x^{3m+res}(1−x)·outer·K dx + b`，t 由目标常数系数定、b 是有理剩余；要求 t≥0 且 b>0 严格，扫描取矩比越过 bound 的最浅 m。参数形 `m`=实级、`n`=flag、`a_val="0"`、`b_val=b`、`au_val/u_val`=t（约分）、`cu_val="1"`、`solution` 恒 `"a = 0, b = 0"`。与站端差异：站端用固定单项式（m=n=0），本族扫 m 取首个命中——更完备且参数如实上报。

## 方向约定与偏差

- '>' 倒数方向沿用 varpi '>' 惯例：目标 `{1: q, pi3_inv: −r}`，印刷 `q − r·π₃^{-1}`。调研文档示例 `π₃>5` 印的是 `1/5 − π₃⁻¹`（按 bound 归一的写法）；本实现按 ϖ 先例印 `1 − 5π₃^{-1}`——同一命题的 5 倍缩放，命中的 m 与可证集相同（m=3 复核命中，系数 (113,156)/28 = 文档 (113,156)/140 的 5 倍）。
- `pi3_u '<'` 沿用 gauss '<' 惯例：目标 `{G3: r, 1: −q}`，印刷 `r·U^{-1} − q`；文档命题 `U < 5/7` 印成 `(5/7)U^{-1} − 1`（= 文档写法 `G₃ > 7/5` 的同命题异缩放），m=2 命中。
- 矩结构里 "1" 符号不会出现在 res=0 配对（残差类 0−1 配对落在 {A,U}/{π₃,A}），`OUTER_FACTOR` 改名表因此不含 "1" 键。

## 搜索行为（实测，均复核过）

- LIMIT = 512 单轴 m 扫描。边比 ~常数+O(1/m)（pi3 '<' 边比 ≈ 0.6/m，'>'-inverse 侧 ≈ 1.2/m）：1e-2 紧界 m≈60 内命中，2e-3 紧界 m≈340 内命中，8e-5 间隙（pi3 < 53/10）需 m≈7000 → 诚实 NoSolution。全 513 级扫描 ~0.4s。
- 覆盖完整性：对三常数双向 0.01–20 网格真命题全命中（主搜索或兜底），假命题全部 NoSolution——**本族 2×2 求解在假命题上从不产出非正 P**（1400+ 样本扫描零 WrongDirection），故核内 `search` 的 WrongDirection 路径实际不可达；管线侧假命题由 certified_cmp 预检统一报 WrongDirection。
- q≤0 系数：q=0 直接 ValueError（命题退化为有理数比较）；q<0 的真命题（如 −π₃<5）给不出非负 P，报 NoSolution——与 varpi 同约定。
- bound=0 的 '>' 命题：主搜索恒给出不定解，落兜底 t=0 分支（`∫0 + q > 0`，命题即其有理部分），恒等与非负均真。

## 注册表合并清单（leader 合并）

- `kernels/EXACT_TYPES += ["pi3", "pi3_u", "pi3_a"]`；`solve.FAMILY` 三型均 → `"dixon"`。
- `integrand.constant_mpf`：`pi3`: `q * mp.beta(mp.mpf(1)/3, mp.mpf(1)/3)`；`pi3_u`: `q * mp.gamma(mp.mpf(2)/3)**2 / mp.gamma(mp.mpf(1)/3)`；`pi3_a`: `q * 2 * mp.pi / mp.sqrt(3)`。
- `integrand.lhs_mpf`：`pi3` '>' → `1 - r/c`（varpi 同款倒数形）；`pi3_u` '<' → `r/c - 1`（gauss 同款）；其余走默认。
- `engine.EXPONENT_LIMIT`：三型均 512（failure_text 用；核内 LIMIT 自定）。
- `server.domain_error`：`kind in ("pi3","pi3_u","pi3_a") and power == 0` → `"系数不能为0"`。
- `certificate._SYMBOL_TEX`：`"pi3": "\\pi_{3}"`、`"pi3_inv": "\\pi_{3}^{-1}"`、`"U": "U"`、`"A": "A"`、`"G3": "G_{3}"`（缺省回落 `\mathrm{sym}` 也能用，只是丑）。

## kernel-spec.md 小节草稿（供合并）

### pi3 / pi3_u / pi3_a（Dixon 核，exact-only）

- 常数：π₃ = B(1/3,1/3) ≈ 5.2999（`pi3`）、U = A/π₃ ≈ 0.6845（`pi3_u`）、A = 2π/√3 ≈ 3.6276（`pi3_a`）。power 是常数系数 q（varpi 惯例）；q=0 拒收（ValueError）。
- 核与基底：K = (1−x³)^{−k/3}（k∈{1,2}），矩 T_k = (1/3)B((k+1)/3, 1−k/3) 按 k mod 3 分三类（递推见 dixon.py docstring）。基底 x^{3m+res}(1−x)(a+bx³) 配对相邻残类成 2 维 span；'>'-倒数方向乘外因子 √3/(2π)=1/A（varpi '/π' 的三次类比）做符号改名。
- 六配置为 (kind,comp)↔(kernel,res,outer) 双射，见上表；span 集 {pi3, pi3_inv, U, A, G3, 1}。
- 参数：m 为真实指数级；n 是方向标志（'>'→0、'<'→1）；cu_val=1 为 '+b' 兜底（t·x^{3m+res}(1−x)·K·outer + b，m 取首个矩比越界者）。
- 搜索：m = 0..512 单轴；耗尽 → '+b' 兜底；仍无 → NoSolution。边比 ~1/m 收敛，~1e-3 以上间隙在预算内；区间内相邻 m 交叠 + 兜底覆盖松界 → 真命题一侧全覆盖。
- 校验：exact_check/dixon.py 按 (kind,comp) 重建基矩（CONFIG 定 kernel/res/outer，n 只作簿记不校验）；恒等 = Moment dict 相等；非负 = poly_nonneg 端点规则（a≥0∧a+b≥0）+ 非零被积守卫；兜底支查 t≥0∧b≥0∧非零。
- 渲染：倒数方向印倒数 LHS（`q − r·π₃^{-1}`、`r·U^{-1} − q`），与站端 varpi '>'/gauss '<' 同例；被积 sympy 直排，兜底为 `∫…dx+b_val` 尾挂。
