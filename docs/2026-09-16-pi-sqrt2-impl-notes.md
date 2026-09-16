# pi_sqrt2 实现笔记（W3，2026-09-16）

`pi_sqrt2` 是 lemniscate 格点上最后一个干净的单核型：常数 `S = π√2 = Γ(1/4)Γ(3/4)`（余元公式，平凡）。规格取自 [[2026-09-16-w3-research-gamma-quarter]] §"lemniscate 族剩余单核型清单"与§"建议核形"；本文记实现决定与回归数据。交付三个文件：`src/attention_calculator/kernels/pi_sqrt2.py`（核，单文件含两方向）、`src/attention_calculator/exact_check/pi_sqrt2.py`（checker）、`tests/test_exact_pi_sqrt2.py`（38 项，全绿）。

## 核形与矩

两方向各用一条核（与 varpi/gauss 共享基形 `x^{4m+res}(1−x)(a+bx⁴)`、固定 `(1−x)` 因子、`a+bx⁴` 非负规则 `a≥0 ∧ a+b≥0`）：

- '>'：`K(x) = (1−x⁴)^{-3/4}`，res=2。`J'_{4j+2} = u_j·S/4`、`J'_{4j+3} = v_j`，槽位 i 的矩为 `{S: u_{m+i}/4, 1: −v_{m+i}}`，`u_j = (3/4)_j/j!`（递推 `u_{j+1} = u_j(4j+3)/(4j+4)`）、`v_j = (1)_j/(5/4)_j`（`v_{j+1} = v_j(4j+4)/(4j+5)`）。
- '<'：`K(x) = (1−x⁴)^{-1/4}`，res=3。`J''_{4j+3} = w_j`、`J''_{4j+4} = z_j·S`，槽位 i 的矩为 `{1: w_{m+i}, S: −z_{m+i}}`，`w_j = (1/3)(1)_j/(7/4)_j`（`w_{j+1} = w_j(4j+4)/(4j+7)`）、`z_j = (1/16)(5/4)_j/(2)_j`（`z_{j+1} = z_j(4j+5)/(4j+8)`）。

矩公式经 mpmath dps=100 复核（J'' 系残差 < 1e-79，J' 系受核奇性限制在 quad 精度 ~1e-26 内为零；调研文档同表已独立核验）。两方向的 2×2 系统均恒非奇异：'>' 行列式 `= −u_mv_m/(4(4m+4)(4m+5))`，'<' 行列式 `= −3w_mz_m/((4m+7)(4m+8))`，故每个 m 都有唯一解，搜索只由定号性决定。

目标向量直证字面命题、无倒数归一（各核 span 自带 S）：'>' 取 `{S: q, 1: −r}`（证 `q·S − r > 0`），'<' 取 `{1: r, S: −q}`（证 `r − q·S > 0`）。矩符号 `"pi_sqrt2"` 与 `"1"`；证书 `_SYMBOL_TEX` 需登记 `"pi_sqrt2": "\\pi\\sqrt{2}"`。

## 参数编码

沿用 emit 九字段；`n` 是方向标志（'>'→0、'<'→1）而非指数，与 beta_family 同约定。`c_val`/`cu_val` 恒 `"0"`（P 只有两系数），无兜底模板，checker 不读 `cu_val`。被积函数由 `(comp, m, au_val, bu_val, u_val)` 唯一确定：'>' 为 `x^{4m+2}(1−x)(au+bu·x⁴)/(u·(1−x⁴)^{3/4})`，'<' 为 `x^{4m+3}(1−x)(au+bu·x⁴)/(u·(1−x⁴)^{1/4})`。integrand.reconstruct 若登记本型可按此直译；render_equation 已实现同式（sympy 排版：'(1-x^4)^{3/4}' 与 '\sqrt[4]{1-x^4}' 分母）。

## 搜索预算与求解行为

单轴 m 扫描，`LIMIT = 256`——沿用 beta 族 exact '<' 档 `EXACT_LT_LIMIT`（调研文档建议"可沿用 256 一档"）。四个验证算例最深 m=13；'<' 深度约为 '>' 的 3–4 倍（实测界距 0.057→m=13、0.007→m=116、0.027→m=29；'>' 0.003→m=95、0.043→m=5），1e-4 级界（如 '<' 4443/1000）超出预算返回诚实 NoSolution，全程约 0.12s。渐近 `u_m/v_m → 4/S`、`z_m/w_m → 1/S` 保证任意真界在足够 m 处可解；预算只截断深度，不产生错误证明。

power 是纯系数、不进核——`check_input` 不设校验。q≤0 沿 varpi/gamma 约定：假命题由上游 certified_cmp 判 WrongDirection；真命题（如 `−S > −5`、`0 < 1`）在每个 m 上解出符号不定的 P——一致非正解在数学上不可能出现（∫P·K≤0 会证出反向不等式，与真命题矛盾），故耗尽后报 NoSolution。q=0 ∧ r=0 由 certified_cmp 判 EqualClaim。核级（绕过 solve.prove 直调 kernels.pi_sqrt2.prove）假命题同样耗尽为 NoSolution——方向判定是上游职责，核内不重复。

## 验证算例的发射恒等式

- `π√2 > 4` → m=0：`∫₀¹ 4x²(1−x)(1−x⁴)^{-3/4}dx = π√2 − 4`（a=4,b=0）。
- `π√2 > 22/5` → m=5：`∫₀¹ x²²(1−x)(2577291/5992448 + 5499385/749056·x⁴)(1−x⁴)^{-3/4}dx = π√2 − 22/5`，与调研文档逐系数一致。
- `π√2 < 5` → m=0：`∫₀¹ x³(1−x)(13/3 + 56/3·x⁴)(1−x⁴)^{-1/4}dx = 5 − π√2`。
- `π√2 < 9/2` → m=13：`a = 2666335599248386194517/81716445659466301440`、`b = 10995247838419146150971/143003779904066027520`。

测试经 monkeypatch 临时接线 FAMILY/EXACT_TYPES/constant_mpf 后跑通完整 `solve.prove(exact=True) → certificate` 链路（含 `verify_cert` 独立复核）；登记后该接线为无害冗余。

## 待 leader 接线清单

- `solve.FAMILY["pi_sqrt2"] = "pi_sqrt2"`（render/exact_check/cert 三路分发共用此表）。
- `kernels.EXACT_TYPES` 追加 `"pi_sqrt2"`。
- `integrand.constant_mpf` 表追加 `"pi_sqrt2": lambda: q * mp.pi * mp.sqrt(2)`（certified_cmp 方向判定依赖）。
- `certificate._SYMBOL_TEX["pi_sqrt2"] = "\\pi\\sqrt{2}"`。
- `engine.EXPONENT_LIMIT["pi_sqrt2"] = 256`（404 文案"指数不超过256"）。
- `integrand.reconstruct` 可选登记（见「参数编码」节公式）；bench/verify 需要时才加。

## kernel-spec.md 补遗草稿（exact-only 型清单条目）

`pi_sqrt2`（kernels/pi_sqrt2.py）：常数 `S = π√2 = Γ(1/4)Γ(3/4)`（余元公式），lemniscate 格点 `(1−x⁴)^{c−1}` 族仅剩的干净单核型（调研 [[2026-09-16-w3-research-gamma-quarter]] §末两节；Γ(1/4) 本体因奇偶引理不可达，走 W6 复合命题）。'>' 核 `(1−x⁴)^{-3/4}`、基 `x^{4m+2}(1−x)(a+bx⁴)`，矩 `{S: u_{m+i}/4, 1: −v_{m+i}}`（`u_j=(3/4)_j/j!`、`v_j=(1)_j/(5/4)_j`）；'<' 核 `(1−x⁴)^{-1/4}`、基 `x^{4m+3}(1−x)(a+bx⁴)`，矩 `{1: w_{m+i}, S: −z_{m+i}}`（`w_j=(1/3)(1)_j/(7/4)_j`、`z_j=(1/16)(5/4)_j/(2)_j`）。目标 `{S:q, 1:−r}` / `{1:r, S:−q}` 直证 `q·S ⋚ r`，无倒数归一、无兜底模板；`n` 为方向标志（同 beta_family）。power 为纯系数不进核；q≤0 免校验（假命题上游 WrongDirection、真命题诚实 NoSolution，同 varpi/gamma）。单轴 m 搜索 LIMIT=256（沿用 EXACT_LT_LIMIT 档）：验证算例 m≤13，'<' 深度约为 '>' 的 3–4 倍；渐近 `u_m/v_m→4/S`、`z_m/w_m→1/S` 保证真界可达。矩符号 `"pi_sqrt2"`/`"1"`。
