# Lean4 接入调研（2026-09-16）

计划文档里"远期可出 Lean 可检形式（非承诺）"的落地调研。结论先行：**架构上用 verified-checker + 数据证书（路线 B），数学上的真工作量集中在"逐族矩引理"，MVP 选 `e` 或 `pi` 族端到端走通**。本文 API 断言均已对 mathlib4 master 源码快照（`/tmp/mathlib4-src`，toolchain `leanprover/lean4:v4.35.0-rc1`）grep 核实。

## 我们的证书在 Lean 里是什么

一条标准证书（`certificate.build`）断言两件事：

- `integrand == target`：`∫₀¹ x^m(1−x)^n·P(x)·K(x) dx = s·(power·C − bound)` 的 ℚ 字典相等；
- `nonneg`：P 在积分域定号（poly_nonneg / cubic QQ(√D) / Sturm / Bernstein，依族而定）且被积函数非恒零。

对应 Lean 命题：`(power·C : ℝ) ⋚ bound`，证明 = `∫f = u−v`（恒等式）+ `∫f > 0`（非负 + 非恒零 → `intervalIntegral.integral_pos`[^intpos]）。

## 架构：三选一的明确答案

社区硬数据（LRAT-Catcher[^lrat-catcher]、mathlib `lrat_proof`[^lratproof]、equational_theories[^etp]、busybeaver-lean4[^bb]）一致指向：

**路线 B —— Lean 内 checker + 一次性 soundness 定理 + 证书纯数据**

```lean
def check (c : Cert) : Bool := ...            -- 复刻 exact_check 的 ℚ 重算
theorem check_sound (c : Cert) : check c = true → statement c := ...
```

每条证书只剩 `theorem pi_lt_22_7 : Real.pi < 22/7 := check_sound cert_xxx (by native_decide)`——定理语句是真数学命题，证明是一行。对照：

- **路线 A（每条证书生成 tactic 证明/proof term）被实测否决**：mathlib `lrat_proof` 在 63MB 证书上 OOM 到 95.7GB（≈证书 1500 倍），proof term 深度还会撞 kernel 栈（`norm_num` 素数证明 >25bit 就 stack overflow）[^lratproof][^normnum]。我们单条证书虽小（几十~几百有理数），但 3000+ 条的项级复杂度总量仍不建议走这条。
- **折中路线（proof-producing elaborator）**：`norm_num`/`linear_combination` 式，零额外公理但证书大就回到项膨胀问题。可作后期"零公理通道"备选。
- **路线 B 的代价**：soundness 定理是真证明工作量（一次到位后证书随便跑）；checker 若要 kernel 可归约（`decide` 通道）须避开 `partial`/wf-递归/HashMap——LRAT-Catcher 为 kernel 模式专门做了可归约替身[^lrat-catcher]。

### 求值引擎三档（v4.15 起统一后端）

| 方式 | TCB | 适用 |
|---|---|---|
| `by decide`（elab+kernel 两遍归约） | 仅 kernel | 小检查 |
| `by decide +kernel`（kernel 单遍） | 仅 kernel | 中等；ETP `decide!` 同思路[^etp] |
| `by native_decide` | kernel+编译器 | 大检查；v4.29 起每次调用一条独立可审计公理 `foo._native...ax_N`[^v429] |

建议形态：**每个生成的 .lean 文件一次 `native_decide` 批量判整批证书**（`certs.all check`），摊销 per-call 编译开销；要零公理就 `decide +kernel` 通道留给小证书，checker 用 Int/Nat 打包（通分避 `mkRat`/`gcd`，kernel 对 Nat 有 GMP 加速；大表仿 ETP `MemoFinOp` 打 Nat 字面量[^memofinop]）。

### 数据导入

主流做法（ETP/`Generated/` 约定[^etp]）：Python 生成 `Generated/CertsNNN.lean` 数据文件——`def certN : Cert := ⟨...⟩`，有理数用 `mkRat num den`（smart constructor，勿用裸结构字面量）；`.gitattributes` 标 generated。备选 `include_str` + Lean 内 JSON 解析（`Lean.Json.parse` 纯函数可用；坑：lake 不追踪 include_str 依赖，kernel 模式下 parser 结果须落成显式数据[^lrat-catcher]）。**checker 数据层可以不 import mathlib**（`Rat`/`Array`/`Json` 全在 core），只有 statement/soundness 需要 mathlib——这让证书数据模块编译极快。

## mathlib API 映射（已核实）

### 积分侧

- 命题写法：`theorem T : (u:ℝ) - v = ∫ x in (0:ℝ)..1, f x`。
- FTC-2：`intervalIntegral.integral_eq_sub_of_hasDerivAt`（uIcc 全可导）/ `integral_deriv_eq_sub'` / **端点奇异走 `integral_eq_sub_of_hasDerivAt_of_tendsto`**（Ioo 可导 + 端点 Tendsto，ln^r x、(1−x⁴)^(−k/4) 这类核用它）。
- IBP：`integral_mul_deriv_eq_deriv_mul` 系列；换元：`integral_comp_mul_deriv` + 仿射一批；线性：`integral_add`/`integral_const_mul`。
- 现成反导数表（`SpecialFunctions/Integrals/Basic.lean`，全是 `@[simp]` 形）：`integral_pow`、`integral_exp`、`integral_exp_mul_complex`（ℂ 参数一条覆盖 e^{qx} 与 sin/cos(qx)）、`integral_log`（0 端 improper 已内建）、`integral_inv_one_add_sq`（= arctan b − arctan a，**pi 核基例现成**）、`integral_sin_pow`/`integral_sin_pow_mul_cos_pow` 降幂系列、`integral_log_sin_zero_pi_div_two`。
- 严格正性：`intervalIntegral.integral_pos`（`a<b` + `ContinuousOn (Icc)` + Ioc 非负 + ∃内点严格正）[^intpos]；端点奇异的用 `integral_pos_iff_support_of_nonneg_ae` 或往内小区间 `integral_mono_interval`。
- 可积性：`ContinuousOn.intervalIntegrable`；log 端点奇异有 `intervalIntegrable_log'` 先例。
- 惊喜：`bound` tactic（`Mathlib/Tactic/Bound.lean`，aesop 封装的不等式递归证明）与 `TrapezoidalRule`（verified quadrature 带误差界）——将来可做证书内嵌数值复核。

### 非负侧（P ⋚ 0 on 区间）

- `positivity` 比预期强：`evalSub` 会扫局部假设——`x ≤ 1 ⊢ 0 ≤ 1−x`，于是 `x(1−x) ≥ 0`、各结构性核（`exp`、`√`、`rpow`、`1/(1+x²)` 的分母正）都能自动收[^pos]。**但混合系数 P（如 4x³−3x+1）会直接失败且无证书通道**——需要显式见证。
- **主路线：Bernstein 证书**。mathlib 有 `bernsteinPolynomial n ν = C(n,ν)X^ν(1−X)^{n−ν}`（`RingTheory.Polynomial.Bernstein`，ℚ 上线性无关已证）与 `bernstein_nonneg`（`SpecialFunctions/Bernstein.lean:71`，配 `@[positivity]` 扩展）。缺口只有换基等式 `P = Σ b_ν • B_ν`——30–60 行公共引理覆盖全部实例；然后 `Finset.sum_nonneg` + 逐项 `positivity`（`Rat.cast` 扩展把 ℚ 系数非负升到 ℝ）。我们 `euler_gamma.py` 的 `_bernstein` 就是同一套算术，Lean checker 可直接复刻。
- **备用：区间 SOS**。`linear_combination` 现已支持 `≤`/`<` 目标且非负侧条件由 `positivity` 自动消——`P = A + x·B + (1−x)·C`（Markov–Lukács 充要形）可一条 tactic 消费[^lc]；NSPI 项目已这么干。社区 sostactic（cvxpy SDP 后端）可离线产证书再烘成 `linear_combination` 调用。
- **兜底：Sturm**。mathlib 无 Sturm/Positivstellensatz；社区有 sorry-free 单文件实现（Zenodo，~1220 行）可 vendor[^sturm]。我们的三次 QQ(√D) 判据和四次 Sturm 判据都能被 Bernstein/升次覆盖，大概率用不到。
- **[0,π]/[0,π/2] 域**（trig_pi/degree/e_pi 族）：换元 x=πt 回 [0,1] 用同一套，或 Schmüdgen 形 `P = A + x·B + (π−x)·C` 直接 `linear_combination`（`π−x ≥ 0` 由 `x ≤ π` + `evalSub` 收）。

### 常数

| 已在 mathlib | 需就地定义 |
|---|---|
| `Real.pi`（π）、`Real.exp 1`（e）、`Real.log`、`arcsin/arccos/arctan`、`sinh/cosh/tanh`、`arsinh/artanh/arcosh`、`riemannZeta`（实值取 `.re`）、`Real.Gamma`、`Real.eulerMascheroniConstant`、`digamma`、`Complex.betaIntegral`（`=∫₀¹x^{u−1}(1−x)^{v−1}`，直接覆盖 x^m(1−x)^n 形）、`Real.agm`（G = (agm 1 √2)⁻¹，收敛已证）、`bernoulli` | `Catalan := ∑' (-1)^n/(2n+1)²` 或 `∫ arctan x / x`；`Li₂ := ∑' x^k/k²`；`Si/Cin := ∫₀^q sinc/(1−cos)/t`（`Real.sinc` 已有）；`ψ′ := deriv digamma` 或级数；`erf/erfi := 2/√π ∫₀^x e^{∓t²}`（`integral_gaussian` 基建在）；`Dawson`；`ϖ := Γ(1/4)²/(2√(2π))` 或直接 `∫ dx/√(1−x⁴)`；`arccot := π/2 − arctan` |

定义工作量小，**但注意**：statement 侧符号（如 ζ(3)=`riemannZeta 3`）和矩引理产出的常数之间要接上——比如 quadlog 族 `∫x^{2k+1}ln²x/(1+x²) = −(3/16)·(3/4)·2ζ(3) + rat`，需要 `∫ln²/(1+x²) = (3/2)ζ(3)` 这类"基例 = 具名常数"的引理，级数展开 + 逐项积分，是各 Tier-2 族里最实的数学工作。

## 逐族矩引理清单与难度

每条发射证明的恒等式部分归结为**基例 + 递推**：例如 pi 族 `J_k + J_{k−2} = ∫x^{k−2}(1+x²)/(1+x²) = 1/(k−1)` 一步代数，基例 `J_0 = π/4`（`integral_inv_one_add_sq` + `arctan 1`）、`J_1 = ln2/2`（换元）。估算：

- **Tier 1（基建全在 mathlib，纯代数递推）**：e/e_q、e_pi（[0,π] sin 基）、pi、arctan_q/arccot_q、ln_q、artanh/arcoth（归约 ln）、sin_q/cos_q/tan/cot、sinh/cosh/tanh/coth、golden、trig_pi 四型（积化和差展开）。
- **Tier 2（ln^r 递推 + 基例要接具名常数，或端点奇异）**：pi_n/catalan/zeta3/zeta5-11/beta4-10（quadlog，基例是 β(r+1)/η(r+1) 级数恒等式）、ln_q_square/cube/quad（log_family 核 + 幂递推）、arcsin/arsinh（IBP 递推 M_k∝M_{k−2}）、gaussint/dawson/erfi（e^{±q²x²} 两步 IBP，寄生常数就地定义为基积分）、si/cin（Taylor 余量阶梯，Si/Cin 本身即积分定义反而省事）。
- **Tier 3（重活）**：varpi/gauss/pi_sqrt2（(1−x⁴)^(−k/4)，B(1/4,·) 格点 + x=1 端点奇异）、pi3 三型（(1−x³)^(−k/3) Dixon）、gamma（1/ln x 端点奇异 + ∫(1/ln x+1/(1−x))=γ 恒等式 mathlib 多半没有）、li2（级数-积分交换）、psi1（ψ′ 级数形接 `deriv digamma`）。
- **第二证明器**：pade（证书含 l_n 有理函数，Lean 侧验证 `D(ln(1+x)−l_n) = s_n` 是 `field_simp`+`ring` 级恒等 + 有理函数定号——中等）；agm（mathlib 有 `Real.agm`，证书是 isqrt 包络的 ℚ 三明治——中等偏易）；euler_gamma（EM 余项定号引理 mathlib 没有，**最重**，可留到最后或用 Bernstein 重证路线）；composite DAG（sqrt_mul/sqrt_div/pos_transfer 规则各一条小引理 + 递归——容易）。

## 工程

- 本机现状：无 lean/elan/lake；磁盘余 23G。elan toolchain ≈2.8GB + mathlib 预编译 cache ≈7GB，够但紧（cache 长期累积 ~10GB）[^cosmos]。
- 建工程：`lake +leanprover-community/mathlib4:lean-toolchain new <name> math`；**pin mathlib tag 且 lean-toolchain 必须与该 rev 一致**（错位 = cache 全 miss + 几百个库内错）[^wiki]。`lake exe cache get` 拉 ~6800 个 `.ltar`，好带宽 <4min。
- CI：`leanprover/lean-action@v1`（pin SHA），自动跑 cache get；`.lake` ~4-7GB 注意 GitHub cache 配额；`leanchecker`（v4.28 起随工具链）可做环境重放审计[^leanchecker]。
- 多 session 注意：`Generated/` 机生文件 + 生成脚本入 `tools/`、`README` 写复现方法（ETP 约定）；与主仓 Python 包平级建 `lean/` 子目录工程。

## 建议分阶段

- **W8a（MVP，约 1 周 Lean 工作量）**：`e` 族端到端——`Cert`/`check`/`check_sound` 骨架 + I_k 递推（`integral_mul_deriv_eq_deriv_mul` + `integral_exp`）+ 线性 P 端点非负（`positivity` 直收）+ `integral_pos` 收尾 + Python 端 `tools/cert_to_lean.py` 生成器；产物是几十条 golden e 记录编成 mathlib 定理。**验证点**：`lake build` 全过、`#print axioms` 只有预期公理、随机抽定理人工对照。
- **W8b**：Bernstein 非负公共引理 + pi 族（quadlog 的 r=0 支路，绕开级数基例）——把 site 29 型里"基建现成"的族先铺满。
- **W8c+**：quadlog 级数基例（β/η 接具名常数）、各 Tier-2/3 族、第二证明器（pade→agm→composite→euler_gamma 按难度排）。

## 价值与边界（诚实版）

- 换来的：信任从自家 exact_check 换成 Lean kernel（+compiler，若 native_decide）；能抓规范级 bug（常数约定、域边界、证书 schema 歧义）——因为命题在 Lean 里是重新陈述的；产物是数学上可引用的形式化定理库。
- 换不来的：**只管每条证明的 soundness，不管搜索完备性**；也不自动覆盖 site-parity 语义（那是另一层）。
- 未证实项：mathlib 对 `riemannZeta 3` 实值的现成引理面、`integral_pos` 在端点奇异核上的顺路程度、`native_decide` 批量摊销的具体数字——MVP 阶段顺手验证。

## 参考文献

[^intpos]: mathlib4. `Mathlib/MeasureTheory/Integral/IntervalIntegral/Basic.lean:1399` `integral_pos`；FTC/IBP 同目录 FundThmCalculus.lean / IntegrationByParts.lean（本机快照 /tmp/mathlib4-src）。
[^pos]: mathlib4. `Mathlib/Tactic/Positivity/Core.lean`（compareHyp 假设匹配）与 `Basic.lean`（`@[positivity]` 扩展表）；Bernstein: `Mathlib/Analysis/SpecialFunctions/Bernstein.lean`、`Mathlib/RingTheory/Polynomial/Bernstein.lean`。
[^lc]: mathlib4. `Mathlib/Tactic/LinearCombination.lean`（`≤`/`<` 支持 PR #16841）；SOS 先例 github.com/mmaaz-git/sostactic、github.com/ruobingzuo66/NSPI。
[^sturm]: Carles Marín Muñoz. Sturm's theorem in Lean 4. Zenodo 10.5281/zenodo.20707348, 2026.
[^lrat-catcher]: LRAT-Catcher: Importing SAT Solver Certificates into Lean 4 by Reflection. arXiv 2607.00815；repo github.com/leansolving/lrat-catcher（native/kernel 双模式、每 leaf 一模块并行、字符串嵌 statement 模式）。
[^lratproof]: mathlib4. `Mathlib.Tactic.Sat.FromLRAT` 文档及 PR #222（proof-term 重放的内存极限记录）。
[^etp]: teorth. equational_theories（`Generated/` 约定、`decideFin!`/`decide!`、Facts 批量、公理白名单）。github.com/teorth/equational_theories.
[^memofinop]: equational_theories. `MemoFinOp.lean`（elab 期预计算 + Nat 字面量打包表）。
[^bb]: mfornet. busybeaver（BB(5) Lean4 移植：native_decide 默认 + kernel-only 证书通道、证书产物 gitignore 只提交生成器）. github.com/mfornet/busybeaver.
[^v429]: Lean 4.29.0 release notes — One Axiom per Native Computation（RFC lean4#12216 / PR #12217）。
[^wiki]: leanprover-community. Using mathlib4 as a dependency（mathlib4 Wiki）；lake 文档 lean4.dev.
[^cosmos]: emilyriehl. infinity-cosmos#206（v4.32.2 实测磁盘/时间数据）。
[^leanchecker]: leanprover. leanchecker（原 lean4checker，v4.28.0 起入工具链）；lean-action GitHub workflow。
[^normnum]: mathlib4. `Mathlib/Tactic/NormNum/Prime.lean`（proof term 深度 kernel 栈溢出注释）。
