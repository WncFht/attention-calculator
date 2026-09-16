# mode=exact 阶段状态汇总（2026-09-16）

方案见 `2026-09-16-math-correctness-plan.md`（W0–W5）；本文件是当日晚间的实际落地快照。site 路径冻结于 `v1.0.0-site-parity`，本文不涉及它。

## W0–W5 落地情况

- **W0 验证层**：`exact_check/` 包按族实现 `check(kind, power, comp, bound, params) -> {identity_ok, nonneg, integrand, target}`，`verify`/`verify_response` 两条入口；既是 mode=exact 的发射自检（失败→InternalError），也是判官裁决器与站端输出的复核器。checker 契约两条硬规则（Moment dict 滤零值键、nonneg 含 `bool(integrand)` 非空守卫）已写进 kernel-spec「exact_check 复核器契约」节——四个 agent 独立踩过同一个坑。
- **W1 失真修复**：trig_pi exact 走 `basis_moment` 真矩（绕开 (1,8) 损坏存式）；beta '<' 删转置档、`EXACT_LT_LIMIT=256` 诚实 NoSolution 兜底；power≠1 令发射参数满足印刷 LHS；ln_q_square q∈{5,7} 崩溃完成行列式级根因（det(1,0)∝−(c−4)(c+1)/72c⁶ 等，预言 q=13 同型崩溃未探）；方向判定换 `certified_cmp`（mpmath 80→240→800→2400dps + 护栏带 `2^(30−dps)·max(1,|c|)`）。
- **W2 正确性判官**：`bench/judge_correct.py` + `cases_correct.py`（commit 89c221b）。ground truth 用 `constant_mpf` 区间比较独立裁决，不信仰站端标签；emitted proof 必须过 exact_check（零容忍假恒等式）；site 对照分歧按已知失真簇归因，`unattributed:*` 是要追的。
- **W4 Padé 第二证明器**：`pade.py` 作 ln_q/arctan_q 的在线兜底（(m,n) 搜索耗尽后，预算 MAX_N=50，commits 3fbd546/371f087）。响应形 `{"success","type","prover":"pade","certificate"}`——**无 parameters**，证书即证明；渲染与判官不得假设 parameters 存在。
- **W5 证书**：`certificate.py` + `tools/verify_cert.py`（schema：`docs/2026-09-16-certificate-spec.md`，含 Padé 变体与 gamma_special 的 proof-DAG 变体）。

## exact-only 型清单（EXACT_TYPES 现 26 型）

| 型 | 落地 | 出处 |
|---|---|---|
| zeta5 / zeta7 / zeta9 / zeta11 | 矩核（quadlog η 支路 + 1−2^{1−s} 换算） | kernels/zeta_odd.py |
| beta4 / beta6 / beta8 / beta10 | 矩核（ln_moment 偶支路，β 直入无换算） | kernels/beta_even.py |
| ln_q_cube | 矩核（4 维 span，三次 P + QQ(√D) 判据） | kernels/ln_pow.py，推导 `2026-09-16-ln-cube-derivation.md` |
| ln_q_quad | 矩核（5 维 span，四次 P + Sturm 奇根计数判据） | kernels/ln_pow.py，推导 `2026-09-16-ln-quad-impl-notes.md` |
| arcsin_q / arsinh_q | 矩核（根号核 + 寄生常数） | kernels/arcsin.py / invhyp.py |
| gaussint_q / dawson_q / erfiint_q | 矩核（erf 缩放常数三连，commit 80e6b74） | kernels/gauss_erf.py |
| pi_sqrt2 | 矩核（lemniscate 余元常数 π√2，两核两基） | kernels/pi_sqrt2.py |
| pi3 / pi3_u / pi3_a | 矩核（Dixon/B(1/3,1/3) 格点六格双射） | kernels/dixon.py |
| li2_q | 矩核（3 维 span + 位移见证族，n×(fam,d) 双轴） | kernels/li2.py |
| psi1_q | 矩核（望远镜核对，一维 m 扫描） | kernels/trigamma.py |
| si_q / cin_q | 矩核（4 维 span，Taylor 余量阶梯 t∈{0,1,2}） | kernels/sicin.py |
| gamma14 / gamma34 / gamma12 | 复合命题型（证书 DAG） | kernels/gamma_special.py |

逐型规格已并入 `kernel-spec.md`「exact-only 型清单」节。已否决：erf 本体（√π 障碍）、Γ(1/4)/Γ(1/3) 字面幂（格点奇偶锁）、lnA/Glaisher（矩空间符号数恒超 P 系数）、arcosh_q（span 无 "1" 方向）。逐篇结论见 `2026-09-16-w3-research-*.md` 头部状态行；实现笔记见 `2026-09-16-*-impl-notes.md`。命名注意：erf 调研文里的 `expint_q` 落地名为 **erfiint_q**。

## 判官结果与修复史

8830 条全量扫描（golden 输入 + 对抗语料）快照：proved 3178、rejected-false 4295、unsolved-true 733、rejected 476、equal-claim 86、**BUG ×62**。两类 BUG 均已修复并复测转阴：

- `BUG:wd-on-true` ×41（hyperbolic q<0 未做奇偶归约）→ commit 460c5b1 奇偶归约到 |q|；
- `BUG:crash` ×21（power=0 触发 Fraction(1,0)）→ 同 commit 退化域拒（e_q/cosh_q q=0、coth_q q≤0 等报 ValueError）。

修复后 a9 复扫 9594 例（含 exact-only 型）：**零 BUG**。覆盖率地板前三弱曾是 varpi 26/113、gauss 28/105、gamma 32/109——已由 W7 第二证明器波次补齐（varpi/gauss 走 `agm.py` 区间包络，gamma 走 `euler_gamma.py` EM 包络；三者深界不再依赖 (m,n) 预算）。W7 落地后 a9 对这三型同输入重跑 660 条的 before→after：**varpi 26→108（+82）、gauss 28→105（+77）、gamma 32→105（+73），零回退、零 BUG/FLAG、零 truth_mismatch**，全部系 unsolved-true→proved；36 条 gauss/varpi 转置簇（真命题、站端假证明）无一漏全证出。残余 unsolved-true：gauss 0 条；gamma 4 条全是 power=0 退化点（euler_gamma 按设计拒"命题退化为有理数比较"，口径上与 varpi q=0 的 rejected-input 分类不统一，待对齐）；varpi 4 条为 ~9e15 分母级的最深连分数界，在 pi 子证地板之外。复测残余 `FLAG:rejected-true-claim` 个位数，全部是退化域拒的真命题（`e_q 0 > 1/2` 即 `1>1/2` 之类）——域校验按设计拒收，属判官口径问题而非求解器缺陷。

## 在途项

- W7 第二证明器波（已收官）：`w7-agm`（AGM 区间证法 commit d5e3000——gauss 双向 ~1e-2000 间隙、varpi ~1e-16 受 pi oracle 地板限）、`w7-euler`（Euler–Maclaurin γ 证法 commit 2d090e7——EM 展开 N=2^t + 定号有界余项 + 两个 ln_q Padé 子证消去 lnN，Bernstein 系数现场重证符号引理，实测地板 1e-200 双向）、`w7-ln4`（ln_q_quad，见上表）、`w7-decomp-ext`（commit 09225b5——decompose_exact 原子表覆盖全部 26 型 + 11 个原不可拼写的 site 型，修掉 arg<0 被误读为倒数的预存 bug）。三弱型地板全被第二证法接管：varpi/gauss 走 agm，gamma 走 euler_gamma。
- `decompose_exact.py`（commit 4856f27）：可证构造的界分配；**已接线**——`/decompose_inequality` 的 `mode=exact` 走 decompose_exact（commit a492d0d），site 路径不变。
- 已修复的精度陷阱：`decompose_exact` 的 slack 曾在 `workdps` 外与 mpf 字面量相乘塌缩成 float64（commit e00e5cd 改 Fraction 侧乘）；`import mpmath as mp` 下 `mp.dps = N` 是模块属性静默无操作（真精度留在 15），须用 `mp.mp.dps` 或 workdps——测试文件已清查。
- 未探明：ln_q_square q=13 的预言崩溃（行列式级根因已定位 q∈{5,7} 同型）；pi3 '<' 的 unsolved-true 残余是搜索强度问题而非正确性。
