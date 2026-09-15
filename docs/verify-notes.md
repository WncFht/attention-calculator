# 验证方法学：bench/verify.py 与 integrand.py

目标：对每条"恒等式声称"（站点 golden 或本方 solver 输出）判定
`claimed_lhs = ∫_a^b f(x)dx > 0` 是否成立、被积函数是否定号。
这是整个复现工程的基准工具：重建错了会在这里暴露，solver 证错了也是。

## 威胁模型

verify 要区分三类情况：

- `valid`：恒等式成立（|∫f − lhs| ≤ 1e-30·max(1,|lhs|)，50dps）且被积函数定号；
- `indefinite-sign`：恒等式数值成立但被积函数在域内变号（证明无效）；
- `false-identity`：积分真值与声称 LHS 不符——本方 solver bug 的证据
  （曾怀疑站点会放假恒等式，后证实是探测字段名笔误导致服务端回退到 pi 型；
  该分类保留用于抓本方错误）。

站点 success=false 的报错记录另外核对 `error_consistent`：
"方向反了"要求真方向确与声称相反；"未找到解"要求真方向与声称一致
（式子为真但预算内证不出）；输入校验类报错不涉及方向，记 None。

## 高精度积分：选了什么、为什么

- **mpmath.quad（tanh-sinh，默认算法），mp.dps=50**。tanh-sinh 天然适应端点
  奇异（1/√(1−x⁴)、ln x 类），从不在端点求值；内部奇异需用 points 切分
  ——本项目的核在域内均无内部奇异（ln 族分母零点 −1/(q−1) 恒在 [0,1] 外），
  无需切分。`error=True` 的误差估计仅作参考（文档明示其本身是启发式，
  有反例显示估计乐观而结果错误）。
- **精度余量**：50dps 下 quad 典型残差 ~1e-45~1e-51，容差取 1e-30 相对值，
  余量约 15~20 个数量级。gamma 型因被积函数结构（见下）实测残差 ~1e-39，
  仍远低于容差；如个别 case 逼近容差，可升 dps 复算。
- **验证收敛性的方法**：同一 case 可在多 dps 下重算比对位数（文档建议
  做法）；本实现对 gamma 已实测 dps 40/50/60/80 结果逐位一致。
- **失败兜底**：`mp_func` 对求值异常（端点舍入到恰好奇异点）向内微扰重试。

## mpmath.iv 与 python-flint（Arb）

- `mpmath.iv` 区间算术对**求值**是严格的（f(v) ⊆ f̂(v)），文档明示可用于
  无理数不等式的机器证明；但 `iv.quad` 复用同一套 tanh-sinh 启发式网格，
  只把舍入误差装进区间——离散化误差没有严格界。**iv.quad 不是认证积分**。
- `python-flint` 的 `flint.acb.integral` 是 Arb 的认证积分（Gaussian 求积
  + 先验误差界，输出含严格半径的复球）。0.9.0 提供 cp310-abi3 与 cp314
  的 macOS arm64 / manylinux x86_64 轮子，本机（3.12）与 devbox（3.14,
  x86_64）均可装。已接入 `--rigorous`：对所有非 valid 记录复核
  `lhs ∈ acb积分区间` 是否成立。实现细节：`flint.ctx.prec=200`
  （~60d，压住 lhs 的 mpf→arb 解析误差），隶属判定用
  `arb.contains(arb(lhs))`（考虑双方半径的集合包含，严格）。
  注意 python-flint 0.9.0 里 `arb.mid()/.rad()` 仍返回 arb 球，
  不能 `float()` 取中值——53bit 精度会淹没 1e-30 量级的判据。
  限制：被积函数需在路径上解析，端点代数奇异靠自适应细分，
  可能较慢。
- 定位：`mp.quad` 是主判定（快、精度足）；`acb.integral` 是对存疑记录的
  终审；`iv` 只用于点值符号判定（目前未用，扫描点值用 mp 足够）。

## 不变号判定

三层，由严到宽：

1. **精确层 `sign_exact`**：被积函数唯一可能变号的因子是待定多项式
   P。各族的 P 都化归 t∈(0,1) 上的 QQ 多项式（t 取 x、x²、x⁴、sin x），
   用 sympy `count_roots`（Sturm 序列）精确数内部根：先剥离端点根
   t、t−1（对应区间端点，允许取零；注意 t−1 在 (0,1) 为负，每剥一次
   要翻转符号），内部无根则 P 定号——这一步是**严格**判定。
2. **扫描层 `sign_scan`**：200 点初始网格 + 端点附近采样 +
   在 |f| 最小处自适应加密三轮（捕捉窄下凹）。覆盖 B·P·K 全体，
   包括 K 因子不显注定号的情形（gamma 的 L/U 核、sin(qx) 若 q>π）。
3. 未覆盖的残余：扫描是启发式，理论上可漏掉极窄变号区间。精确层
   对多数类型已兜底；gamma 第一段核这类超越函数组合只受扫描保护
   （其定号性在数学上有经典结论支撑，见 kernel-spec）。

## gamma 核的数值稳定性（关键实现细节）

L(x)=1/(1−x)+1/ln x−1/2、U(x)=(2−x)/2−1/(1−x)−1/ln x 在 x→1 都是
±∞ 对消出有限值。分立式写法在 tanh-sinh 深节点（x−1 ~ 1e-25）上把
1/(1−x) 与 1/ln x 两个 ~1e25 的量相减，结果完全是舍入噪声；
mpmath 的误差估计仍报告 1e-69，属系统性错误（实测偏差 ~0.013，
恰好吞掉子证明项）。合并成单分式（y=1−x）后所有项都是 O(y) 量级，
对消只损失 ~log(1/y) 位：
U = ((y²+y−2)·ln x − 2y)/(2y·ln x)，L = ((2−y)·ln x + 2y)/(2y·ln x)。
注意保持 (1−x) 字面形式——sympy 展开成 2x−2 会把对消搬回 O(1) 项。

## 参数编码反推记录（/~50 次 API 探测）

- **ln_q / ln_q_square / artanh_q / arcoth_q 的分母幂 s 不回传**。
  多数样本 s=n，但 arcoth_q 3 实测 s=1、n=0 → s 是独立搜索量。
  重建把 s 留作自由符号，verify 枚举 0..15 数值反解。
- **gamma 的核指数 k 即参数 cu_val**（全量 golden 渲染式实证：
  cu_val=4→x⁴、5→x⁵、16→x¹⁶，与早前"ln(k+1) 子证明界"反推一致——
  k+1 正是子证明要证的 ln 参数）。完整结构：
  `f = x^{cu}·K_dir + x^m(1−x)^n(au+bu·x)/(u·(1+cu·x)^e)`，
  分母幂 e：'<' 取 n，'>' 取 cu_val（后者仅 r=57/100 一单样本，
  残余不确定度在此）；u_val=0 时无子项，改为加性常数 a_val。
  子证明分子用 au_val,bu_val/u_val 的**原始整数比**——注意
  a_val≠au/u 普遍成立（如 937/56 vs 2811/168 差了千分位），
  不能用 a_val 代替。verify 的 k 枚举代码留作兜底但不再触发。
- **gauss '>' 在 a_val=0 时切换模板**（仅 r=0 一单样本）：
  `∫₀¹ au·(1−x)·√(1−x⁴)/π dx + b_val`，核是分子 √ 非 1/√，
  且带加性常数项——与标准模板并存，按 a_val==0 判别。
- **站点会返回数值不成立的"恒等式"**（全量 golden 实证 34 例）：
  trig_pi 四型在 m=1,n=8 档共 30 例（渲染方程两侧差 ~1e-6~1e-7，
  如 sin(2π/5) vs 5266/5537 的界证明里积分真值是 lhs 的 ~40 倍），
  gauss '<' 大系数档 3 例（差 ~4000 倍）。
  **trig_pi 30 例的根因已定位**：站点的 (m,n)=(1,8) j=0 存储公式比真值多
  δ(α)·(C−1)，δ(α) 为有理函数（见 kernel-spec.md "trig_pi 四型"节，
  33 点插值恢复并留出验证）。站点 (a,b) 满足的是含 δ 的方程，
  真实积分两边不等，差 a·δ(α)·(1−C)；不等式本身仍成立
  （bound_ok=True、被积函数定号）→ verify 记 verdict=false-identity。
  gauss 3 例根因另查。
- **trig_pi 四型（sin_pi_q、cos_pi_q、sin_q_degree、cos_q_degree）的
  c_val 不是多项式系数**，是核频率：sin 型核 sin((1−2q)x)、
  cos 型核 sin(2qx)。P=a+b·sin x 只用 a,b。
- **artanh_q/arcoth_q 的 c_val 是归约后的 ln 参数 q'**
  （artanh q→ln((1+q)/(1−q))/2；arcoth q→ln((q+1)/(q−1))/2）。
- **pi_n 分数幂**：power=u/v 时站点证 (π^{u/v})^v ⋚ r^v 即 π^u ⋚ r^v，
  核为 ln^{u−1}(1/x)/(1+x²)。基指数 = 2m + [u 偶]（奇偶须与 u 相反：
  偶 u 取奇指数得 η(u)=Q·π^u；奇 u 取偶指数得 β(u)∈Q·π^u）。
- **varpi/gauss 的 n 不进被积函数**：基固定 x^{4m+r}(1−x)，
  r：varpi '<'→3、'>'→1（核含 1/π，LHS 1−r·ϖ⁻¹）；
  gauss '>'→0（含 1/π，LHS G−r）、'<'→2（无 π，LHS r·G⁻¹−1）。
- **tan_q/cot_q/tanh_q/coth_q 外层因子**分别为 1/cos q、1/sin q、
  1/cosh q、1/sinh q，乘在积分号外。
- **pi/e 的 power 是常数倍率**（3π、e/2），核不变。
- m,n 编码各型不同：pi 的 m→x^{2m}，varpi 的 m→x^{4m+r}，
  pi_n 的 m→x^{2m+[u偶]}，zeta3 的 m→x^{2m+1}，e_pi 的 m→sin^m x。

## 残余风险

- 自由参数反解只剩 ln 族的 s（枚举 ≤15，界过紧会误判
  unresolved-param；观测值 ≤3）；gamma 的 k 已改用 cu_val 直读。
  gauss '>'/a_val=0 变体与 gamma '>' 子项分母幂=cu_val 均为单样本
  规则，后续新样本出现偏差时应首先复查这两处。
- 恒等式判定的 1e-30 相对容差对 |lhs| 极小（紧贴真值的界）可能过紧，
  失败时应看 abs_deviation 数量级人工复核。
- 扫描层对极端窄变号可能漏检；精确层 + acb 复核缓解。
- acb.integral 对端点代数奇异可能慢；rigorous 层目前只对非 valid 跑。
