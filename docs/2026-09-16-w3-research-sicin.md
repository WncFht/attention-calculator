# W3 调研：Si(q) 与 Cin(q) 能否成为新 type 的目标常数

> **落地状态（2026-09-16）**：已落地（`kernels/sicin.py`）——`si_q` 与 `cin_q` 各成一型注册进 `EXACT_TYPES`（四维 span、三次 `poly_nonneg` 判据）；实现笔记见 `docs/2026-09-16-sicin-impl-notes.md`。

日期：2026-09-16。范围：只做调研，不改实现。调研对象为有理参数型常数 $\mathrm{Si}(q)=\int_0^q \frac{\sin t}{t}dt$ 与 $\mathrm{Cin}(q)=\int_0^q \frac{1-\cos t}{t}\,dt$（后者即 $\gamma+\ln q-\mathrm{Ci}(q)$，定义见 DLMF §6.2[^dlmf-sici]）。所有数值断言均用 `.venv/bin/python` 在 `mp.dps = 150` 下复核，参考值取自收敛幂级数而非 `mp.quad`（quad 只有 ~1e-16）。注意 mpmath 的坑：`mp.workdps = N` 是无效赋值——`workdps` 是方法不是属性，静默不改精度；正确写法是 `mp.dps = 150`（限 `from mpmath import mp` 拿到 context 对象的形态——`import mpmath as mp` 拿到模块时此句同样静默无效，须 `mpmath.mp.dps = 150`）。

## 结论

**Si(q) 与 Cin(q) 都可行**，各成一个 type（或一族两核），但有三条硬结构事实决定了形态：

1. **span 是四维的**：$\{C_q,\ \sin q,\ \cos q,\ 1\}$（$C_q$ = si_q 或 cin_q），故 P 必须是三次（四系数），不是 trig_q 族的二次；`poly_nonneg` 目前只覆盖 ≤2 次，需要新写三次精确判据。
2. **Si/Cin 只从 $T_0$/$U_0$ 进入矩**：分部积分在 $k\ge1$ 上自我循环、永远落不回 Si，因此任意基矩里 si_q 的系数恒等于多项式因子在 0 点的值 $F(0)=P(0)=a_0$。这条"钉死"同时解释了 '<' 方向为什么必须换互补核、以及 Si/Cin 为什么不可能共核。
3. **基不能用 $x^m(1-x)^n$**：$m\ge1$ 时 si_q 行全零、方程组恒奇异（实测 121 组 (m,n) 中 110 组奇异，命中只能出现在 m=0 一列，即退化成 1-D 族）；改用 $B_{m,n}=(1+x)^m(1-x)^n$ 后恢复二维搜索，且对大 q 是实质性的（cin(10) '<' 可证边从 1-D 族的 8.13 提到 3.10）。

此外每个方向都存在**无穷 Taylor 阶梯核族**：sin 与 1−cos 的交错截断是全局单侧界（由 1−cos t ≥ 0 逐次积分即得），每个级别 t 给一枚新核、矩公式只在同一 span 内多添纯有理项。q=1 处 t=2 级在 $(m,n)\le10$ 预算内把可证边推到与真值 12 位有效数字重合。大 q 是诚实的短板：q=20 时 '>' 方向在预算内证不出任何正界、'<' 方向只有弱界（详见 §覆盖实测）。

**Si 与 Cin 不能共核**：混合核 $a\sin(qx)/x+b(1-\cos qx)/x$ 的 si、cin 系数分别为 $aF(0)$、$bF(0)$，要消去其一必令 $F(0)=0$ 而两者同灭。落地形态应为"一族两核、按目标常数选核"——与 varpi/gauss 按方向选核同构。**Ci(q) 则结构性不可行**：Ci(q)=γ+ln q−Cin(q)，框架只证"C ⋚ 有理数"，而 Ci(q)⋚r 等价于 Cin(q)⋚γ+ln q−r，右端非有理；把 γ、ln q 加进符号表也不行——它们只在 $U_0$（即 j=0 列）出现且与 Cin 锁成固定比例，{ci_q, euler, ln_q} 三行永远线性相关，方程组恒奇异。不建议单设 ci_q 型。

## 判据回顾

新 type 需要：(a) 区间上定号且零点有限的核 K(x)；(b) 矩 $\int x^j(1+x)^m(1-x)^n K\,dx$ 落在有限维 ℚ-span{常数符号…} 内；(c) 该 span 含目标常数。`engine.solve_moment` 要求方阵（符号数 = P 系数数），`mn_order` 枚举 (m,n)，命中后 `exact_check` 用 combine==target 的字典相等加 `poly_nonneg` 复核。

## 矩结构

记（指标比 trig_q 的 `sin_moments` 平移 1：$T_k=S_{k-1}$，$k\ge1$）

$$T_k=\int_0^1 x^{k-1}\sin(qx)\,dx,\qquad C_l=\int_0^1 x^{l-1}\cos(qx)\,dx,\qquad U_l=\int_0^1 x^{l-1}(1-\cos qx)\,dx.$$

基例与递推（q∈ℚ，全部精确成立、已数值核验至 <1e-120）：

$$T_0=\mathrm{Si}(q),\quad T_1=\frac{1-\cos q}{q},\quad T_2=\frac{\sin q}{q^2}-\frac{\cos q}{q},\quad T_k=-\frac{\cos q}{q}+\frac{(k-1)\sin q}{q^2}-\frac{(k-1)(k-2)}{q^2}T_{k-2}\ (k\ge3);$$

$$C_1=\frac{\sin q}{q},\qquad C_l=\frac{\sin q}{q}-\frac{l-1}{q}T_{l-1}\ (l\ge2);$$

$$U_0=\mathrm{Cin}(q)=\gamma+\ln q-\mathrm{Ci}(q),\qquad U_l=\frac1l-C_l\ (l\ge1).$$

$T_k$ 落在 span{si_q（仅 k=0）, sin_q, cos_q, 1}，$U_l$ 落在 span{cin_q（仅 l=0）, sin_q, cos_q, 1}。闭式也存在（核验过，备用）：

$$T_k=(k-1)!\Big[q^{-k}\sin\tfrac{k\pi}{2}-\sum_{j=0}^{k-1}\frac{q^{j-k}}{j!}\sin\big(q+\tfrac{(k-j)\pi}{2}\big)\Big],\quad C_l=(l-1)!\Big[q^{-l}\cos\tfrac{l\pi}{2}-\sum_{j=0}^{l-1}\frac{q^{j-l}}{j!}\cos\big(q+\tfrac{(l-j)\pi}{2}\big)\Big].$$

组合矩（裸核 '>’方向；'<’与阶梯形只差纯有理加项，见下节）

$$M_j(m,n)=\int_0^1 x^j(1+x)^m(1-x)^n K\,dx=\sum_{l=0}^m\sum_{i=0}^n\tbinom{m}{l}(-1)^i\tbinom{n}{i}\,(\text{T 或 U})_{l+i+j}.$$

三个结构性事实：

- **$F(0)$ 钉死**：si_q/cin_q 只在 $l=i=j=0$ 项出现、系数为 1，故组合矩里特殊常数的系数恒为 $\delta_{j0}$，进而 $\int F K$ 的 si_q 系数 $=a_0=P(0)=F(0)$。'>' 证明要求系数 +1 即 $P(0)=1$，与 P≥0 自洽；'<' 证明要求 −1 即 $P(0)=-1$，在裸核 K≥0 下被积函数必在 0 附近取负——裸核 '<' 方向原理性不可证，必须换互补核（与 varpi '<' 需要伴随核同理）。
- **$m\ge1$ 的 $x^m$ 基恒奇异**：$B=x^m(1-x)^n$ 时 si_q 仅在 $m=i=j=0$ 出现，m≥1 则 si_q 行全零。换成 $(1+x)^m$ 后 l=0 项永远在、si_q 行恒为 $(1,0,0,0)$，非退化恢复。
- **(m,n)=(0,0) 已非退化**：{1, sin_q, cos_q} 三行在 j=1,2,3 列的行列式 $=1/q^4$（六个有理 q 上 Fraction 精确验证），配 si_q 行 $(1,0,0,0)$ 得满秩——(0,0) 处即可解，实测 Si(1)>9/10 与 Cin(1)<1/4 都在 (0,0) 命中。

## 方向裂解与 Taylor 阶梯核族

由 $1-\cos u\ge0$ 逐次积分得全局交错界（$u\ge0$ 严格成立）：$\sin u\le u,\ \ge u-\frac{u^3}{6},\ \le u-\frac{u^3}{6}+\frac{u^5}{120},\ \ge\cdots$；$1-\cos u\le\frac{u^2}{2},\ \ge\frac{u^2}{2}-\frac{u^4}{24},\ \le\cdots$。记下侧截断 $T_{4t+3}(u)=\sum_{r=0}^{2t+1}(-1)^r\frac{u^{2r+1}}{(2r+1)!}\le\sin u$、上侧截断 $T_{4t+1}(u)\ge\sin u$；cos 侧下截断 $S_{2t+2}(u)=\sum_{r=1}^{2t+2}(-1)^{r+1}\frac{u^{2r}}{(2r)!}\le1-\cos u$、上截断 $S_{2t+1}(u)\ge1-\cos u$。四个核族（$t=0,1,2,\ldots$，均在 $(0,1]$ 严格正、零点仅 $x=0$）：

$$K^{>}_{si,t}=\frac{\sin qx-T_{4t+3}(qx)}{x},\quad K^{<}_{si,t}=\frac{T_{4t+1}(qx)-\sin qx}{x},\quad K^{>}_{cin,t}=\frac{(1-\cos qx)-S_{2t+2}(qx)}{x},\quad K^{<}_{cin,t}=\frac{S_{2t+1}(qx)-(1-\cos qx)}{x}.$$

t=0 即 $\big(\sin qx-qx+\frac{q^3x^3}{6}\big)/x$、$\big(qx-\sin qx\big)/x$、$\big(1-\cos qx-\frac{q^2x^2}{2}+\frac{q^4x^4}{24}\big)/x$、$\big(\frac{q^2x^2}{2}-(1-\cos qx)\big)/x$。矩公式 = 特殊组合 ∓ 纯有理部，仍在原四维 span 内（$\beta(m,n,s)=\sum_{l,i}\binom ml(-1)^i\binom ni\frac{1}{s+l+i+1}\in\mathbb Q$ 为 $\int x^s(1+x)^m(1-x)^n dx$）：

$$M^{>,si}_{j,t}=\sum_{l,i}c_{li}T_{l+i+j}-\sum_{r=0}^{2t+1}(-1)^r\frac{q^{2r+1}}{(2r+1)!}\beta(m,n,j{+}2r),\qquad M^{<,si}_{j,t}=\sum_{r=0}^{2t}(-1)^r\frac{q^{2r+1}}{(2r+1)!}\beta(m,n,j{+}2r)-\sum_{l,i}c_{li}T_{l+i+j},$$

$$M^{>,cin}_{j,t}=\sum_{l,i}c_{li}U_{l+i+j}-\sum_{r=1}^{2t+2}(-1)^{r+1}\frac{q^{2r}}{(2r)!}\beta(m,n,j{+}2r{-}1),\qquad M^{<,cin}_{j,t}=\sum_{r=1}^{2t+1}(-1)^{r+1}\frac{q^{2r}}{(2r)!}\beta(m,n,j{+}2r{-}1)-\sum_{l,i}c_{li}U_{l+i+j},$$

其中 $c_{li}=\binom ml(-1)^i\binom ni$。si_q/cin_q 的系数仍是 $\pm\delta_{j0}$（'>' 取 +、'<' 取 −），'<' 侧 $P(0)=+1$ 与 P≥0 相容——互补核把整个方向救活。

**裸核的地位**：si '>' 的裸核 $\sin(qx)/x$ 只在 $0<q\le\pi$ 定号，且在 q=1 的可证边 0.86300 显著差于阶梯 t=0 的 0.94608，还会在大 q 产假证明（见下），**不建议进核族**；cin '>' 裸核 $(1-\cos qx)/x$ 全 q 定号但同样被 t=0 全面压制（0.23662 vs 0.23981）。四族建议统一用阶梯形，裸核仅作参考。

## 覆盖实测

$(m,n)\le10$ 预算内扫"可证边"（'>' 取最大可证 r、'<' 取最小可证 r；目标方程 {C:±1, 1:∓r} 精确解出后数值求 P 在 [0,1] 的允许 r 区间）。真值：Si(1)=0.9460830703672、Cin(1)=0.2398117420006、Si(7/2)=1.8331、Cin(10)=2.9253、Si(20)=1.5482。

| 常数与方向 | 核/级别 | 可证边 | 命中 (m,n) | 与真值差距 |
|---|---|---|---|---|
| Si(1) '>' | 裸核 | 0.8630001924837 | (0,10) | 9% |
| Si(1) '>' | t=0 | 0.9460827 | — | 4e-7 |
| Si(1) '>' | t=1 | 0.9460830704 | — | 3e-11 |
| Si(1) '>' | t=2 | 0.946083070367 | — | 饱和（12 位） |
| Si(1) '<' | t=0 | 0.9462348711304 | (0,10) | 1.5e-4 |
| Si(1) '<' | t=1 | ≈真值 | — | ~1e-9 |
| Si(1) '<' | t=2 | 0.9460830703672 | — | 饱和 |
| Cin(1) '>' | 裸核 | 0.2366210499891 | (0,10) | 1.3% |
| Cin(1) '>' | t=0 | 0.2398117234226 | (0,10) | 1.9e-8 |
| Cin(1) '>' | t=1 | 0.2398117420005 | (0,10) | ~1e-13 |
| Cin(1) '<' | t=0 | 0.2398193274811 | (0,10) | 7.6e-6 |
| Cin(1) '<' | t=1 | 0.2398117420413 | (0,10) | 4e-11 |
| Cin(1) '<' | t=2 | 饱和 | — | — |
| Si(7/2) '>' | t=0 | 1.8329 | — | 2e-4（q>π 仍可证） |
| Cin(10) '>' | 裸/t=0/t=1 | 2.7687 / 2.8387 / 2.9176 | — | 单调逼近 2.9253 |
| Cin(10) '<' | t=0 / t=1 | 3.104 / 2.9614 | — | 逼近中 |
| Si(20) '>' | t=1 / t=3 / t=5 | −5614 / −3048 / −59 | — | 预算内证不出 |
| Si(20) '<' | t=0 / t=1 / t=2 | 8.30 / 1884 / 7175 | (2,0) 等 | t 升高反而劣化 |

两条非对称规律值得写明：'>' 阶梯把核质量压成 Taylor 尾（矩里特殊常数列相对变大、P 退化近常数），t 升高单调变锐，但大 q 下尾项系数 $\sim q^{4t+3}/(4t+3)!$ 要 $t\gtrsim q/(2\pi)$ 量级才压得住，q=20 在预算内失败；'<' 互补核的积分质量 $\approx T_{4t+1}(q)-\mathrm{Si}(q)\sim q$（常数本身只有 ~π/2），证明要从大量里抠小差，且 t 升高注入 $q^{4t+1}$ 级有理压舱物使 P 系数爆炸——大 q 下 '<' 宁可用 t=0 加高 (m,n)。q=20 级别的深振荡区是本族的真实短板。

端到端样例（Fraction 精确求解 → 独立 quad 复核，误差 ≤1e-152）：

$$\mathrm{Si}(1)-\tfrac9{10}=\int_0^1(1+48x-32x^2+16x^3)\,\frac{\sin x-x+\frac{x^3}{6}}{x}\,dx,\qquad \tfrac14-\mathrm{Cin}(1)=\int_0^1\frac{\frac{x^2}{2}-(1-\cos x)}{x}\,dx,$$

$$\tfrac{19}{10}-\mathrm{Si}\!\left(\tfrac72\right)=\int_0^1(1+x)(1-x)^2\Big(1-\tfrac{2856709}{1554672}x-\tfrac{28829}{63456}x^2+\tfrac{98785}{63456}x^3\Big)\Big(\tfrac72-\tfrac{\sin\frac{7x}{2}}{x}\Big)dx.$$

第三条即 q>π 的 '<' 实例（Pmin=0.022>0）。

## Si/Cin 共核与 Ci(q) 的判定

**共核**：设 $K=a\frac{\sin qx}{x}+b\frac{1-\cos qx}{x}$，其组合矩的 si_q、cin_q 系数分别为 $a\delta_{j0}$、$b\delta_{j0}$，即 $aF(0)$、$bF(0)$。要证 Si 目标须 cin_q 系数为 0 ⟹ $F(0)=0$ ⟹ si_q 系数同灭。除非 b=0（回到纯 si 核），不存在同时服务两目标的单核。落地即"一族两核"。

**Ci(q)**：Ci(q)=γ+ln q−Cin(q)。三条路都堵死：(i) 目标语法 "C ⋚ 有理数 r" 下 Ci(q)⋚r ⟺ Cin(q)⋚γ+ln q−r，右端不是有理数（q 有理时 ln q 只有 q=1 为有理即 0）；(ii) 把 {ci_q, euler, ln_q} 列进符号表，则 $U_0=\{$euler$:1,$ ln_q$:1,$ ci_q$:-1\}$——三者只在 l=0（即 j=0 列）出现、三行严格线性相关，方程组恒奇异，能解出的永远是锁定组合 Cin=γ+ln q−Ci 本身；(iii) q=1 时 Ci(1)=γ−Cin(1)，仍含 γ 非有理。结论：不单设 ci_q 型；若站点确有 Ci 展示需求，只能印成 Cin 的组合式，不是新 type。

## 工程备注

- **三次非负判据**：`poly_nonneg`/`poly_nonpos` 现仅支持 ≤2 次（端点 + 判别式）。本族必须扩到三次：建议精确做法为 Sturm 序列在 (0,1) 内实根计数（端点值 ≥0 且内部实根处值 ≥0 的等价判定），或先用 Bernstein 系数全正作充分条件快筛、失败再走 Sturm。这是本族落地的唯一新判据工程量。
- **参数位**：P 四系数，params 需新增 `du_val`（现 au/bu/cu 之外）。
- **矩表复用**：$T_k=S_{k-1}$（trig_q 的 `sin_moments` 平移一位），$U_l=1/l-C_l$ 由 T 表生成，β 纯有理——无需新矩机器，只是同机多符号。
- **搜索维**：(m, n, t)，t 建议集合 {0,1,2}（小 q 已饱和）；WrongDirection/defer 机制沿用。
- **q 域**：四阶梯族全 q>0 合法；若保留 si '>' 裸核则 0<q≤π 必须在**输入校验**强制——q>π 时搜索仍会产出"恒等式精确成立但定号断言为假"的假证明（已实测 q=10 命中 (10,0)），靠搜索侧拦截不干净。q<0 经奇偶性归约：Si 奇（翻方向）、Cin 偶（保方向），同 sin_q 族处理；q=0 退化拒收（递推含 1/q）。

## 建议核形

两型共用一台机器，域 x∈[0,1]，基 $(1+x)^m(1-x)^n$，P 一律三次 $a+bx+cx^2+dx^3$，搜索 (m,n,t) 上限建议 10 与 {0,1,2}。

| 型（建议名） | 常数 | 方向 | 核 K(x)，t=0,1,2,… | 定号性 | q 合法域 |
|---|---|---|---|---|---|
| `si_q` | Si(q) | '>' | $x^{-1}\big[\sin qx-\sum_{r=0}^{2t+1}(-1)^r\frac{(qx)^{2r+1}}{(2r+1)!}\big]$ | (0,1] 严格正，x=0 唯一零点 | q>0（裸核 sin(qx)/x 限 0<q≤π 且不建议用） |
| `si_q` | Si(q) | '<' | $x^{-1}\big[\sum_{r=0}^{2t}(-1)^r\frac{(qx)^{2r+1}}{(2r+1)!}-\sin qx\big]$ | 同上 | q>0 |
| `cin_q` | Cin(q)=γ+ln q−Ci(q) | '>' | $x^{-1}\big[(1-\cos qx)-\sum_{r=1}^{2t+2}(-1)^{r+1}\frac{(qx)^{2r}}{(2r)!}\big]$ | 同上 | q>0 |
| `cin_q` | Cin(q) | '<' | $x^{-1}\big[\sum_{r=1}^{2t+1}(-1)^{r+1}\frac{(qx)^{2r}}{(2r)!}-(1-\cos qx)\big]$ | 同上 | q>0 |

矩公式与符号表：si 侧 $M_j$ 见上节，span $\{$`si_q`, `sin_q`, `cos_q`, `1`$\}$；cin 侧 span $\{$`cin_q`, `sin_q`, `cos_q`, `1`$\}$；目标向量 {'>'}: $\{$C_q$:1,$ 1$:-r\}$，'<' 反号。递推基例 $T_0=$si_q、$T_1=(1-\cos q)/q$、$T_2=\sin q/q^2-\cos q/q$、$T_k$ 二阶递推；$U_0=$cin_q、$U_l=1/l-C_l$、$C_l=\sin q/q-\frac{l-1}{q}T_{l-1}$；$\beta$ 纯有理。q<0 经 Si 奇 / Cin 偶归约，q=0 拒收。

## 数值核验记录

核验脚本 `/tmp/w3_sicin_verify.py`、`/tmp/w3_sicin_comp.py`、`/tmp/w3_sicin_solve2.py`（mpmath dps=150，`.venv/bin/python`，参考值用收敛幂级数），累计 351+ 断言全部 <1e-120：

- $T_0=\mathrm{Si}(q)$ 对 mp.si、$U_0=\mathrm{Cin}(q)=\gamma+\ln q-\mathrm{Ci}(q)$ 对 mp.ci 换算与直接 quad：q∈{2/5,1,3,7/2,5,20} 全过；Si 奇性、Cin 偶性核验过。
- $T_1,T_2,C_1,U_1$ 基例与 $T_k$（k=3..10）、$C_l,U_l$（l=2..10）递推、$T_k,C_l$ 闭式：六 q 值全过。
- 组合矩公式（$(1+x)^m(1-x)^n x^j$ 对 T/U 表）与互补核、阶梯核的矩公式：多组 (m,n,j,q) 对全被积函数级数核验 <1e-80。
- (0,0) 处 {1,sin,cos} 块行列式 $=1/q^4$：六个有理 q 上 Fraction 精确相等。
- 四族阶梯核定号性：q 至 20、t 至 5，网格最小值 >0（如 si '>' q=20 t=5：3.6e-58；cin '>' q=1 t=0：4.3e-20）。
- $x^m$ 基 m≥1 奇异性：121 组 (m,n) 中 110 组 gauss_solve 报奇异、无一命中。
- 端到端：三条样例恒等式 Fraction 解 + 独立 quad 复核误差 ≤1e-152；覆盖表全部边值由精确解 + 1500 点数值 Pmin 扫描给出。
- 假证明实测：si '>' 裸核 q=10 在 (m,n)=(10,0) 给出恒等式精确成立、P≥0 但被积函数变号的"证明"，确认 q≤π 必须进输入校验。

### 参考文献

[^dlmf-sici]: NIST Digital Library of Mathematical Functions. Sine and Cosine Integrals, §6.2（Si(z)、Cin(z)=γ+ln z−Ci(z) 定义与级数）. [dlmf.nist.gov/6.2](https://dlmf.nist.gov/6.2)
