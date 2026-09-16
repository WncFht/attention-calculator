# W3 调研：li2_q（Li₂(q)，q 有理）的矩空间推导与可行性（2026-09-16）

> **落地状态（2026-09-16）**：已落地（`kernels/li2.py`）——`li2_q`（q∈ℚ，q<1，q≠0）已注册进 `EXACT_TYPES`；实现偏离与参数编码见 `docs/2026-09-16-li2-impl-notes.md`。

结论：可行，建议注册为 exact-only 型 `li2_q`，目标常数 $\mathrm{Li}_2(q)=\sum_{i\ge1}q^i/i^2$，合法域 $q\in\mathbb Q$、$q<1$、$q\neq0$。核取 $\varphi(x)=-\ln(1-qx)/x$，质量矩 $\int_0^1\varphi=\mathrm{Li}_2(q)$；矩恰闭于三符号 $\mathrm{span}\{\mathrm{Li}_2(q),\,\ln(1-q),\,1\}$，与 ln 族同构地要求 $P=a+bx+cx^2$ 三系数。本型区别于已有类型的结构性约束是：$\mathrm{Li}_2$ 只出现在零阶矩中，因此基被强制为 $(1-x)^n$（$m\equiv0$）、$P$ 的常数项恒被钉死为 $+1$，剩余 $(b,c)$ 由 $\{\ln,1\}$ 二阶方程组解出。证明方向不靠翻转 $P$ 而由核的符号因子吸收：'>' 用劣超子核 $\varphi-L$，'<' 用优超子核 $L-\varphi$；按 $\varphi$ 的凸性（$q>0$ 严格凸、$q<0$ 凹）四象限各有一族系数落在 $\mathbb Q+\mathbb Q\ln(1-q)$ 的多项式 $L$。$q>0$ 时截断次数 $d$ 是第二条搜索轴，间隙随 $d$ 几何衰减。全部矩公式经 mpmath dps=100 对照 `mp.polylog(2,·)` 与 `mp.quad` 核验（误差 $\sim10^{-98}$），端到端恒等式核验至 $\sim10^{-102}$。注册侧与 `docs/kernel-spec.md` 的 mode=exact 清单、`docs/2026-09-16-math-correctness-plan.md` 的 W3 条目对应。

## 积分表示与核 φ

$\mathrm{Li}_2(q)=-\int_0^q\frac{\ln(1-t)}{t}dt$（定义级数逐项积分即得[^lewin]；特殊值见 DLMF §25.6[^dlmf-dilog]）。代换 $t=qx$（$q\neq0$，不依赖 $q$ 的符号）得

$$\mathrm{Li}_2(q)=\int_0^1\varphi(x)\,dx,\qquad \varphi(x)=-\frac{\ln(1-qx)}{x}=\sum_{i\ge1}\frac{q^i}{i}\,x^{i-1}.$$

$\varphi$ 在 $[0,1]$ 上连续，$\varphi(0)=q$ 取极限值，$\varphi(1)=-\ln(1-q)$；单调递增（$q>0$ 时由级数逐项显然，$q<0$ 时经网格核验多个 $q$ 值）。定号随 $q$：$q>0$ 时 $\varphi>0$，$q<0$ 时 $\varphi<0$（$\varphi(x)=-\ln(1+|q|x)/x$），故裸核只能承载一个方向；$\mathrm{Li}_2(q)$ 本身随 $q$ 变号，无需奇性归约。凸性是关键性质：$q>0$ 时 $\varphi''>0$ 逐项成立、严格凸；$q<0$ 时 $\varphi$ 凹（$q\in\{-1/2,-2,-9/10\}$ 等网格核验）。凸性决定哪一侧存在多项式界：凸函数的弦在上方、切线在下方，凹函数反之。

## 矩闭合推导

记 $A_k=\int_0^1x^k\varphi\,dx$（$k\ge0$）。$A_0=\mathrm{Li}_2(q)$ 即目标常数。$k\ge1$ 时 $A_k=-\int_0^1x^{k-1}\ln(1-qx)\,dx$；对 $\int x^j\ln(1-qx)\,dx$ 分部积分一次得

$$\int_0^1x^j\ln(1-qx)\,dx=\frac{\ln(1-q)}{j+1}+\frac{q}{j+1}\,J_{j+1},\qquad J_l=\int_0^1\frac{x^l}{1-qx}\,dx.$$

由 $x^l/(1-qx)=-x^{l-1}/q+x^{l-1}/\big(q(1-qx)\big)$ 递推（$J_0=-\ln(1-q)/q$）得闭式 $J_l=-q^{-l-1}\ln(1-q)-\sum_{i=0}^{l-1}q^{i-l}/(i+1)$，代回即得一般矩公式（$k\ge1$）：

$$A_k=\frac{q^{-k}-1}{k}\,\ln(1-q)+\frac{q^{-k}}{k}\sum_{l=1}^{k}\frac{q^l}{l}.$$

要点：$k\ge1$ 的 $A_k$ 只含 $\ln(1-q)$ 一次项与有理数——无 $\ln q$、无 $\pi^2$、无更高阶 polylog；$\mathrm{Li}_2$ 仅活在 $A_0$。因此对任意多项式位移 $L\in\big(\mathbb Q+\mathbb Q\ln(1-q)\big)[x]$，核 $\sigma(\varphi-L)$ 的矩空间恒为三维 $\mathrm{span}\{\mathrm{Li}_2(q),\,\ln(1-q),\,1\}$，位移只改后两列。

基元 $x^{m+j}(1-x)^n$（$j=0,1,2$ 三行对应 $P$ 三系数）对核 $\sigma(\varphi-L)$ 的矩通式为

$$M_{m,n,j}=\sigma\sum_{t=0}^{n}(-1)^t\binom nt\Big[A_{m+j+t}-\sum_{i}\frac{L_i}{m+j+t+i+1}\Big],$$

其中 $L_i$ 为 $L$ 的系数（可带 $\ln(1-q)$ 符号分量），$\{m,n,j\}$ 的依赖全部显式。

## 结构事实：Li₂ 只在零阶矩

$A_k$ 的 $\mathrm{Li}_2$ 分量为 $[k=0]$：$x^{m+j}(1-x)^n$ 展开后零阶项仅存于 $m=j=t=0$。三条推论：其一，$m\equiv0$ 被强制——$m\ge1$ 的基整列无 $\mathrm{Li}_2$ 分量，$x^m$ 不再是搜索轴；其二，$j=0$ 行的 $\mathrm{Li}_2$ 列恰为 $\sigma$，方向方程 $a\sigma=s$（$s=+1$ 为 '>'、$-1$ 为 '<'）与端点条件 $f(0)=a\,\sigma\,(q-L(0))\ge0$ 联立，强制 $\sigma=s$ 且 $a=+1$——四象限一律 $P(0)=1$；其三，$(b,c)$ 由 $j=1,2$ 两行在 $\{\ln,1\}$ 列上的 $2\times2$ 方程组解出，该行列式在全部测试核与 $n\le15$ 上非零（个别奇异 $n$ 跳过即可）。同一"$k=0$ 独占"还给出两条负结论：核乘 $x^t$（$t\ge1$）会彻底杀死 $\mathrm{Li}_2$ 分量——$\varphi$ 中的 $1/x$ 是本型命门，核形没有可乘之余地；把型写成 $[0,q]$ 上核 $-\ln(1-x)/x$ 经 $x=qu$ 化为同族（仅差 $q^{m+n}$ 因子），不构成新型。

方向—核配对随之确定：'>' 要 $+\mathrm{Li}_2$ 质量故 $\sigma=+1$、需要劣超子 $L\le\varphi$；'<' 要 $-\mathrm{Li}_2$ 质量故 $\sigma=-1$、需要优超子 $L\ge\varphi$。是否存在这样的多项式 $L$ 完全由 $\varphi$ 的凸性决定。

## 方向—核配对：四象限核族

- '$>$'、$q>0$（$\varphi$ 凸）：Taylor 截断劣超子 $L=p_d=\sum_{i=1}^{d}q^ix^{i-1}/i$。尾项 $\varphi-p_d=\sum_{i>d}q^ix^{i-1}/i$ 逐项非负（解析成立，非仅数值）；$d=0$ 即裸核 $\varphi$，$d=1$ 为 $\varphi-q$，$d=2$ 即 0 点切线。$d$ 是第二搜索轴，间隙随 $d$ 几何衰减。可再叠 $x^d q^{d+1}/(d+1)$ 项（尾项 $\ge$ 其首项）作微紧化。$x=1$ 处切线系数含 $\ln(1-q)$、仍在 span 内，是合法备选；内点切线会引入 $\ln(1-qx_0)$ 第四符号，不进首批。
- '$<$'、$q>0$（$\varphi$ 凸）：弦优超子 $L=\mathrm{chord}(x)=q+\big(-\ln(1-q)-q\big)x$。凸函数弦在上方，且 chord 是**最优一次优超子**：任一一次 $L\ge\varphi$ 必满足 $L(0)\ge q$、$L(1)\ge-\ln(1-q)$，而 $L-\mathrm{chord}$ 为一次式、两端非负故处处非负，即 $L\ge\mathrm{chord}$ 逐点成立。高次族 $L=p_d+C'_d\,x^d$、$C'_d=q^{d+1}/\big((d+1)(1-q)\big)\in\mathbb Q$：由尾项 $\le x^d q^{d+1}/\big((d+1)(1-q)\big)$ 得（$x^{i-1}\le x^d$、$1/i\le1/(d+1)$ 逐项放缩），符号自由且带 $d$ 轴几何收敛，$d\ge2$ 起实测优于弦；另有尾弦族 $p_d+\big(-\ln(1-q)-H_d\big)x$（凸尾 $\le$ 自身弦，$d\le2$ 时退化为 chord），实测介于两者之间。
- '$>$'、$q<0$（$\varphi$ 凹）：弦劣超子 $L=\mathrm{chord}$——镜像论证给出最优一次劣超子（任一一次劣超子逐点 $\le$ chord）。升级：$L=\mathrm{chord}+t\,x(1-x)$，有理 $0<t\le t^*(q)=\min_{x\in(0,1)}(\varphi-\mathrm{chord})/(x(1-x))$（$q=-1/2$ 时 $t^*\approx0.02240$）。
- '$<$'、$q<0$（$\varphi$ 凹）：切线优超子 $L=\tan0=q+q^2x/2$（凹函数在切线下方；在锚定 $(0,q)$ 的一次优超子中斜率最小）。升级：$L=\tan0+t\,x^2$，$t^*(q)\le t$，$t^*(q)=\sup_x(\varphi-\tan0)/x^2<0$——即向下收紧（$q=-1/2$ 时 $t^*\approx-0.03047$）；$x=1$ 处切线（系数含 $\ln(1-q)$）同样有效，实测略差。

统一图景：四族都是"$\varphi$ 的首一阶多项式逼近"，由 $\mathrm{sign}(q)$ 决定的凸性选择用 0 点切线还是端点弦，由方向决定取上界还是下界；$q>0$ 两方向都有以截断次数 $d$ 为轴的几何收敛族。

## q = 1/2 特检

Euler 恒等式 $\mathrm{Li}_2(1/2)=\pi^2/12-\ln^2\!2/2$[^dlmf-dilog]。$q=1/2$ 时 $\ln(1-q)=-\ln2$，span 为 $\{\mathrm{Li}_2(1/2),\ln2,1\}$，判据 $\mathrm{Li}_2(1/2)\lessgtr r$ 等价于 $\pi^2/12-\ln^2\!2/2\lessgtr r$——pi_n 型（span $\{\pi^2,1\}$）与 ln_q_square 型（span $\{\ln^2q,\ln q,1\}$）各自都无法表达这个混合常数，即使只做 $q=1/2$ 一个值本型也覆盖了新地面。示例（裸核、$n=0$、$r=4/9$ 恰为该点最优界）：

$$\mathrm{Li}_2(1/2)-\frac49=\int_0^1\Big(1-\frac{8x}3+\frac{16x^2}9\Big)\Big(-\frac{\ln(1-x/2)}x\Big)dx\ge0,$$

被积函数逐点非负（$P$ 顶点 $x=3/4$ 处判别式 $4ac-b^2=0$ 贴零），mpmath 核验左右差 $3.6\times10^{-102}$。

## q 合法域与退化

建议注册域 $q\in\mathbb Q$、$q<1$、$q\neq0$。$q<0$ 任意幅度合法（$\mathrm{Li}_2$ 由积分实值定义，四象限核表覆盖）。$q=1$ 时 $\ln$ 列恒为零、span 塌缩为 $\{\pi^2/6,1\}$（$A_k=H_k/k$ 纯有理，$\mathrm{Li}_2(1)=\pi^2/6$）：'>' 仍工作但与 pi_n 冗余，'<' 本质不可能（$\varphi(1)=\infty$，无多项式优超子）；建议排除或文档化为退化别名。$q>1$ 时 $1-qx$ 在 $x=1/q\in(0,1)$ 过零、支点落入域内，排除；$q=0$ 核塌缩为零，排除。实用瓶颈在 $q\to1^-$：'>' 的截断轴需 $d\sim\log(\mathrm{gap})/\log q$（$q=9/10$ 时 $d\ge20$ 才有非平凡界），'<' 各优超子族可行集变稀疏岛。

## 覆盖数值实验

扫描以方向（按上表配对）$\times\,n\times$ 核族参数（$d$ 或 $t$）为轴；gap 指可证最优 $r$ 与 $\mathrm{Li}_2(q)$ 之差。

- '$>$' $q=1/2$、$K=\varphi-p_d$：固定 $d=2$ 时间隙随 $n$ 振荡下降（偶 $n$ 系统性更优），$n=16$ 达 $5.4\times10^{-6}$；固定 $n$ 随 $d$ 几何衰减，$(n,d)=(16,6)$ 达 $6.6\times10^{-11}$、$(8,8)$ 达 $1.4\times10^{-10}$。裸核（$d=0$）仅 $\sim O(\log n/n)$（$n\cdot\mathrm{gap}$ 从 $n=4$ 的 0.22 缓增至 $n=32$ 的 0.39）——截断轴是覆盖的关键。$n=0$ 可解出闭式界 $r^*=\min\!\big(4q^2/(1+q)^2,\;q/(2(1-q))\big)$（$a=1,\,b=-2r(1+q)/q,\,c=4r$；顶点入域 iff $q>1/3$），$q=1/2$ 给 $4/9$，与扫描 $r^*$ 精确一致。
- '$<$' $q>0$：弦核在 $q=1/2$ 达 gap $2.8\times10^{-4}$（$n=8$）、$2.1\times10^{-4}$（$n=10$），下降缓慢；$x^d$-lift 族同 $q$ 下 $(d,n)=(6,0)$ 达 $2.4\times10^{-6}$、$q=2/3$ 下 $(6,0)$ 达 $1.7\times10^{-4}$，为该方向推荐族。注意 '$<$' 可行集可为**有界岛**而非射线（如 $q=2/3,n=2$ 的可行区间约 $[0.84,1.72]$）：求解器逐 $r$ 判定不受影响，但扫描"最优可证界"须向下探岛缘。
- '$>$' $q<0$、$K=\varphi-\mathrm{chord}$：$q=-1/2$ 时 gap $2.8\times10^{-3}$（$n=0$）$\to1.3\times10^{-4}$（$n=12$）；二次劣超子 $\mathrm{chord}+t\,x(1-x)$（$t=0.0224<t^*$）改进至 $4.4\times10^{-4}\to3.0\times10^{-5}$。
- '$<$' $q<0$、$K=\tan0-\varphi$：$q=-1/2$ 时 gap $4.0\times10^{-3}\to2.6\times10^{-5}$（$n=0..12$，约 $1/n^2$）；二次优超子 $\tan0+t\,x^2$（$t$ 取 $t^*\approx-0.0305$ 上侧有理）$n=8$ 达 $1\times10^{-5}$。
- $q=9/10$（近 1 退化示范）：'>' 在 $d\le10$、$n\le11$ 无非平凡界，$d=20$ 起恢复；'<' 弦核 $n\ge5$ 全灭，$x^d$-lift $d=6$ 仅 $n\le3$ 有岛（最优 gap $0.027$）。结论：$1-q$ 不宜太小，或实现上让 $d$ 上限随 $q$ 放宽。
- 双参数基 $(1-x)^n(1+x)^p$ 扫描全部回退 $p=0$——无增益，已排除。

## P 形与定号规则

$P=a+bx+cx^2$ 且 $a\equiv+1$（结构推论，非额外约束）；非负判据完全沿用站点二次规则——$c\le0$ 查两端点 $a\ge0$、$a+b+c\ge0$，$c>0$ 且顶点 $-b/2c\in(0,1)$ 时查 $4ac-b^2\ge0$（见 `docs/kernel-spec.md` 非负性判据节），本型不需要新判据。仅当引入第四符号（内点切线扩展）时 $P$ 升四次，才需借用 ln_q_cube 的三次判据（`docs/2026-09-16-ln-cube-derivation.md`）。

## 负面结果（已排除方向）

除前文已述（核乘 $x^t$ 杀死目标常数、$[0,q]$ 写法同型、$(1+x)^p$ 无增益、内点切线引入第四符号）外：ln_q 式分母幂 $(1+cx)^s$ 进入核会破坏闭合——部分分式把矩带到平移参数处的 $\mathrm{Li}_2$，越出三符号 span；2 系数 $P=a+bx$ 对三符号超定、generic 无解；'<' 在 $q=1$ 本质不可能（无多项式优超子），属结构性限制而非搜索不足。

## 建议核形

单一类型 `li2_q`（exact-only，站点 29 型中无对应）。搜索轴：comp $\times\,n\ (0\le n\le16)\times$ 核族参数（$d$ 或 $t$；$q$ 接近 1 时放宽 $d$）。

- **K(x)**：$K(x)=s\big(\varphi(x)-L(x)\big)$，$\varphi(x)=-\ln(1-qx)/x$（$x=0$ 处取极限 $q$），$s=+1$（'>'）/ $-1$（'<'）。$L$ 按下表：

| comp | q 范围 | 基础 L | 升级轴 |
|---|---|---|---|
| > | $0<q<1$ | $p_d=\sum_{i=1}^{d}q^ix^{i-1}/i$（$d\ge0$，$d=0$ 为裸核） | 增大 $d$；可叠 $x^dq^{d+1}/(d+1)$ 微紧化 |
| < | $0<q<1$ | $p_d+C'_dx^d$，$C'_d=q^{d+1}/\big((d+1)(1-q)\big)\in\mathbb Q$（$d\ge1$；一次最优简化版 chord $=q+(-\ln(1-q)-q)x$） | 增大 $d$ |
| > | $q<0$ | chord（最优一次劣超子） | $\mathrm{chord}+t\,x(1-x)$，$0<t\le t^*(q)$ |
| < | $q<0$ | $\tan0=q+q^2x/2$ | $\tan0+t\,x^2$，$t^*(q)\le t<0$ |

- **domain**：$x\in[0,1]$；基 $(1-x)^n$（$m\equiv0$ 为结构推论而非搜索轴）；$j=0,1,2$ 三行对应 $P$ 三系数。
- **矩公式**：$M_{n,j}=\sigma\sum_{t=0}^n(-1)^t\binom nt\big[A_{j+t}-\sum_iL_i/(j+t+i+1)\big]$；$A_0=\mathrm{Li}_2(q)$；$k\ge1$：$A_k=\frac{q^{-k}-1}{k}\ln(1-q)+\frac{q^{-k}}{k}\sum_{l=1}^kq^l/l$。
- **目标符号**：$\mathrm{Li}_2(q)$（注册 `constant_mpf = mp.polylog(2, mpf(q))`）。
- **span 符号表**：$\{\mathrm{Li}_2(q),\ \ln(1-q),\ 1\}$，恰三维；$L$ 系数的 $\ln(1-q)$ 分量只改后两列、不扩维。
- **定号性**：各核 $\ge0$ 由解析论证（凸性/弦切关系/尾项放缩）加网格核验成立，$K(0)=0$ 或 $K(1)=0$ 贴端点为常态；$P$ 沿用站点二次非负判据。
- **求解结构**：$\mathrm{Li}_2$ 行钉 $a=+1$；$\{\ln,1\}$ 列 $2\times2$ 解 $(b,c)$，$n\le15$ 实测行列式恒非零。
- **合法 q 域**：$q\in\mathbb Q$、$q<1$、$q\neq0$；建议排除 $q=1$（退化、与 pi_n 冗余且 '<' 本质不可能）。
- **恒等式形态**：'>' 发 $\mathrm{Li}_2(q)-r=\int_0^1(1-x)^n(1+bx+cx^2)\,\big(\varphi-L\big)\,dx\ge0$；'<' 发 $r-\mathrm{Li}_2(q)=\int_0^1(1-x)^n(1+bx+cx^2)\,\big(L-\varphi\big)\,dx\ge0$；被积函数逐点非负，与站点输出形态一致。

## 数值核验记录

核验脚本在 /tmp（li2_moments.py、li2_poly.py、li2_lt.py、li2_final.py 等，`.venv/bin/python`，mpmath dps=60–100）：

- $A_k$ 闭式：$q\in\{1/2,1/3,2/3,-1/2,-2,9/10\}$、$k=0..10$，对照 `mp.quad` 直接积分与级数 $\sum_iq^i/\big(i(k+i)\big)$ 两路，最大误差 $\sim10^{-98}$。
- 核定号：四类核在 2000–4000 点网格 $\min K\ge-5\times10^{-61}$（零为端点贴合）；尾项界 $\mathrm{tail}\le x^dq^{d+1}/\big((d+1)(1-q)\big)$、$\mathrm{tail}\le x\cdot\mathrm{tail}(1)$ 及两侧 $t^*$ 定义式均经网格核验。
- $2\times2$ 行列式：四象限代表核 $n=0..15$ 全部非奇异。
- 端到端（Fraction 精确求解 $\to$ 二次非负判据 $\to$ `mp.quad` 对照 LHS，dps=100）：$\mathrm{Li}_2(1/2)-4/9=\int_0^1(1-\frac83x+\frac{16}9x^2)\varphi\,dx$（裸核 $n=0$，差 $3.6\times10^{-102}$）；$\mathrm{Li}_2(1/2)-29/50$（$\varphi-p_2$、$n=2$、$P=1+\frac1{10}x+\frac45x^2$，差 $2.2\times10^{-102}$）；$59/100-\mathrm{Li}_2(1/2)$（$\mathrm{chord}-\varphi$、$n=2$、$P=1+\frac{439}{130}x-\frac{89}{65}x^2$，差 $4.4\times10^{-102}$）；$\mathrm{Li}_2(-1/2)+9/20$（$\varphi-\mathrm{chord}$、$n=4$，差 $1.7\times10^{-102}$）；$-2/5-\mathrm{Li}_2(-1/2)$（$\tan0-\varphi$、$n=4$，差 $4.5\times10^{-103}$）。
- $n=0$ 闭式界 $\min\big(4q^2/(1+q)^2,\,q/(2(1-q))\big)$ 与扫描 $r^*$ 一致（$q=1/2$ 得 $4/9$，顶点贴零）。

### 参考文献

[^dlmf-dilog]: NIST Digital Library of Mathematical Functions. Dilogarithms, §25.6（积分表示、特殊值 Li₂(1)=π²/6 与 Li₂(1/2)=π²/12−ln²2/2）. [dlmf.nist.gov/25.6](https://dlmf.nist.gov/25.6)
[^lewin]: Lewin. Polylogarithms and Associated Functions. North-Holland 1981（Li₂ 的积分表示与函数方程）.
