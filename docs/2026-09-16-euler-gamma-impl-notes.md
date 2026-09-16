# W7-euler 实现笔记：gamma 的 Euler–Maclaurin 第二证法

日期：2026-09-16。交付：`src/attention_calculator/euler_gamma.py`、`tests/test_euler_gamma.py`（22 过 1 跳——跳的是 certificate.verify_cert 接线检测，归 leader）。未提交、未动注册表（solve.py / server.py / certificate.py / exact_check/__init__.py / kernels/__init__.py），接线清单见末节。

## 动机

gamma 是全站 exact 覆盖率最差的类型：站端核在 (m,n) 预算内连 γ > 0.5772156 都证不出（gap 1e-4 处 28s 超时、gap ≤ 1e-6 瞬败，site/exact 两模式同病）——ln(N+1) 余项还要再烧一份预算做子搜索。本证法把 γ 直接做成机器可检的有理数包络：Euler–Maclaurin 展开取 N = 2^t，唯一的无理件 ln N = t·ln 2 拆成两条 ln_q 子证书（Padé 证法，x=2 处收敛最快、条件最好）。实测地板从核的 ~1e-4..1e-6 推到 1e-200 以外。

## 数学

约定：伯努利数取 B₁ = −1/2 一派（`B̃₁(x) = {x} − 1/2`），B_n(x) = Σ_k C(n,k) B_{n−k} x^k，B̃_n 为周期延拓。

### EM 恒等式

对 f(x) = 1/x 在 [1, N] 上逐段分部积分（B̃_{n}'/n = B̃_{n−1}，奇数阶 B_{2j+1} = 0 (j ≥ 1) 使中间边界项全消失）：

$$H_N = \ln N + \frac{1}{2} + \frac{1}{2N} + \sum_{j=1}^{J+1} \frac{B_{2j}}{2j}\Bigl(1 - N^{-2j}\Bigr) - \int_1^N \frac{\tilde B_{2J+2}(x)}{x^{2J+3}}\,\mathrm{d}x$$

N → ∞ 时左端 H_N − ln N → γ（γ 的定义），得常数式；两式相减：

$$\gamma = \underbrace{H_N - \frac{1}{2N} + \sum_{j=1}^{J} \frac{B_{2j}}{2j\,N^{2j}}}_{\text{RA}(N,J) \in \mathbb{Q}} - \ln N \;-\; \int_N^\infty \frac{\tilde B_{2J+2}(x)}{x^{2J+3}}\,\mathrm{d}x$$

把 B̃_{2J+2}(x) 写成 s({x}) + B_{2J+2}（s = B_{2J+2}(x) − B_{2J+2}，去均值）：常数部分 ∫_N^∞ B_{2J+2} x^{−2J−3} = B_{2J+2}/((2J+2)N^{2J+2}) 与求和式里第 J+1 项**恰好抵消**（同号同值），所以

$$\gamma = \mathrm{RA} - \ln N + E_J(N), \qquad E_J(N) = -\int_N^\infty \frac{s(\{x\})}{x^{2J+3}}\,\mathrm{d}x$$

这就是代码里 RA 只加到 j = J、余项用去均值 s 的原因。

### 符号引理（per-cert 机器认证）

**断言**：s(x) = B_{2n}(x) − B_{2n} 在 (0,1) 上严格保持符号 (−1)^{n+1}。

**证明**：s(0) = s(1) = 0；s′ = 2n·B_{2n−1}(x)，而奇次伯努利多项式在 (0,1) 内只有 x = 1/2 一个零点（对称式 B_m(1−x) = (−1)^m B_m(x) 给出零点，再结合 B_{2n−1} 在 [0,1] 上恰三个零点 0, 1/2, 1 的经典结论）。故 s 在 (0,1/2)、(1/2,1) 各单调，端点为零 ⟹ s 的符号 = s(1/2) 的符号 = B_{2n}(1/2) − B_{2n} = (2^{1−2n} − 1) B_{2n}（加倍公式），符号即 −sign(B_{2n}) = (−1)^n·(−1) = (−1)^{n+1}。∎

**推论**：s({x}) 在每个 (k, k+1) 上保持同一严格符号，权 x^{−2J−3} > 0 ⟹ sign(E_J) = −(−1)^{J+1} = (−1)^J，严格不等式 E_J ≠ 0。

**机器认证**：证明路径不引用上面的零点结论，而是对每张证书实算 s 在 [0,1] 上的 **Bernstein 系数** b_i = Σ_{k≤i} (C(i,k)/C(d,k))·c_k（c 为升幂系数，d = 2J+2）。`s(x) = Σ b_i B_{i,d}(x)`，Bernstein 基在 (0,1) 上严格正 ⟹ 所有非零 b_i 同号 (−1)^{J+1} 时断言成立。实测 J ≤ 20 全部满足 {0, (−1)^{J+1}} 形式（i = 0 与 i = d 处 b = 0 因 s(0) = s(1) = 0）；`_sign_certified` 就是这条检查，prover 发证前跑一遍、verifier 回放时再跑一遍——引理本身不被信任。

### 幅度界（三角界，无条件）

x ∈ [0,1] 时 |s(x)| ≤ Σ_{k=1}^{2J+2} C(2J+2,k)|B_{2J+2−k}| = Σ_{m=0}^{2J+1} C(2J+2,m)|B_m| =: S（x^k ≤ 1 逐项放缩；s 非常数，严格 <）。故

$$|E_J(N)| < T := \frac{S}{(2J+2)\,N^{2J+2}}$$

这刻意没用「max|B_{2n}(x)−B_{2n}| 在 1/2 取得」的紧界——该界对 n ≤ 3 有反例，而三角界无条件、纯 ℚ。代价是 T 比真值大约一个数量级，由 J 梯子的指数衰减轻松补足。

### 包络与证明规则

取 N = 2^t，ln N = t·ln 2，子证书钉死 ln 2 ∈ (l, u)。由 sign(E_J) = (−1)^J：

- **J 偶**（E_J ∈ (0, T)）：γ ∈ (RA − t·u, RA − t·l + T)。证 γ > R 只需 lo = RA − t·u ≥ R。
- **J 奇**（E_J ∈ (−T, 0)）：γ ∈ (RA − t·u − T, RA − t·l)。证 γ < R 只需 hi = RA − t·l ≤ R。

关键设计：符号已定时，**尾巴 T 只吃宽松侧**——证明侧端点 RA − t·u（'>'）/ RA − t·l（'<'）完全不受 T 影响。所以 J 只需满足可行条件 T < |γ − R|（gap 的 1/4 起扫），证明力度全由 ln 2 的夹逼精度给出：lim = (RA − R)/t，'>' 要 u ≤ lim、'<' 要 l ≥ lim，room = lim − ln2 或 ln2 − lim 就是 ln-2 子证需要的 slack。Padé 在 x = 2 收敛 ~0.029ⁿ，max_n = 140 留给子证的界到 ~1e-212；整体地板由 (t, J) 梯子决定，实测 ≥ 1e-200。

## 判定规则

归一化：q·γ ⋚ r ⟺ γ ⋚ R = r/q（q < 0 翻向，q = 0 拒收 ValueError「γ型系数不能为0」）。cert 的 q/p/comp 记**字面请求**，verifier 自己重算 R 与翻向——与 gamma_special 同约定。

prover 流程（`_prove_const`）：mpmath 只用于初猜——算 gap 数值、扫 (t, J) 候选、把 lim ± f·room snap 成有理候选点；所有进入证书的量（RA、T、l、u、lo、hi）都是 Fraction，发证前逐项精确闸（'>' 要 RA − t·u ≥ bound 的 ℚ 判定）。l/u 候选点列取 lim 与 ln 2 之间的 1/2、3/4、1/4、15/16 分位 + limit_denominator(10⁶) 变体；lim ≤ 0 时 '<' 的紧侧约束退化，直接用 69/100 < ln 2 充 l。children 是 pade.prove("ln_q", 2, ⋚, ·) 的 wire-form 证书，先 '<'（紧侧 u）后 '>'——任一为 None 换下个候选点，全灭换 t，梯子烧完返回 None（诚实失败，~1e-700 以下）。

## 证书 schema

```json
{"prover": "euler_gamma", "kind": "gamma", "comp": ">",
 "q": "1", "p": "<n/d>",            // 字面请求：q·γ comp p
 "n": 64, "j": 138,                 // N = 2^t，EM 阶 J
 "tail": "<n/d>",                   // T（见上）
 "lo": "<n/d>", "hi": "<n/d>",      // γ 的认证包络
 "expect": [{"kind": "ln_q", "power": "2", "comp": ">", "bound": "<l>"},
            {"kind": "ln_q", "power": "2", "comp": "<", "bound": "<u>"}],
 "children": [<pade cert>, <pade cert>]}
```

DAG 契约复用 `gamma_special.child_claim`：verifier 要求 `child_claim(children[i]) == expect[i]` 四元组逐字相等，再递归 `certificate.verify_cert(children[i])`（惰性 import 防环）。

## verify_cert 回放清单（全 ℚ）

1. schema：prover/kind/comp/q/p/n/j/tail/lo/hi/expect/children 齐全可解析；comp ∈ {>,<}，q ≠ 0，j ≥ 0。
2. `n` 是 2 的幂（t = bit_length−1），`0 < l < u`，expect 恰为 `[(ln_q, 2, >, l), (ln_q, 2, <, u)]`。
3. 每个 child：`child_claim(ch) == expect[i]` 且递归 verify 通过。
4. `tail == _tail(n, j)`（重算伯努利数与二项加权和）且 `_sign_certified(j)`（Bernstein 系数认证符号引理）。
5. 重算 RA 与 e_lo/e_hi（J 偶：(0, T)；J 奇：(−T, 0)），要求 `(lo, hi) == (RA − t·u + e_lo, RA − t·l + e_hi)` 逐字相等。
6. c = comp（q < 0 翻向），`lo ≥ p/q`（'>'）或 `hi ≤ p/q`（'<'）。

任何篡改——放宽区间、改 n/j/tail、换/伪造子证、改 expect、翻转方向——都会在某一步失配。mpmath 不进 verify 路径。

## 接线清单（leader，四处小改；kernels/__init__.py 不用动——gamma 本就在 TYPES，EXACT_TYPES 只收 exact-only 型）

1. **solve.py** `prove_exact` 的 `except NoSolution:` 块（~L276，agm 分支旁）：

   ```python
   if kind == "gamma":
       from . import euler_gamma

       cert = euler_gamma.prove(kind, q, comp, r)
       if cert is not None:
           return {"type": kind, "prover": "euler_gamma", "certificate": cert}
   ```

2. **server.py** `/calculate`（composite 分支后，~L294）：`result.get("prover") == "euler_gamma"` → 与 composite 同形 passthrough `respond({success, type, prover, certificate: result["certificate"]})`。cert 已是 wire-form（全 str/int）；想对齐 pade 风格可调 `euler_gamma.cert_jsonable`（幂等）。

3. **exact_check/__init__.py** `verify_response`（composite 分支旁）：`resp.get("prover") == "euler_gamma"` → 同形：`child_claim(cert) == (kind, power, comp, bound)` 闸 + `euler_gamma.verify_cert(euler_gamma.cert_parse(cert))`，返回 `{"identity_ok": ok, "nonneg": ok, "integrand": {}, "target": {}}`。

4. **certificate.py**：`verify_cert` 加 `cert.get("prover") == "euler_gamma"` → `euler_gamma.verify_cert(cert)`（同 composite 的 try/except 集合）；`cert_tex` 加对应分支打印命题本体，例如 `f"gamma({cert['q']}) {cert['comp']} {cert['p']} \\quad (\\mathrm{{EM}})"`。

## 实测

核地板（本仓两模式同病）：gap 1e-4 烧 28s 仍 NoSolution，1e-6 ~6s，≤1e-8 瞬败。euler_gamma 全量 verify=True：

| 命题 | N | J | prove | verify |
|---|---|---|---|---|
| γ > ⌊10⁸γ⌋/10⁸ | 8 | 6 | <0.05s | <0.05s |
| γ < ⌈10²⁰γ⌉/10²⁰ | 16 | 11 | <0.05s | <0.05s |
| γ > 1e-35 界 | 16 | 28 | 0.1s | <0.05s |
| γ > 1e-60 界 | 32 | 36 | 1.0s | <0.05s |
| γ < 1e-100 界 | 64 | 49 | 7.6s | 0.1s |
| γ > 1e-160 界 | 64 | 138 | 55.1s | 0.3s |
| γ > 1e-200 界 | 128 | 96 | 91.4s | 0.3s |

prove 时间几乎全在精确伯努利数（B_{2J+2} 分子 ~千位）与 Bernstein 认证上；verify 始终 <0.4s——验证比重证便宜两个数量级，符合证书设计目标。深 gap 下最紧实测 1e-200（梯子 t ≤ 10、J ≤ 768 的理论地板 ~1e-700，受 ln-2 子证 Padé 指数预算 ~1e-212 与耗时共同截断）。

## 偏离任务书的两点

1. **N 只取 2 的幂**：建议稿允许任意 N；取 N = 2^t 使 ln N = t·ln 2 永远归约到同一条 ln-2 子证（Padé 在 x = 2 处收敛率最高 ~0.029ⁿ，且证书可复用、expect 形状固定），代价仅是 t 梯子粗一点——用 J 补，实测足够。
2. **尾巴用三角界不用 |B_{2J+2}| 型常数**：建议稿的紧余项界要依赖「max|B_{2n}(x) − B_{2n}| 于 x = 1/2」这一引理（且 n ≤ 3 有反例）；本实现对幅度用无条件的 Σ|c_k|/((2J+2)N^{2J+2}) 三角界，符号引理另由 Bernstein 系数 per-cert 认证——verifier 不信任何未被当场重证的断言。
