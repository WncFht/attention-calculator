# beta4 / beta6 落地笔记（2026-09-16）

W7 实现笔记。Dirichlet β 的偶数点：β(4) = Σ (−1)^k/(2k+1)⁴ ≈ 0.98894455174110533611（无已知闭式，Catalan β(2)=C 的偶数伴侣），β(6) ≈ 0.99868522221843813544 作为同构附带型一并实现。结构与 zeta_odd 完全镜像，差别只在奇偶支路的选择。

## 文件与类型

- `src/attention_calculator/kernels/beta_even.py` — 一族两型：`beta4`（power·β(4)）、`beta6`（power·β(6)），power 是系数 q（zeta_odd/varpi 约定）。
- `src/attention_calculator/exact_check/beta_even.py` — `check(kind, power, comp, bound, params)`。
- `tests/test_exact_beta_even.py` — 46 通过 + 2 skip（注册后转绿）。

## 矩推导与支路选择

quadlog `ln_moment(k, r, odd)` 的两条支路（docs/kernel-spec.md 闭式）：

- 偶次幂 x^{2k}：`(−1)^k r!·(β(r+1) − Σ_{i<k} (−1)^i/(2i+1)^{r+1})` —— 符号列直接是 β(r+1)。
- 奇次幂 x^{2k+1}：`(−1)^k r!/2^{r+1}·(η(r+1) − …)` —— zeta_odd 走这条，η 乘 1−2^{−r} 归一到 ζ。

β(4) 要求 r+1=4 且落在 β 支路，即 **r=3、偶次幂**（奇次幂 r=3 落 η(4)=7ζ(4)/8 ∈ Q·π⁴，表达不了 β(4)）。同理 β(6) 是 r=5、偶次幂。**无归一化系数**——cc 本身就乘在 β(r+1) 上，factor=1（zeta_odd 的 η→ζ 换算在这里没有对应物）。

核形状：`x^{2m}(1−x²)^n (a+b·x²) ln^r(1/x)/(1+x²)`，r=3（beta4）/ 5（beta6），(m,n) 走 `mn_order(LIMIT)`。r 为奇时 ln^r(1/x) ≥ 0 保持核非负；渲染端必须印 `ln^r(1/x)` 而非 `ln^r(x)`——奇 r 下两者差一个负号，印错会让 "> 0" 尾自相矛盾（zeta_odd 的 r 全偶无此问题）。

### 数值复核（60dps，mpmath 1.3.0）

β(4) = `mp.dirichlet(4, [0,1,0,−1])`、β(6) = `mp.dirichlet(6, [0,1,0,−1])`（Hurwitz 分解 4^{−s}(ζ(s,¼)−ζ(s,¾)) 交叉验证一致）。`ln_moment` 偶支路对数值求积：

| k | r=3: cc | r=3: rat | quad − (cc·β4+rat) | r=5: cc | quad − (cc·β6+rat) |
|---|---------|----------|--------------------|---------|--------------------|
| 0 | 6 | 0 | −7e-59 | 120 | −2e-54 |
| 1 | −6 | 6 | −4e-61 | −120 | 2e-60 |
| 2 | 6 | −160/27 | 6e-61 | 120 | −6e-60 |
| 3 | −6 | 100162/16875 | −4e-61 | — | — |

全在求积精度内，矩式属实。

## 参数编码

复用 emit 九槽：`{m, n, a_val, b_val, c_val, au_val, bu_val, cu_val, u_val, unified_form}`——P=(au+bu·x²)/u 乘 `x^{2m}(1−x²)^n ln^r(1/x)/(1+x²)`，c_val/cu_val 恒 0，u 为系数分母 lcm。目标向量 `{kind: s·power, "1": −s·bound}`（s=+1 为 '>'）。

## LIMIT 与实测深度

LIMIT=10，与 zeta_odd 及全族默认一致。实测覆盖（CF 渐近分数作界）：

| 命题 | err | 解出 (m,n) |
|------|-----|-----------|
| β(4) < 98309/99408 | −5.2e-12 | (9,9) |
| β(4) > 5725/5789 | +1.7e-9 | (6,7) |
| β(4) > 805/814 | +1.1e-6 | (4,4) |
| β(6) < 559814/560551 | −4.1e-13 | (8,9) |
| β(6) > 4048587/4053917 | +2.6e-14 | (9,10) |
| β(6) > 129889/130060 | +1.3e-11 | (7,7) |

深度边界：β(4) > 1873596/1894541（err 7.7e-14）首个解在 (11,11)，恰超预算 → 诚实 NoSolution（测试断言此行为）。'>' 侧渐近分数（下界）比 '<' 侧（上界）同深度压得更紧，与 zeta_odd 观察同型。假命题的行为：中等紧度（1e-5..1e-9）一律 mid-scan WrongDirection；超紧度假命题（'<', 1873596/1894541，err +7.7e-14）在预算内搜不出定号 P，耗尽 NoSolution——同为诚实回答。

## q=0 与退化输入

沿用 zeta_odd 约定：**不做域拒**。power=0 退成纯有理命题（0·β(4) ⋚ r ⇔ 0 ⋚ r），核照常解常数恒等式；0-vs-0 在 exact 管线里被 certified_cmp 的 sign=0 → EqualClaim 先行拦掉，核内则解出 P≡0、checker 的 `bool(integrand)` 非零守卫判 nonneg=False。负 q 无折叠——β(4) 无奇偶对称可言，−q 只是另一系数，目标向量自然携带。

## kernel-spec.md 草稿（「exact-only 型清单」追加）

```markdown
- `beta4` / `beta6`（kernels/beta_even.py）：`power·β(4)`、`power·β(6)`，power 是系数。复用 quadlog ln_moment 的**偶支路** r=3/5：矩 (−1)^k r!·(β(r+1)−partial) 的符号列即 β(r+1) 本身，无归一化系数（与 zeta_odd 的 η→ζ 换算不同——那里 cc 乘的是 η，要补 1−2^{−r}）。核 x^{2m}(1−x²)^n(a+bx²)ln^r(1/x)/(1+x²)；r 为奇，渲染必须印 ln^r(1/x)（ln^r(x) 差负号）。LIMIT=10：'<' 侧证到 ~5e-12（β4 (9,9)），'>' 侧 ~3e-14（β6 (9,10)）；β(4) '>' 7.7e-14 需 (11,11) 超预算。β(8)+ 同构可扩（r=7,9,… 仅加 BETA 表行），未实装。
```

## 注册所需改动（leader 侧，本 session 未动）

- `kernels/__init__.py`：`EXACT_TYPES += ["beta4", "beta6"]`。
- `solve.py`：`FAMILY += {"beta4": "beta_even", "beta6": "beta_even"}`。
- `integrand.py` `constant_mpf`：`"beta4": lambda: q * mp.dirichlet(4, [0, 1, 0, -1])`、`"beta6": lambda: q * mp.dirichlet(6, [0, 1, 0, -1])`（mpmath 1.3.0 实测可用；备选 `mp.mpf(4)**(-s) * (mp.zeta(s, mp.mpf(1)/4) - mp.zeta(s, mp.mpf(3)/4))`，两者 60dps 一致）。
- `certificate._SYMBOL_TEX`：`"beta4": "\\beta(4)"`、`"beta6": "\\beta(6)"`（缺省回退 \mathrm{beta4} 也能印，但 \beta(4) 更规范）。
- 无需 EXPONENT_LIMIT 条目（默认 10 即核内 LIMIT）。
- `docs/kernel-spec.md`「exact-only 型清单」并入上方草稿节；清单计数与「11 型」表述需随登记同步。
