# W7-agm 实现笔记：gauss / varpi 的 AGM 第二证法

日期：2026-09-16。交付：`src/attention_calculator/agm.py`、`tests/test_agm.py`（35 项全过）、`kernels/gamma_special.py` 增量（新规则 `pi_div_agm` + kind 白名单加 `varpi`，属本任务交付的非注册表文件）。未提交注册表改动（solve.py / server.py / certificate.py / exact_check/__init__.py 归 leader 合并，清单见末节）。

## 动机

lemniscate 二型是全站可证界地板最粗的族（实测见末节表）：varpi '>' 连 ϖ > 2.6 都证不出（(m,n) 搜索 + `sqrt_bound_proof` 兜底在 bound ≥ 8/5 处整体失效），gauss '>' 止于 ~0.8253。AGM 把两类地板分别推到 ~1e-4000（受 isqrt 标度限制）与 ~1e-16（受 pi 子证地板限制）。

## 数学

M = AGM(1, √2)：a₀ = √2, b₀ = 1，迭代 a_{k+1} = (a_k+b_k)/2、b_{k+1} = √(a_k·b_k)。常数关系（数值已复核）：

$$G = \frac{1}{M} \approx 0.8346268417, \qquad \varpi = \frac{\pi}{M} = \pi G \approx 2.6220575543$$

**AGM 夹逼（证明依赖的引理）**：a₀ > b₀ > 0 时，对一切有限 n 严格有 b_n < M < a_n。归纳：a_{n+1} − b_{n+1} = (√a_n − √b_n)²/2 > 0（有限步永不收敛，a_n ≠ b_n 保持），且 M ≥ b_{n+1} = √(a_n b_n) > b_n、M ≤ a_{n+1} < a_n。严格性让端点等号也可用：`1/a_hi ≤ R ⟹ G > R`（因为 M < a_hi 严格）。

**isqrt 区间引理**：对 x = p/q ∈ ℚ≥0 与十进制标度 D，`isqrt(D²p/q)` 给出最大 k 使 (k/D)² ≤ x；对称地有最小 k 使 (k/D)² ≥ x。两者都是精确 Fraction。

**b-迭代近似为何严格**：不直接算 √ 的真值，而是对每条真迭代量维护含它的有理区间。a-迭代是区间端点算术均值（精确）；b-迭代取 `[sqrt_floor(a_lo·b_lo, D), sqrt_ceil(a_hi·b_hi, D)]` ⊇ [√(a_k b_k) 的真值范围]。归纳得真 b_n ∈ [b_lo, b_hi]、真 a_n ∈ [a_lo, a_hi]，故 M ∈ [b_n, a_n] ⊆ [b_lo, a_hi] —— 记 (lo, hi) = (b_lo, a_hi)，恒有 lo < M < hi 严格。

**收敛**：AGM 间隙 ~exp(−c·2ⁿ)（实测 n=5 → 2e-43、n=6 → 3e-80），每步 isqrt 漂移 O(1/D)，n 步累计 O(n/D)。所以地板只由标度 D 决定：digits=64 时宽度 ~4·10⁻⁶⁴。

## 判定规则（全 ℚ）

归一化：q·C ⋚ r ⟺ C ⋚ R=r/q（q<0 翻向，q=0 拒收）。

- **gauss**：'>' 证成 ⟺ hi·R ≤ 1（即 1/hi ≥ R，R ≤ 0 时自动成立）；'<' 证成 ⟺ lo·R ≥ 1。驳斥：'>' 在 lo·R ≥ 1（M > 1/R ⟹ G < R）时即停，'<' 对称 hi·R ≤ 1。
- **varpi**：ϖ = π·G 拆成两条子证。'>' 取 B = 1/hi（agm 子证 G > B）、A = R·hi（pi 子证 π > A），恰使 A·B = R；<' 对称 B = 1/lo、A = R·lo。'>' 且 R ≤ 0 时取 A = 1（任意正 A 都让 A·B > 0 ≥ R）。子证严格蕴含乘积严格，故 A·B == R 的边界形即可。

pi 子调用 `solve.prove("pi","1",c,str(A),exact=True)` 作 oracle（gamma_special 同款），证书逐字进 children。witness 选取不需要 mpmath——最优 B 就是包络端点、A = R/B 直接是最松的 pi 命题；整个 prove 路径零 float、零 mpmath。

## 证书 schema

gauss（顶层与 varpi 子证共用；记号对齐 child_claim 回退分支的字段名）：

```json
{"prover": "agm", "kind": "gauss", "power": "1", "comparison": ">",
 "bound": "417/500", "agm_iter": 2, "agm_digits": 64, "lo": "<n/d>", "hi": "<n/d>"}
```

`lo`/`hi` 是 M 的严格包络端点；`agm_iter`=重放步数 n，`agm_digits`=十进制标度（D = 10^digits）。verify 逐项重放 `agm_enclosure(n, 10^digits)` 并比对 (lo, hi) 逐字相等，再复核归一化后的单边判定——篡改 n/digits/端点/方向/界任一项都失配。n ≤ 64、digits ≤ 4096 双帽挡 DoS 形伪造。

varpi（gamma_special DAG 复用，prover="composite"）：

```json
{"prover": "composite", "kind": "varpi", "comp": ">", "q": "1", "p": "131/50",
 "rule": "pi_div_agm", "witness": {"A": "…", "B": "…"},
 "expect": [{"kind": "pi", "power": "1", "comp": ">", "bound": "A"},
            {"kind": "gauss", "power": "1", "comp": ">", "bound": "B"}],
 "children": [<pi kernel cert>, <agm cert wire-form>]}
```

`gamma_special.verify_cert` 增量：`_RULES += "pi_div_agm"`、kind 白名单 +`varpi`、新增规则块——校验 kind=="varpi"、A,B > 0、'>' 要 A·B ≥ R / '<' 要 A·B ≤ R、expect 恰为 `[(pi,1,c,A),(gauss,1,c,B)]`。每条规则各自绑定 kind（sqrt_mul→gamma14、sqrt→gamma12、sqrt_div→gamma34、pi_div_agm→varpi、pos_transfer 自指子证），白名单放行 varpi 不会让别的规则为它背书。

## 文件归属理由

放 `agm.py` 顶层而非 `kernels/`：与 `pade.py` 同为按名分派的第二证法，不走 FAMILY 矩核入口、不产 (m,n,P) 参数、无站端对应——kernels/ 留给站端可渲染的矩族。

## leader 接线清单（精确编辑）

`solve.py` `prove_exact` 的 `except NoSolution` 块，加在 pade 分支旁：

```python
if kind in ("varpi", "gauss"):
    from . import agm

    cert = agm.prove(kind, q, comp, r)
    if cert is not None:
        return {"type": kind, "prover": cert["prover"], "certificate": cert}
```

（varpi 的 cert["prover"]=="composite" 直接走既有 composite 返回通道；gauss 的是 "agm"。）

`server.py` `/calculate`，加在 pade 分支旁（varpi 的 composite 响应走既有分支、免改）：

```python
if result.get("prover") == "agm":
    from . import agm

    return respond({
        "success": True,
        "type": result.get("type", kind),
        "prover": "agm",
        "certificate": agm.cert_jsonable(result["certificate"]),
    })
```

`exact_check/__init__.py` `verify_response`，加在 pade 分支旁：

```python
if resp.get("prover") == "agm":
    from .. import agm

    ok = agm.verify_cert(agm.cert_parse(resp["certificate"]))
    return {"identity_ok": ok, "nonneg": ok, "integrand": {}, "target": {}}
```

`certificate.py` `verify_cert`，加在 composite 分支旁（数据标记分派，与 "serr" 同款）：

```python
if isinstance(cert, dict) and "agm_iter" in cert:
    from . import agm

    try:
        return agm.verify_cert(agm.cert_parse(cert))
    except (KeyError, TypeError, ValueError, ZeroDivisionError, AttributeError):
        return False
```

`certificate.py` `cert_tex`，加在 serr 分支旁（agm/composite 证书无 check 块）：

```python
if "agm_iter" in cert:
    return f"{cert['kind']}({cert['power']}) {cert['comparison']} {cert['bound']} \\quad (\\mathrm{{AGM}})"
```

可选：`decompose_exact.py` 两处 `s.proof.get("prover") == "pade"` 旁可考虑同样放行 "agm"/"composite"（本任务未测该路径）。

## 地板实测（exact 模式，power=1，前后对照）

before 用 denom-10⁴ 二分 + 十进制逐位扫描（核单独，含 exact '<' 的 m≤256 宽预算）：

| kind | comp | 核地板（最紧可证界） | 核余量 | agm 后 |
|---|---|---|---|---|
| gauss | > | 0.8253 | ~9.3e-3 | ≥1e-2000（实测），~1e-4000 顶 |
| gauss | < | 0.8351 | ~4.7e-4 | ≥1e-2000 同上 |
| varpi | > | 2.5934 | ~2.9e-2 | ~1e-16 |
| varpi | < | ~2.630（decade-2 扫描） | ~8e-3 | ~9e-17 |

gauss 的极限只由标度阶梯决定（digits=4096 实测证到 1e-2000 间隙，11 步迭代）；varpi 由 pi oracle 地板决定（quadlog pi 双向 ~1.2e-16，换算 ϖ 侧 ~ϖ/π·π地板 ≈ 0.83×pi 地板）。varpi 的 pi 子证余量 = π·(ϖ−R)/ϖ，与 R 侧间隙同阶，oracle 失败即诚实 None。

## 副产物

gamma14 的 sqrt_mul 规则经 `_child_cert("varpi", 2, c, s1)` 调 varpi——接线后 s1 可压进 AGM 覆盖的紧界（实测 2ϖ > 5.244115 成证），gamma14 的 '<' 可证窗口随之加宽。

## 偏差与注意

- `pi_cert` 对 oracle 用 `except Exception` 兜底（gamma_special._child_cert 同款）：子证崩溃只意味着"证不出"，诚实回落。
- gauss 证书内部持 Fraction（pade 惯例，cert_jsonable 上线路）；varpi 复合证书字段全字符串（gamma_special 惯例）——两种内部表示各随其兄。
- verify_cert 约定 parsed 形入参（certificate.py 分派先做 cert_parse）；直接调 wire 形会 TypeError，由分派处的 try/except 收为 False——与 pade 同契约。
- `converged` 的 64/D 饱和阈值只管紧度不管正确性：包络在任意 n 都严格成立，早停只是选标度支持的最紧对。
