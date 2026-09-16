# zhuyidao.net API 实测协议

## POST /calculate

`Content-Type: application/x-www-form-urlencoded`，字段 `type, power, comparison, rational`。

- `type`：29 种之一：`pi, e, pi_n, e_q, ln_q, ln_q_square, sin_q, cos_q, tan_q, cot_q, sin_q_degree, cos_q_degree, sin_pi_q, cos_pi_q, arctan_q, arccot_q, sinh_q, cosh_q, tanh_q, coth_q, artanh_q, arcoth_q, gamma, golden, catalan, zeta3, e_pi, varpi, gauss`
- `power`：系数/参数，整数或 `n/d` 分数（如 `1`, `3/2`）；pi/e/常数类传 `1`（系数场景下是"几倍"）
- `comparison`：`>` 或 `<`
- `rational`：右侧有理数 `n/d` 或整数

成功返回（实测）：

```json
{"success": true, "type": "pi",
 "equations": {"solution": "a = 47/120, b = -13/120"},
 "parameters": {"m": 3, "n": 3,
   "a_val": "47/120", "b_val": "-13/120", "c_val": "0",
   "au_val": "47", "bu_val": "-13", "cu_val": "0", "u_val": "120",
   "unified_form": {}}}
```

- `a_val,b_val,c_val` = 待定多项式系数（分数）；`au_val,bu_val,cu_val` = 通分后整分子；`u_val` = 公分母。
- `m,n` 为基函数指数（各核族定义见 kernel-spec）。

失败返回裸 `{"error": "..."}`（无 `success` 键；HTTP 按错误类别分：格式校验 400、域检查与搜索失败 404、内核异常一律 catch-all 500——详见 behavior-notes §1/§10.3）。已观察错误：

- `要证明的式子不等号方向反了`（方向反，服务端可判定真方向）
- `在指数不超过{cap}的范围内未找到<方向的解`（搜索预算内无解；cap 按 type 插值：e/pi=30，其余=10）

## GET /get_integral_image

参数：`/calculate` 返回的 `parameters` 全部字段 + `type, coef(=power), comparison, rational`。

返回 `{"equation": "<LaTeX，无定界符>"}`，如：

```
\dfrac{22}{7} - \pi = \int_0^1 \frac{x^{6} \left(1 - x^{2}\right)^{3} \cdot \left(47 - 13 x^{2}\right)}{120 x^{2} + 120} \mathrm{d} x > 0
```

## POST /decompose_inequality

字段 `problem`（如 `pi^2+8*pi>35`、`e*pi+phi+sin(1)<11`）。

```json
{"success": true, "problem": "...", "normalized_latex": "...",
 "direct_basic": false, "basic_count": 2,
 "decomposition_latex": "\\pi^{2}+8\\pi>\\dfrac{227}{23}+\\dfrac{578}{23}=35",
 "steps": [{"type": "pi_n", "label": "\\pi^{2}", "coefficient": "1",
            "comparison": ">", "bound": "227/23", "bound_latex": "\\dfrac{227}{23}",
            "equation": "<完整积分式>"}]}
```

`direct_basic: true` 时输入本身即基本类型，无 decomposition。

## 其他端点

- `/convex`（凸函数不等式计算器）、`/health`（健康计算器）为姊妹应用，已复现——协议细节见 `docs/sibling-apps.md`。
- 前端：服务端渲染 + MathJax + KaTeX + html2canvas；中英双语；sessionStorage 历史 10 条。

## 本地扩展（非站端行为）

以上全是站端实测协议。本实现另加以下扩展，站端不存在：

- `POST /calculate` 表单加 `mode=exact`：走数学正确性路径（HTTP 入口 `solve.prove(exact=True)`；`prove_exact` 是其下的 `(module,kind,q,comp,r)` 内部函数）。语义差异：方向判定经 `certified_cmp` 递增精度认证（不做 float64 预检/兜底）；接受 `EXACT_TYPES` 里的 exact-only 类型（site 模式按站端口径 400）；响应附 `certificate` 字段（机器可检证书，格式 `docs/2026-09-16-certificate-spec.md`，离线复核 `tools/verify_cert.py`）；假恒等式/失真簇不发——发不出就 `未找到解`。缺省或 `mode` 为其他值时行为与站端逐字节一致。
- exact 模式下 (m,n) 搜索耗尽有四种第二证法兜底，响应统一为 `{"success","type","prover","certificate"}` 包络——无 `parameters`/`equations`，证书即证明：`prover:"pade"`（ln_q/arctan_q 回落 Padé 插值，schema 见证书规格 §Padé）、`"agm"`（gauss 的 AGM 区间包络）、`"euler_gamma"`（gamma 的 Euler–Maclaurin 包络，子证为 ln2 的 pade 证）、`"composite"`（varpi 的 pi+AGM 复合，以及 gamma14/34/12 的 Γ-型子证明 DAG）。
- `POST /decompose_inequality` 同样认 `mode=exact`：走 `decompose_exact` 的可证界分配，响应形状与站端不同（`steps`/`term_bounds`/`all_proved`/`certifies`/`failures`/`slack` 等字段）；分解器非 ValueError 异常报 500 `组合证明生成失败，请稍后再试`（站端无此错误形状）。
- `GET /demo`：本地演示页——与 `/` 同页，但把 MathJax/KaTeX/html2canvas 三个 CDN URL 改写为 `/static/vendor/` 本地副本（vendor 资产由 `scripts/fetch-vendor-assets.sh` 拉取）；`/`、`/en` 保持站端逐字节不动。
- exact-only 类型的 equation 渲染不经 `/get_integral_image`（该端点仍只认站端 29 型），由核的 `render_equation` 程序化产出；gamma14/34/12 的 composite 证书无 parameters，不产渲染式。
