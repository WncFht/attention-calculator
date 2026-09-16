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

失败返回 `{"success": false, "error": "..."}`（HTTP 200 或 4xx）。已观察错误：

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
