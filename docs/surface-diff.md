# zhuyidao.net 站点表面审计

2026-09-15 对 https://zhuyidao.net/ 全表面抓取，与 `src/attention_calculator/` 复刻版逐项对比。线上栈：waitress（`server: waitress`）经 Caddy 反代（`via: 1.1 Caddy`，HTTP/2，http→https 308）。

## 站点全貌

线上是三个子应用挂同一域名，注意力计算器只是其中之一：

- `/`（注意力计算器）：我们复刻的对象。同一应用还挂在 `/attention` 前缀下（HTML 相同但静态资源路径改写为 `/attention/static/...`），`/convex`、`/health` 疑似 Caddy 按前缀路由到各自后端。
- `/convex`（凹凸不等式计算器）：独立页面 14501 B，调 `POST /convex/prove`，静态在 `/convex/static/`。子域 `convex.zhuyidao.net` 也指它（308 跳转）。错误体格式 `{"error":"...","ok":false}`，108 B，与主站不同。
- `/health`（健康计算器）：独立页面 53521 B，`/health/` 相同；`/health/en` 是**真正的英文页**（55657 B，服务端渲染，非客户端 i18n）。调 `POST /health/calculate`。`/health/en/`（带尾斜杠）500。`/health/static/*` 路径返回 500 而非 404。

我们只复刻了 `/`，`/convex`、`/health` 两个子应用整体缺失。要不要克隆它们由你定。

## 端点清单

| 端点 | 线上行为 | 本地行为 | 状态 |
|---|---|---|---|
| `GET /` | 200，99007 B（含 UTF-8 BOM） | 200，73138 B（无 BOM） | 已实现，字节级不一致（见下） |
| `GET /en` | 200，与 `/` **逐字节相同**（md5 一致）；语言靠前端 `currentLang = pathname === '/en'` 判定 | 200，同模板 | 一致（无需 en 模板） |
| `GET /en/` | 404 JSON | 200 页面 | **不同**：我们多挂了 `/en/` |
| `GET /attention`、`/attention/`、`/attention/en` | 200，99037 B（`/static` → `/attention/static`） | 404 | **缺失** |
| `POST /calculate` | 正常 | 正常 | 已实现（细微差异见下） |
| `GET /calculate` | 500 `{"error":"服务器内部错误，请稍后再试"}` | 405 HTML | **不同** |
| `POST /decompose_inequality` | 正常 | **500 HTML——`from . import decompose` ImportError，包里没有 `decompose` 模块（只在 `bench/decompose_model.py`），该端点整体坏掉** | **已挂但坏** |
| `GET /decompose_inequality` | 500 JSON | 405 HTML | **不同** |
| `GET /get_integral_image` | 有序 400 校验（见下） | 一切失败都 500 | **不同** |
| `POST /get_integral_image`、`POST /`、`POST /en` | 500 JSON catch-all | 405 HTML | **不同** |
| 任意未知路径 | 404 `{"error":"请求的页面不存在"}`（application/json，61 B） | 404 HTML（Flask 默认页，207 B） | **不同**：缺 errorhandler(404) |
| `/favicon.ico` | **204 No Content**，`content-type: text/html`，0 B | 404 HTML | **不同**（站端显式吐 204） |
| `/robots.txt`、`/sitemap.xml` | 404 JSON | 404 HTML | 一样不存在，仅响应格式不同 |
| `/static/bg1.png` `bg2.png` `title.png` | 200 | 200 | **逐字节一致**（md5 全对上） |

线上还有这些细节：路由精确匹配（`/En`、`/EN` 均 404）；`HEAD /calculate` 也 500；`OPTIONS /` 200。

## `/calculate` 行为差异

校验顺序两端一致（type→comparison→rational→power 内部再 格式→分母→上限），但**缺省处理不同**：

- 线上：字段**缺席**时 `type` 默认 `pi`、`comparison` 默认 `>`；字段**存在但为空**时按 400 报错。实测 `POST rational=3` → `400 {"error":"左侧系数格式无效"}`（comparison 缺席被默认掉了，直接落到 power 校验）。
- 本地：`request.form.get("comparison", "")` → `''` → `400 无效的不等号方向`。即同样的请求线上走更远、报错文案不同。

## `/get_integral_image` 行为差异

线上是逐字段有序校验，全部 400；本地是 `coerce_params` 一把梭，失败即 500。

线上校验顺序与文案（全部实测）：

1. `type`：缺席→默认 `pi`；`type=` 或 `type=zzz` → `400 无效的证明类型`
2. `comparison`：缺席→默认 `>`；空或非法 → `400 无效的不等号方向`
3. `m`、`n`、`u_val`、`au_val`、`bu_val`、`cu_val`：必须整数，否则 `400 {k}必须是整数`；越界 `400 {k}过小`/`{k}过大`。实测上限 `m ∈ [0,30]`、`n ∈ [0,50]`；`u/au/bu/cu_val` 只见下限（-1 → 过小），200 仍 200
4. `a_val`、`b_val`、`c_val`：分数格式（`1/2`、`-1/2` 均可），否则 `400 {k}格式无效`；分母 0 → `400 {k}分母不能为0`；分子无上限
5. `coef`、`rational`：**完全不校验**，原样回显进 LaTeX（`coef=x` → `x\pi - ...`，`rational=x` → `\pi - x`）。本地反而用 `render.wire_fraction` 先解析一遍，不可解析 → 500，比线上严

另外本地在 `type` 缺席时 `kind=''` → 500，线上默认 `pi` 正常出式子。

合成参数 `type=pi&m=1&n=1&a=b=c=u=au=bu=cu=1&coef=1&rational=3` 两端输出文本不同：线上 `\pi - 3 = \int_0^1 x^{2} \cdot \left(1 - x^{2}\right) \mathrm{d} x > 0`，本地 `\pi - 3 = \int_0^1 - x^{4} + x^{2} \mathrm{d} x > 0`（数学等价、因式分解形态不同；真实 /calculate 回参走 golden.jsonl 校准过的路径，此属边角输入的渲染分支差异）。

## `/decompose_inequality` 行为差异

- 线上空/缺 `problem` → `400 {"error":"请输入一个只包含一个 > 或 < 的不等式"}`；`problem=foo`、`x>y>z` 同此文案。本地文案是 `请输入一个组合不等式。`，且非空输入直接 ImportError → 500 HTML 页。**端点当前完全不可用，且 500 返回的是 HTML 不是 JSON**。
- 线上正常返回字段：`basic_count, decomposition_latex, normalized_latex, problem, steps[]`，steps 含 `bound, bound_latex, coefficient, comparison, equation, label, type`。

## HTML 页面 diff（live 99007 B vs 模板 73139 B）

`/en` 与 `/` 逐字节相同，语言切换全在前端（`data-i18n` + `UI_TEXT` 字典 + `tr()`），无需单独 en 模板。本地模板是把线上 HTML"清理"过的版本，差异如下（行号指线上文件）：

1. **BOM**：线上文件以 UTF-8 BOM（EF BB BF）开头，本地没有。
2. **CSS**：
   - 线上第 ~299 行 `width: 5rem;    统一宽度` ——注释符丢了，`统一宽度` 以裸文本留在 CSS 里（非法声明，浏览器丢弃）；本地是正常注释 `/* 统一宽度 */`。
   - 线上多出 `.pi-input-group`、`.pi-container`（重复定义）、`.pi-fraction-input`、`.fraction-input input, .pi-fraction-input input` 四段规则；本地缺。
   - 多处行尾空格差异（如 `margin: 0 auto; `）。
3. **HTML 注释**：线上保留大量被注释掉的遗留选项/符号（`arcsin_q`、`arccos_q`、`arctan_q`、`arccot_q`、`sin_pi_q`、`cos_pi_q`、`Gamma_1_3_2_3`、`psi`、`zeta`、`erf`），本地全部删除。
4. **选项文案空白**：线上 `\(\zeta(3)\)`、`\(\varpi\)`、`\(\sin q^\circ\)(角度)`，本地多了空格 `\( \zeta(3) \)`、`\( \varpi\)`、`\(\sin q^\circ\) (角度)`。
5. **结果区结构**：线上初始 DOM 有**两个重复 `id="result"`** 的 div——第一个含 `#showTexBtn`，第二个含 `#exportBtn`（`showResult()` 重写 innerHTML 后才归一）。本地合并为单个 `.result-header`/`.result-actions`。
6. **JS 逻辑**：
   - 系数校验：线上是 ~28 个 `else if (selectedType === '...')` 平铺链；本地重构为 `POWER_PROMPTS` 查表。提示文案逐条核对一致（含线上 bug：`pi_n` 分数模式提示误写"请输入e的系数的分子和分母"，本地如实复刻）。等价但非逐字节。
   - 线上多出一个内联 `Fraction` 分数库（`!function(a,b){...}`，挂 `window.Fraction`，页面内未见调用）、`window.load` 里的 `MathJax.typesetPromise` 处理（引用已被注释掉的 `#proofType` → `select` 为 null → `.catch` 吃掉并在 console 打 `MathJax error:`）、`#eProof` 的 toggle 绑定（`#eProof` 不在 DOM 中，死代码）。本地均删。
   - **历史记录**：线上 `sessionHistory` 只在内存，**不写 sessionStorage**；本地加了 `HISTORY_KEY`/`loadHistory()`/`sessionStorage.setItem` 持久化——行为差异：线上刷新页面历史即清空，本地会恢复。
   - `showResult` 内部线上多取 `showTexBtn`/`exportBtn` 局部变量等零碎差异。
7. **fetch 目标、字段名、错误显示路径完全一致**：`POST /calculate`（form: type/power/comparison/rational）、`POST /decompose_inequality`（form: problem）、`GET /get_integral_image?m&n&a_val&b_val&c_val&u_val&au_val&bu_val&cu_val&type&comparison&coef&rational`。MathJax 3 + KaTeX 0.16.9 + html2canvas 引用一致。

## 结论清单

缺失端点：`/attention`（含 /attention/、/attention/en）、`/convex` 子应用全套、`/health` 子应用全套、`/favicon.ico` 204。

行为不一致：`/en/`（应 404）、JSON 404/500 错误体（vs Flask HTML 405/404/500）、`/calculate` 缺席字段默认值、`/get_integral_image` 逐字段 400 校验 + type/comparison 默认值 + coef/rational 免校验回显、`/decompose_inequality` 错误文案且本地整体 ImportError 坏掉、历史记录 sessionStorage 持久化（本地比线上多）。

字节差异：BOM、CSS 残注释与缺失规则、HTML 注释、选项空白、双 `#result`、if-else 链 vs 表驱动、内联 Fraction 库、load/eProof 死代码、若干行尾空格。

资源：bg1/bg2/title.png 与线上逐字节一致，无需补下载；无可下载的 favicon（线上 204 空响应）。
