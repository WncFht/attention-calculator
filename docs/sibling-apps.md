# 姊妹应用克隆侦察：/convex 与 /health

zhuyidao.net 上除注意力计算器外还有两个同域应用：`/convex`（凹凸不等式计算器）和 `/health`
（健康计算器）。本文档是克隆它们的侦察记录：端点契约、JS 调用图、数学规模评估与成本建议。
全部为实测（2026-09-15，约 30 次请求，≤1 req/s）；捕获页面见 `bench/data/site-convex.html`、
`site-health.html`、`site-health-en.html`，当日重新抓取与捕获**字节一致**。

## 0. 共性

- 两个应用的响应包络与主站一致：`json.dumps(sort_keys=True, ensure_ascii=True)` 风格——
  键按字典序排列、紧凑分隔符、非 ASCII 一律 `\uXXXX` 转义、正文以 `\n` 结尾。
  `server: waitress`，前置 `via: 1.1 Caddy`，与主站同一栈。
- 错误包络统一为 `{"error":"...","ok":false}`（注意：**带 `ok` 字段**，与 `/calculate` 的
  裸 `{"error":...}` 不同——姊妹应用是另一套（更新的）包络约定）。
- 对计算端点发 GET 均返回 **500**（方法未处理）。
- `/convex/static/bg1.png` 与本仓库 `src/attention_calculator/static/bg1.png`
  sha256 完全相同（`0634677e…`）——背景图可直接复用。
- 健康页脚里出现 `http://convex.zhuyidao.net/` 子域链接；主域 `/convex/` 路径同样可用。

## 1. /convex —— 凹凸不等式计算器

### 路由

| 路由 | 方法 | 说明 |
|---|---|---|
| `/convex/`、`/convex/en` | GET | 返回**同一份字节**；语言由 JS 按 `location.pathname.endsWith('/en')` 在客户端切换（`TEXT` 字典 + `data-i18n`），不走服务端渲染 |
| `/convex/prove` | POST | 唯一计算端点 |
| `/convex/static/*` | GET | 静态资源（bg1.png 等） |

### 请求

表单字段（`FormData` 提交；实测 **multipart 与 urlencoded 均接受**）：

- `inequality`（必填）：单变量不等式，如 `exp(x)-log(x)-261/112>0`。
- `line`（选填，页面上没有对应 UI，是隐藏参数）：用户自供分离直线，仿射式 `m*x+b`，
  用与 `inequality` 相同的原子语法解析（见下）。

语法是 **Python `ast` 解析**——错误信息直接泄漏 AST 节点 repr：
`foo(x)>0` → 400 `{"error":"unsupported atom Call(func=Name(id='foo', ctx=Load()), args=[Name(id='x', ctx=Load())], keywords=[])","ok":false}`；
`x+y>0` → 400 `unsupported atom Name(id='y', ctx=Load())`。

已证实的语法事实：

- `^` 与 `**` 等价（`x^2` 与 `x**2` 响应字节相同）——先替换再 `ast.parse`。
- 支持原子：`exp(x)`（输出记作 `e^x`/`e^{x}`）、`log(x)`（输出 `ln x`/`\ln x`）、
  `sqrt(x)`（输出 `√x`/`\sqrt{x}`）、`x^a`（**有理指数**，`x^(3/2)` 规范化为 `x^1.5`，
  latex 输出 `x^{\frac{3}{2}}`）、`x`、整数/分数/小数常数、顶层 `+ - *` 线性组合与系数
  （`2*exp(x)` 合法）。
- `exp/log/sqrt` 的参数**必须是裸 `x`**：`log(17/30)` → 400 `log only supports argument x`
  （在 `line` 字段上触发）。
- 不等号 `>`、`>=`、`<`、`<=` 均接受；`<`/`<=` 被取反归一为 `>0` 方向的 difference。
- 缺 `inequality` 字段 → 400 `{"error":"请输入一个不等式。","ok":false}`。

### 响应

200 成功包络：`{"ok":true,"result":{...}}`。`result` 全字段（键序即线上输出序）：

```json
{"curvature":{"difference":"convex","left":"convex","right":"concave"},
 "minimum":{"value":8.981904537908036e-06,"value_text":"8.98190453791e-06",
            "x":0.567143290409784,"x_text":"0.56714329041"},
 "normalized":{"difference":"e^x - ln x - 261/112",
   "difference_latex":"e^x - \\ln x - \\frac{261}{112}",
   "left":"e^x","left_latex":"e^x",
   "right":"ln x + 261/112","right_latex":"\\ln x + \\frac{261}{112}"},
 "ok":true,
 "proof":{"formula_latex":"e^x > (\\frac{30}{17}x - 1 + \\ln \\frac{17}{30}) + \\frac{261}{112} \\ge \\ln x + \\frac{261}{112}",
   "left_gap_min":8.004855962417956e-06,"left_gap_min_x":0.5679840376059393,
   "line_latex":"(\\frac{30}{17}x - 1 + \\ln \\frac{17}{30}) + \\frac{261}{112}",
   "line_text":"(30/17*x - 1 + ln(17/30)) + 261/112",
   "tangent_at":"17/30","tangent_at_latex":"\\frac{17}{30}"},
 "provided_line":null,"reason":"","status":"proved"}
```

字段语义：

- `curvature.{difference,left,right}` ∈ `convex|concave|affine`（二阶导数符号分类）。
- `normalized`：difference 拆成 `left > right`。已观察的拆分规则：带正系数的凸/仿射项留在
  左侧；负系数项取负后若为凹/仿射则移到右侧；**常数按符号归属**——正常数留在左
  （`x^2+1>log(x)` → left `x^2 + 1`），负常数移右（`-261/112` → right `ln x + 261/112`）。
- `minimum`：difference 在定义域上的数值最小值（`x_text`/`value_text` 为约 12 位有效数字
  字符串）。定义域从 `x≈1e-08` 起扫（ln/sqrt 要求 x>0），最小值点如 `0.567143290409784`
  是典型的 scipy 标量优化输出。方向归一失败（left 非凸/仿射或 right 非凹/仿射）时
  `minimum` 为 `null`。
- `proof`：仅在内置切点搜索成功时非空。`tangent_at` 为有理切点 `x0`；`line_*` 是**凹侧**
  在 `x0` 处的切线（`ln x` 在 `17/30` 处切线斜率 `30/17`，与样本一致）；`left_gap_min*` 是
  `left - line` 的数值最小值及取到位置，`formula_latex` 据此写 `left > (line) \ge right`
  （左段严格 `>` 因为数值 gap>0）。
- `provided_line`：传了 `line` 时为
  `{"b":..,"left_gap_min":..,"left_gap_min_x":..,"m":..,"ok":..,"right_gap_min":..,"right_gap_min_x":..}`，
  即把直线归约为浮点 `(m,b)` 并数值检验两侧 gap；`ok:true/false` 仅作诊断，**不回灌 proof**
  （`line=1` 对 `x^2+1>0` 验证通过但 `proof` 仍为 null）。
- `status` ∈ `proved | inconclusive | failed`，`ok` 与 `status=="proved"` 一致。

已观察的 `reason` 原文（决定渲染分支）：

| status | proof | reason 原文 |
|---|---|---|
| proved | 非空 | `""` |
| proved | null（left 严格凸，搜索未中） | `不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。` |
| proved | null（left 为仿射，跳过搜索） | `不等式数值上已通过，但当前情形没有生成中间直线证明。` |
| inconclusive | null | `整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。` |
| failed | null | `数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。` |

值得记录的怪癖（字节级复现时要注意）：

- `e^x >= x+1`（切线 `x+1` 在 x0=0 处与右侧完全重合）也报"没找到漂亮的有理切点直线"——
  候选搜索对零 gap 或边界切点不认可，说明候选集/容差有特定排除规则。
- `line_latex` 不做约简：斜率 1 印作 `1x`（`(1x - 1 + \ln 1)`），`ln 1` 不化简。
- 文本态用 `ln`、`√x`（Unicode 根号），latex 态用 `\ln`、`\sqrt{x}`、`\frac{p}{q}`。
- 页面 JS 只用 `result.ok`/`reason`/`proof.formula_latex` 渲染；
  `TEXT` 里的 `normalized/curvature/minimum/tangentAt` 键和 `.meta-grid` CSS 是**死代码**
  （旧版界面的残留），说明 result 的完整字段曾全部上屏——克隆时按全字段透传实现即可。

### JS 调用图

`proofForm.submit` → `fetch(POST /convex/prove, FormData)` → `data.ok && data.result` →
`result.ok ? normalCard : failedCard` → MathJax `typesetPromise` 渲染 `\[formula_latex\]`；
"显示TeX代码"按钮 `alert(formula_latex)`。语言切换纯前端（`history.replaceState`）。

## 2. /health —— 健康计算器

### 路由

| 路由 | 方法 | 说明 |
|---|---|---|
| `/health/` | GET | 中文页，`IS_EN=false` |
| `/health/en` | GET | 英文页，与中文页仅差静态文案和 `IS_EN=true`（diff 已确认） |
| `/health/calculate` | POST | 唯一计算端点，**仅接受 `Content-Type: application/json`** |

### 请求

JSON body，字段与表单一一对应：

```json
{"nickname":"TestUser","sex":"male","age":35,"height_cm":175,"weight_kg":70,
 "waist_cm":82,"hip_cm":96,"resting_hr":62,"pal":1.55,
 "exercise_type":"jogging","exercise_minutes":30,
 "bp_context":"clinic","systolic_bp":118,"diastolic_bp":76,
 "sleep_hours":7.5,"website":""}
```

- 必填：`nickname`（≤24 字符，中文/字母/数字/`_`/`-`）、`sex`（`male|female`）、
  `age`（18–120 整数）、`height_cm`（1–250）、`weight_kg`（1–500）、`pal`（页面只给
  1.55/1.85/2.20 三档）。
- 选填：`waist_cm`、`hip_cm`、`resting_hr`、`exercise_type`+`exercise_minutes`（**必须成对**）、
  `bp_context`+`systolic_bp`+`diastolic_bp`（三者成组）、`sleep_hours`。
- `website` 是蜜罐字段，非空即拒绝。
- 实测错误（全部 400 + `{"error","ok":false}`）：非 JSON/坏 JSON → `提交内容格式不正确。`；
  蜜罐非空 → `提交未通过校验。`；缺 sex → `请选择性别，以便使用对应公式。`；
  age=17 → `年龄应在 18～120 之间。`；昵称含非法字符 → `昵称只能包含中文、字母、数字、下划线或短横线。`；
  只有运动类型缺时长 → `运动类型和运动时长需要一起填写，或都留空。`。

### 响应

```json
{"ok":true,"record_id":220,"results":{...}}
```

- `record_id` **逐次递增**（220→221→222）——服务端把每次提交落库，克隆时自建计数器即可，
  线上绝对值无法也不必对齐。
- `results.profile`：`{age, height_cm, nickname, sex, weight_kg}` 回显。
- `results.{body, metabolism, heart, nutrition, exercise, vitals}`：六个分组的**原始数值**
  字典（如 `body` 含 `bmi/bmi_status/body_fat_pct/whr/whr_limit/rfm_pct/bai_pct/absi/bri/
  conicity/bsa_m2/janma_lbm_kg/fat_mass_kg/ffm_kg/watson_tbw_l/body_water_pct/
  blood_volume_l/target_weight_low/high/weight_adjustment_*` 等），未填写的选填项对应
  字段为 `null` 或整个分组为空。
- `results.items`：渲染用卡片数组，每项
  `{category_key, direction, display_name, display_value, formula_version, key, metadata,
  numeric_value, reference, sort_order, source_url, status, unit}`。
  全量 36 个 key：`bmi, target_weight, weight_adjustment, body_fat_cun_bae,
  body_fat_deurenberg, fat_mass, ffm, janma_lbm, bsa_mosteller, waist_assessment, whtr, whr,
  rfm, bai, absi, bri, conicity, watson_tbw, body_water_pct, blood_volume_nadler,
  mifflin_ree, harris_ree, cunningham_ree, tdee, resting_hr, max_hr_tanaka, moderate_hr,
  vigorous_hr, karvonen_hr, heart_rate_reserve, estimated_vo2max, protein_range,
  exercise_kcal, systolic_bp, diastolic_bp, blood_pressure_category, sleep_duration`。
  必填-only 提交返回 19 项；`whr`/`bai` 依赖 `hip_cm`，Karvonen/HRR/VO₂max/resting_hr
  依赖 `resting_hr`，血压三项依赖血压组，以此类推。
- `results.tips` / `tips_en`：条件触发的提示语（中英双语都下发；item 层面的 `display_*`/
  `status`/`unit`/`reference` 只发中文，英文页靠 `ITEM_NAMES_EN`/`STATUS_EN`/`UNIT_EN`/
  `REFERENCE_EN` 查表翻译——所以线上表已经把所有 `status`/`display_name` 中文串枚举在
  `site-health-en.html` 里，这是**完整的状态串清单**，克隆时直接照抄）。
- 精度按指标分别取整：bmi/bai/体重类 1 位小数，whr/whtr/absi/conicity 3 位，bri/bsa/
  血容量 2 位，心率类整数（`184.0`）。`display_value` 是字符串、`numeric_value` 是数值，
  两者同源。

### JS 调用图

`healthForm.submit` → 前端 `reportValidity` + 身高范围复查 → 组 JSON →
`fetch(POST /health/calculate)` → `showResults(data.results)`：hero 卡（bmi/body_fat_cun_bae/
tdee/max_hr_tanaka 四项）+ 按 `category_key` 分组的 `metric-grid` + `tips` 列表 +
固定免责段；英文经 `translatedItem` 查表替换。无服务端历史、无图像端点。

## 3. 数学规模评估

### /health：纯"表单 → 算术 → 模板"，无求解器

全部 ~20 个公式连同出处**逐字印在页面上**（公式与出处一节 1–11）：BMI/目标体重反算、
CUN-BAE 与 Deurenberg 体脂、腰高比/腰臀比/RFM/BAI/ABSI、Janmahasatian 瘦体重、
Mosteller 体表面积、Mifflin–St Jeor/Harris–Benedict/Cunningham 代谢、TDEE（PAL 三档）、
Tanaka 最大心率、中/高强度目标心率、Karvonen、心率储备、Uth 法 VO₂max
（`15.3×HRmax/HRrest`）、运动热量（`MET×3.5×kg/200×min`，14 个运动类型的 MET 表需从
探测或文献汇编补齐）、Watson 体水分、Nadler 血容量、BRI、锥度指数、血压/睡眠阈值分类。
需要逆向的只有：每项的舍入模式、阈值边界（≥ vs >）、tips 触发规则、MET 映射、
`record_id` 持久化——全是可探测的确定行为。**成本：约半天到一天**（含采集 golden 验证
舍入/边界），风险低。

### /convex：真求解器，符号 + 数值混合

页面自述与实测行为指向"凹凸分离切线法"：

1. `ast` 解析受限文法（原子：exp/log/sqrt/x^a/x/常数的顶层线性组合）；
2. 归一 `left > right`（按项的符号与凹凸性拆边，常数按符号归属）；
3. 对 left/right/difference 做凹凸分类（convex/concave/affine，按二阶导符号；
   分类不过 → `inconclusive`）；
4. 网格 + 标量优化求 difference 数值最小值（不达要求 → `failed`）；
5. left 为仿射则跳过证明；否则在**有限有理切点候选集**上枚举 x0，取凹侧切线，
   数值验证 `left - line > 0`，成功则组装 `formula_latex`。

克隆需要复刻：与 Python `ast` 完全一致的解析与报错串（错误信息直接 `repr` AST 节点，
等于公开了实现栈——照用 `ast` 即可对齐）、项级凹凸分类器、与线上**逐位一致**的数值
最小化（`x_text` 形态指向 scipy `minimize_scalar`/网格扫描，版本敏感）、有理切点候选
枚举（候选集与排除规则需大量探测反推，`e^x>=x+1` 这类边界行为是重点）、以及
`line_text`/`line_latex`/`formula_latex` 的不约简格式化怪癖。
**成本：约 2–5 天**，沿用主站的"采集 golden → 对齐输出"打法可行；难点在数值路径的
逐位对齐与切点候选集的反推，而非算法本身。

## 4. 建议

先做 /health（低成本、纯确定性算术，一天内可交付 byte-parity）；再做 /convex
（真求解器移植，按主站相同方法学：先建 probe/golden 管线再写求解器）。
两个应用的响应字段都超出其前端渲染所需，按"全字段透传 + 原样格式化"实现即可覆盖
将来前端升级的兼容面。
