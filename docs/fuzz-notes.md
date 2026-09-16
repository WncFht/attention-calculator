# 随机差分模糊测试记录

> **已收口（2026-09-16）**：199 条分歧全部修复，`bench/replay_fuzz.py` 离线重放 922/922 全绿。本文件保留作采集快照与归因记录。

2026-09-15，`bench/fuzz.py --seed 1`：对 zhuyidao.net 与本机克隆做同请求差分。共 918 条探针（/calculate 756、/get_integral_image 162），站端限速 1 req/s，全程无 传输失败、无本机超时，**199 条分歧全部是真实行为差异**，无一 harness 工件。 探针原文在 `bench/data/fuzz-probes.jsonl`（gitignored）。

结论：站端在三个层面和我们不同——/calculate 的字段校验细节（空 comparison、 空白字符）、/get_integral_image 的逐字段校验与越界渲染、若干核路径（退化零解、 q=0 崩溃不对称、度数类型回显名）。

## /calculate 分歧

### 缺失或空 comparison 被站端容忍（41 例）

我们要求 `comparison ∈ {">","<"}`，空值即 400 `无效的不等号方向`；站端只在 comparison 非空且非法时报 400，**缺失/为空时跳过该校验继续走**（随后按 rational → power → 定义域 → 求解推进）。

- `POST {type:gamma, power:6/4, rational:116/204}`（无 comparison） → 站端 200 `{"solution":"a = 0, b = 0", ... u_val:"0", c_val:"1"}` 退化解； 我们 400 `无效的不等号方向`。
- `POST {type:pi, power:12/12, rational:1/0}`（无 comparison） → 站端 400 `右侧有理数分母不能为0`；我们 400 `无效的不等号方向`。
- `POST {type:ln_q, power:00, rational:4}` → 站端 404 `请在ln后输入一个大于1的数`。
- 校验序佐证：`{power:10, comparison:"><", rational:"∞"}` 站端报 `无效的不等号方向`——非空非法 comparison 仍在 rational 之前检查，与我们一致。

query 串带参（`POST /calculate?type=...&...=`）41 例中的 21 例同属此类：站端 同样只读 POST body（query 里 rational 合法仍报 `右侧有理数格式无效`，证明它 看到的是空表单），只是它先炸在 rational、我们先炸在 comparison。传输语义本身 无分歧。

### 数字字段空白字符规则不同（17 例）

站端把 ASCII 空格当透明字符（`"3 /4"`、`"1 / 2"`、`"1 2"`→12、`"  7  "` 全部 按去空格后解析），还容忍尾随换行（`"1\n"`→1）；但拒绝前导 tab（`"\t1"` 400）。 我们只 strip 两端、拒绝一切内部空白，而 strip 会吃掉两端 tab——双向都有分歧：

- `{type:cosh_q, power:"3 /4", comparison:<, rational:10/4}` → 站端 200 正常解 （`a = 17/8, b = 9/4, c= -9/8`）；我们 400 `左侧系数格式无效`。
- `{type:gauss, power:"3/ 4", comparison:>, rational:3/4}` → 站端 404 `方向反了`（power 已解析）；我们 400 `左侧系数格式无效`。
- `{type:cos_q, power:1571/1000, comparison:>, rational:"3/ 4"}` → 站端 404 域 错误（rational 已通过）；我们 400 `右侧有理数格式无效`。
- `{type:catalan, power:"\t1", comparison:>, rational:1707/1864}` → 站端 400 `左侧系数格式无效`；我们 200 正常解（strip 成 1）。
- `{type:tanh_q, power:6, comparison:<, rational:"\t1"}` → 站端 400 `右侧有理数格式无效`；我们 200。

### GET /calculate 不 405（10 例）

站端对 GET 返回 JSON 500 `{"error":"服务器内部错误，请稍后再试"}`（路由不区分 method 进 catch-all）；我们 Flask `@app.post` 返回 405 HTML 错误页。

### 响应 type 回显内核名而非请求名（4 例）

站端成功响应的 `type` 字段回显内核族名：`sin_q_degree`→`"sin_pi_q"`、 `cos_q_degree`→`"cos_pi_q"`；我们回显请求原值。golden.jsonl 全部 degree 成功 记录本就如此，parity.py 只比 `parameters`+`solution` 从未比 `type`，所以一直 没暴露。

### 退化零解与核路径不对称（7 例 + 2 例双 200 参数不同）

站端在我们搜不到解的位置返回 `solution:"a = 0, b = 0"` 的退化成功，参数内部 不自洽（`a_val:"0"` 却 `au_val` 非零、`u_val:"0"`、`cu_val:"1"` 或 `"2"`）。 `sqrt_bound_proof` 已复现 `cu_val:"1"` 变体，但还有 `cu_val:"2"` 变体和 gauss/varpi 的其它退化形状未覆盖：

- `{type:varpi, power:4/12, comparison:<, rational:20/8}` → 站端 200 `{a_val:"0",au_val:"7",b_val:"4/3",bu_val:"0",cu_val:"2",u_val:"1"}`； 我们 404 `未找到<方向的解`。
- `{type:varpi, power:0/1, comparison:<, rational:8258/3149}` → 站端 200 同型； 我们 404。
- `{type:gauss, power:"00", comparison:<, rational:2/9}` → 站端 200 `{au_val:"10",b_val:"5/27",cu_val:"2",u_val:"9"}`；我们 500—— `transposed_lt_proof → lhs_mpf` 里 `r / c` 对 c=0 抛 ZeroDivisionError。
- `{type:gauss, power:8473, rational:313/373}`（无 comparison）→ 站端 200 退 化（`au_val:25419=3·power`）；我们 400。
- `{type:gauss, power:7/12, comparison:<, rational:3583/4293}` → 双 200：站端 `a=0,b=0`（`au_val:17915,b_val:5779/51516`），我们正常大解 `a = 556560875/199386, b = 580641061/244224`。
- `{type:gauss, power:10^15/9999999999999999, comparison:<, rational:3174/7078}` → 双 200：站端退化零解且参数疑似 64 位溢出垃圾 （`b_val:"19371999999999997355/70779999999999992922"`），我们精确大数解。

双曲核 q=0 崩溃方向相反：

- `{type:cosh_q, power:0, comparison:>, rational:9}` → 站端 404 `方向反了`； 我们 500（我们的 cosh 核在 q=0 崩）。
- `{type:sinh_q, power:0, comparison:<, rational:01/02}` → 站端 500；我们 404。
- `{type:tanh_q, power:0, comparison:<, rational:7762}` → 站端 500；我们 404。
- `{type:cosh_q, power:10^15/1, comparison:<, rational:00}` → 站端 500；我们 404。

## /get_integral_image 分歧

### 站端逐字段 400 校验，我们笼统 500 或直接渲染（108 例，最大类）

站端校验序：type（`无效的证明类型`）→ comparison（`无效的不等号方向`）→ 整数字段（`m/n/au_val/bu_val/cu_val/u_val必须是整数`，缺失值同罪）→ 分数字段 （`a_val/b_val/c_val格式无效`、`c_val分母不能为0`）→ 范围（`m/n/u_val过小`、 `m过大`；已见 m=-3、u_val=0、n=-10^16 触发过小，m=10^19 触发过大）。

- 我们 `coerce_params` 一把 int()/Fraction()，任何字段坏 → 500 `服务器内部错误`（80 例）。
- comparison 非法时我们根本不校验、照渲染 200：如 `comparison:"≥ "` 或 swap 进 junk 值 → 站端 400 `无效的不等号方向`，我们 200（15 例）。
- `u_val:"0"` → 站端 400 `u_val过小`，我们渲染 200；`n:-10^16` → `n过小` 同理（13 例）。

### 站端渲染我们 500 的输入（4 例）

- `rational:"\\frac{0}{0}"` → 站端原样 `\frac{0}{0}` 回显进式子；我们 `wire_fraction` 除零 500。
- `rational:"1/"` → 站端渲 `\dfrac{1}{}`（空分母回显）；我们 `Fraction("1/")` 炸。
- 缺 `type` → 站端默认 pi 渲 `3\pi - \frac{-6}{0} = ...`；我们 `kind=""` → 500。
- `type:pi_n, coef:"\\frac{-9}{8}"` → 站端渲 `(\dfrac{311}{10})^{8} - (\pi^\frac{-9}{8})^{8} = ...(\ln(1/x))^8` （分母当外层指数与 ln 指数）；我们 `spec()` 拒绝分子越界 500。

### 双 200 但 equation 文本不同（6 例）

- `type:tan_q, coef:"\\dfrac{3}{1}"` → 站端函数内层剥壳写 `\cos(3/1)`、表面 `\tan\dfrac{3}{1}` 回显；我们内层仍 `\dfrac{3}{1}`。
- `type:e_q, coef:"\\dfrac{6}{-4}"` → 站端被积 `e^{\frac{3x}{2}}`（负号丢 失/取绝对值）；我们 `e^{-\frac{3x}{2}}`。
- `type:ln_q, coef:"\\frac{-4}{2}"` → 站端分母 `(2x+2)^6`；我们 `(-3x+1)^6` （负 ln 参数的分母推导不同）。
- `type:tan_q, coef:-10^29` → 站端 `\dfrac{1}{\cos(-10^{29})} - x\sin(10^{29}x)`； 我们 `... x - \sin(10^{29}x)`（负参数的因子排布不同）。
- `type:sin_q_degree, au_val:-10^23` → 站端 `(...)^{3} \left(...` 无 `\cdot`； 我们有 `\cdot`。

## 无分歧但值得记录的面

- Unicode 十进制数字（`３`、`٣`、`１２/３`、`٣/٤`）两端都按 Nd 类数字接受并 解析为真值（站端 200 正常解）；`²`、`½`、`∞`、`NaN`、`inf` 两端都拒绝。
- 缺失 `type` 字段两端都默认 pi；`type` 非法值两端都 400 `无效的证明类型`。
- query 串参数两端都忽略（站端不是 request.values 合并语义）。
- 右侧 `n/d` 上限 10^16、零分母、常见畸形（`1.5`、`1e3`、`+1`、`-1`、`0x10`、 `1_0`）两端判定一致。

## 复现

```bash
.venv/bin/python bench/fuzz.py --seed 1 --calculate 750 --image 160 \
    --max-minutes 45
```

同 seed 生成同序列探针；`fuzz-probes.jsonl` 已按 (endpoint, request) 去重， 中断后重跑自动续采。分歧重放：从 jsonl 取 `request.form` POST 给两端比对即可。
