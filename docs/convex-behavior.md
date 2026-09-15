# /convex/prove 行为规范

zhuyidao.net `/convex/prove` 的完整行为规格，供字节级复现。全部结论来自实测：2026-09-15
共 333 次请求（平均约 0.125 req/s）；原始语料 `bench/data/convex-probes.jsonl`（333 行，每行含
完整请求、verbatim raw 响应、tag），求解路径子集 `bench/data/convex-golden.jsonl`
（323 行 = 全部 POST `/convex/prove` 记录：proved 135 / failed 65 / inconclusive 25 /
400 错误 98，覆盖全部分支与全部错误串模板；其余 10 条为传输层路由/方法探针）。本文档
取代 `docs/sibling-apps.md` §1 的细节（其中两处实测修正见 §12）。

站点行为与作者开源参考实现 `lianghuatiaojiushi/ConvexConcaveProver` 的 `scripts/prove.py`
逐位一致（文章例题本地复算 `left_gap_min=8.004855962417956e-06`、`x_min=0.567143290409784`、
切线式 `(30/17*x - 1 + ln(17/30)) + 261/112` 与站点输出完全相同）[^prover-repo]。所有常数
（EPS、容差、连分数参数、二分上界、渲染分母上限）均已实测验证。站点 = 参考实现 +
JSON 包络 + latex 渲染层。

语料统计：333 条记录（200: 227，400: 99，500: 7）。

## 1. 传输层

- 唯一计算端点：`POST /convex/prove`。`multipart/form-data` 与
  `application/x-www-form-urlencoded` 均接受；JSON body、text/plain body **不解析**
  （字段读不到 → 400 `请输入一个不等式。`）。query string 不参与 POST 参数。
- 表单字段：`inequality`（必填）、`line`（隐藏参数，页面无 UI）、`domain`（隐藏参数）。
  其余字段一律忽略（`foo`、`lang`、`no_line`、`no-line`、`noline` 均实测无效果——
  参考实现有 `--no-line` 开关，站点未接线）。
- 路由与方法：对本路由发 GET/HEAD/PUT、POST 到 `/convex/prove/`（尾斜杠）、
  `/convex/Prove`（大小写敏感）、`/convex/` → 均 500
  `{"error":"服务器内部错误，请稍后再试。","ok":false}`（带句号、带 `ok` 字段；
  HEAD 响应体为空属正常）。OPTIONS → 200 空响应。`GET /convex/en` → 200 HTML 页面。
- 包络：`json.dumps(sort_keys=True, ensure_ascii=True)` 风格——键字典序、紧凑分隔符、
  非 ASCII `\uXXXX` 转义、正文以 `\n` 结尾。成功 `{"ok":true,"result":{...}}`；
  失败 `{"error":"...","ok":false}`。`result.ok` 恒等于 `status=="proved"`。
- 前置校验（只作用于 `inequality`，按 strip 后文本）：空 → 400 `请输入一个不等式。`；
  `len > 500` → 400 `输入过长，请输入一个较短的单变量不等式`（恰好 500 通过；
  `"x"*500 + " "` 即 raw 501/strip 500 → 通过，证明长度按 strip 后计）。`line`、`domain`
  无长度限制（`line` 送 3001 字符进入了解析）。
- 求解器内部抛出的异常一律 400 `str(exc)`，包括非 ValueError：`x/0>0` →
  `float division by zero`（ZeroDivisionError）、`x^1e309>0` →
  `cannot convert Infinity to integer ratio`（OverflowError）、`line` 超长递归 →
  `maximum recursion depth exceeded`（RecursionError）。**未观察到任何由合法不等式
  触发的 500**；500 只出现在传输层（上一条）。

## 2. 输入语法

处理流水线：`strip` → 字符串改写 → `ast.parse(mode="eval")` → 项收集。

改写（按序）：`^` → `**`；`\bln\s*\(` → `log(`；`\be\s*\*\*\s*x\b` → `exp(x)`
（只匹配裸 `x`，不带括号）。推论：

- `ln(x)`、`e^x`、`e**x` 合法；`e^(x)` **不合法**（正则不匹配，剩下 `e**(x)` →
  `unsupported atom BinOp(left=Name(id='e', ...)`）；`ln x`、`exp x`（缺括号）→
  `invalid syntax`；`E^x`、`EXP(x)`、`Log(x)` 大小写不符 → `unsupported atom`。
- `x ^ 2`、`  x^2  >  0  `（任意空白）合法——`ast` 忽略空白。

不等号在改写后的文本上搜索，顺序 `>=`、`<=`、`>`、`<`、`=`，首个命中即按子串切一刀
（`split(op, 1)`）。`<`/`<=` 交换左右后按 `>`/`>=` 处理（`log(x)<exp(x)` 归一为
`e^x > ln x`）。**无不等号** → 右端 `0`、按 `>=` 处理（裸 `x` → `x>=0`，proved）。
`=` 按非严格处理（`x^2-2*x+1=0` 输出 proof，`formula_latex` 中显示 `\ge`）。
`==`、`!=`、`>>` 等 → `invalid syntax`；`x>0>1` 第二段残留的 `Compare` →
`unsupported atom Compare(...)`。

原子（`parse_atom`，`x` 是唯一变量）：

| 形态 | 结果 |
|---|---|
| `x` | linear |
| `exp(x)` / `log(x)` / `sqrt(x)` | exp / log / power(0.5) |
| `x^num` | power(num)，`num` 为数值常量表达式 |

`Call` 检查顺序：先查参数个数（≠1 → `{name} expects one argument`，任意函数名都走到这，
如 `max(x,1)` → `max expects one argument`），再查参数是否裸 `x`（否 →
`{name} only supports argument x`，如 `log(2)`、`exp(2*x)`、`log(log(x))`），最后按名分派；
名不在 `{exp,log,sqrt}` → `unsupported atom Call(func=Name(id='foo', ...)`。
`x^x`、`x/(x+1)`、`1/x` 等 → `expected numeric constant, got {ast.dump}`。
其余一切 → `unsupported atom {ast.dump(node)}`：`y`、`x,x`、`[x]`、`{x:1}`、`nan`、
`2^x`、`(-x)^2`、`(x+1)^2`、`exp(x)^2`、`x*(x+1)` 等。

数值常量（`parse_number`）：int/float 字面量（含 `0x10`、`1_000`、`1e3`）、一元负号、
`a/b`、`a**b`（即 `a^b`，**只用于常量**）、`sqrt(const)`（`sqrt(4)` → 2；
`x^sqrt(2)` → 指数 1.41421356…）。`1e309` → inf → 渲染时
`cannot convert Infinity to integer ratio`；`x^(1/0)`、`x/0` → `float division by zero`；
`nan`/`pi` 等名字不是常量也不是 `x` → `unsupported atom Name(...)`。

项收集（`collect_terms`）：`a+b`/`a-b` 按符号递归；`a*b` 任一侧是数值常量则把系数乘进去
继续递归（`2*(x+1)` → `2x+2`，`(x+1)*3` → `3x+3`——**常数系数对和式可分配**）；两侧都
非常量 → `parse_product`：常数×原子折叠，原子×原子 →
`products of two non-constant atoms are unsupported`。`a/b` 分母为数值常量则系数除以它
（同样可分配：`(x+1)/2` → `0.5x+0.5`），否则报 `expected numeric constant`。一元 `-`
翻转符号（`--0.5` → `+0.5`，双负号可用）。`x^-1/2` 按 `x**(-1)/2` 解析 → `1/2*x^-1`。

## 3. 归一化输出 `normalized`

差式项 = 左侧项 + 右侧项取负，`combine_terms` 按 `(原子种类, 指数)` 合并：系数相加、
丢弃 `|c| ≤ 1e-9`（**含等号**：`1e-9*exp` 整个消失，`1.1e-9*exp` 保留）、
**保留首次出现的顺序**（dict 插入序：`x^2+x` 与 `x+x^2` 的 left 文本不同——
前者 `x^2 + x`，后者 `x + x^2`）。然后 `split_positive_negative` **纯按系数符号**拆分：
`coeff ≥ 0` 留左，`coeff < 0` 取负移右。与曲率无关（修正 `sibling-apps.md` §1 的
"若为凹/仿射则移右"说法：`-x^2+10>0` → left `10` right `x^2`；正系数凹项照样留左——
`x^2+sqrt(x)>0` → left `x^2 + √x`，left 判为 mixed）。差式全正时 right 为 `0`。

字段：`{left, left_latex, right, right_latex, difference, difference_latex}`。
渲染规则（站点自带 fmt，非常数一律经 `Fraction(float).limit_denominator(1000000)`）：

- 文本态：常数 `p/q` 或整数；系数 `p/q*原子`（`0.123456789` → `10/81*x`，
  `1e-6` → `1/1000000*x`，`1.1e-9` → 显示 `0` 但项仍在内部生效）；系数 1 省略、
  -1 显示 `- `。原子：`x`、`e^x`、`ln x`、`√x`（`a=0.5` 一律显示根号，含
  `x^0.5`、`x^(1/2)`、`x^(2/4)`、`sqrt(x)`）、`x^{a:g}`（`g` 格式 6 位有效数字：
  `x^0.666667`、`x^1.41421`、`x^-0.5`、`x^1000`；`x^0` **不化简**，原样输出）。
- latex 态：分数 `\frac{p}{q}`、系数并置（`2x`、`\frac{1}{10}x`）；原子 `x`、`e^x`、
  `\ln x`、`\sqrt{x}`、幂 `x^{...}`：指数 `limit_denominator(1e6)` 分母为 1 → `x^{-1}`、
  `x^{2}`、`x^{1000}`，否则 `x^{\frac{p}{q}}`（负号进分子：`x^{\frac{-1}{2}}`；
  `sqrt(2)` → `x^{\frac{665857}{470832}}`）。
- 注意不对称：**文本指数用 `{a:g}`（`x^0.666667`），latex 指数用 Fraction
  （`x^{\frac{2}{3}}`）**；系数两侧都用 Fraction。

## 4. 曲率分类 `curvature`

每项符号 = `sign(coeff × 原子二阶导符号)`，`|值| ≤ 1e-9` 判 0（EPS 作用在乘积上：
`x^1.0000000001` → `a(a−1)=1e-10` 判 affine；`x^1.000000001` → `1.000000001e-9` 刚越界
判 convex）。原子二阶导符号：const/linear → 0，exp → +1，log → −1，power → `a(a−1)`
（对 `x>0` 恒定号）。集合分类：`affine`（全 0）、`convex`（只含 {0,+1}）、`concave`
（只含 {0,−1}）、`mixed`（两号并存——**第四个取值，已实测**，如 left `x^2 + √x`）。

输出 `{left, right, difference}`，分别对拆分后左项、右项、合并差式计算。
模板规则：`left ∈ {convex, affine}` 且 `right ∈ {concave, affine}`，否则
`status=inconclusive`（`x^2+1<2x` 翻转后 left=`2x` affine、right=`x^2+1` convex →
inconclusive）。参考实现还有第二道 `difference ∈ {convex,affine}` 检查，但模板通过时
差式必为 convex/affine——**该分支不可达**，全语料无一命中。

## 5. 最小值 `minimum`

模板通过时计算，否则 `null`。对**差式**在 domain 上求最小——**导数变号二分，不是
scipy**（修正 `sibling-apps.md` §1）：

- `lo_eff = max(lo, 1e-8)`。有 `hi`：`d(lo) ≥ 0` → 边界 `(lo, f(lo))`；`d(hi) ≤ 0` →
  `(hi, f(hi))`；否则在 `[lo, hi]` 上二分。无 `hi`：`d(lo) ≥ 0` → 边界；否则 `b` 从
  `1.0` 起倍增直到 `d(b) > 0` 或 `b ≥ 1e8`；若 `b ≥ 1e8` 时仍 `d(b) ≤ 0` → 返回
  `(b, f(b))`——**`b = 2^27 = 134217728`**（`x^-1>0` 的 `x_text` 实测 `"134217728"`）。
  否则 160 次二分折半后取中点。
- `x_text`/`value_text` = `%.12g`；JSON 数值为原生 float repr（`1e-08`、`134217728.0`、
  `0.75`）。`x=0` 这种单调式最小值落在 `1e-08`（下界钳制）。
- 反向区间合法：`domain="10,0"` 时 `d(10) ≥ 0` → 直接报边界 `x=10, value=90`。

## 6. 判定与容差

`strict` = 原不等号为 `>`/`<`（`<` 翻转后仍 strict）。`ok = f_min > 1e-8`（strict）或
`f_min ≥ -1e-8`（非严格，**含等号**）。容差阶梯实测：

| 输入 | f_min | 结果 |
|---|---|---|
| `x^2-2x+1>0` | 0 | failed |
| `x^2-2x+1>=0` | 0 | proved |
| `…+1.000000005>0` | 5e-9 | failed |
| `…+1.00000001>0` | 1e-8 | failed（严格需 `>1e-8`） |
| `…+1.00000002>0` | 2e-8 | proved |
| `…+0.999999999>=0` | −1e-9 | **proved（假命题宽容通过）** |
| `…+0.99999999>=0` | −1e-8 | proved（边界含等号） |
| `…+0.999999989>=0` | −1.1e-8 | failed |

`status`：`ok` → `proved`（无论有无 proof）；否则 `failed`；模板不符 → `inconclusive`。

## 7. 切线候选搜索（核心模型）

触发条件：`ok` 且 `left_curv == "convex"`（**严格 convex，affine 不搜**）且
`right_curv ∈ {concave, affine}`（模板已保证）。affine-left 情形即使数值成立也不搜
（`x+1-ln x > 0`：left `x + 1` affine、right `ln x` concave、差式 convex → proved 无
proof，reason=「没有生成中间直线证明」）。

候选集 = `minimum.x` 的**连分数渐近分数**，`max_terms=12`、`max_denominator=10000`：

```
y = x_min; 至多 12 项: a = floor(y); 若 |y-a| < 1e-14 停, 否则 y = 1/(y-a)
convergent_i = 前 i 项折叠值; 保留 denom ≤ 10000 且未出现过的, 按生成顺序
```

x* < 1 时序列以 `0/1` 开头（`0;1,1,3,4,…` → `0/1, 1/1, 1/2, 4/7, 17/30, 38/67, …`）。
`x*=1e-8` → 只剩 `[0/1]`（`1/1e8` 分母越限被滤）。`x*=134217728` → `[134217728/1]`。

逐个候选 `x0`：**`x0 ≤ domain_lo` 或 `x0 ≥ domain_hi` 跳过（用未钳制的原始 domain
界）**；取凹侧（right）在 `x0` 的切线 `h`；接受条件为 `min_domain(left − h)` 通过
`> 1e-8`（strict）或 `≥ -1e-7`（非严格——**注意是 −1e-7**，比整体判定的 −1e-8 松
10 倍）。**首个成功者胜出**，搜完无果 → `proof: null`。

由此解释全部观测：

- `exp(x) >= x+1` 无证书：x* 落在下界 `1e-8` → 候选仅 `0/1` → `0 ≤ domain_lo`
  （默认 1e-8；显式 `0,inf` 则 `0 ≤ 0`）→ 候选集空。这就是「零 gap 切点不认」的真相。
- 分母上界实测含等号：`x^2+10000x−ln x>0` → `tangent_at=1/10000`（denom=1e4 通过）；
  `+10002x` → 1/10002 被滤 → 无证书。
- `domain="-1,inf"`：原始 `lo=-1` → 候选 `0` 不被排除 → affine right 时
  `tangent_at="0"` 出证书；含 `√x` 的 right 在 `x0=0` 求切线时 `0^-0.5` →
  400 `0.0 cannot be raised to a negative power`。
- 首个命中序的悬崖证据：`exp(x)-ln(x)-c>0` 随 c 增大 `tangent_at` 走
  `1 → 1/2 → 4/7 → 17/30 → 38/67`（沿渐近分数序推进，非最优 gap）。
- `>=` 与 `>` 选出不同切点：`exp(x)-log(x)-2>=0` 在 x0=1 的 gap −1.1e-16 经 −1e-7
  容差被接受 → `tangent_at=1`；同式 `>0` 版本 x0=1 被拒（需 >1e-8）→ 落到 `1/2`。

`proof` 字段：`{tangent_at, tangent_at_latex, line_text, line_latex, left_gap_min,
left_gap_min_x, formula_latex}`。`tangent_at` 是分数标签（`1`、`17/30`、`0`、
`134217728`），latex 为 `\frac{p}{q}` 或整数。只报告 `left − h` 的数值 gap——
`h ≥ right` 由凹性保证，不输出右 gap。

`line_text` = right 各项在切点的切线表达式以 ` + ` 连接（`+ -` → `- `），系数经
`limit_denominator(1e6)`，1 省略、−1 为 `-`：

| right 原子 | line_text 模板（切点 r） | line_latex |
|---|---|---|
| 常数 c | `p/q` | 同 |
| c·x | `{c}*x`（c=1 → `x`） | `{c}x` |
| c·ln x | `{c}*({1/r}*x - 1 + ln({r}))` | `{c}({1/r}x - 1 + \ln {r})` |
| c·√x（a=0.5） | `{c}*((x + {r})/(2*sqrt({r})))` | `{c}\frac{x + {r}}{2\sqrt{r}}` |
| c·x^a（a≠0.5） | `{c}*(({r})^{a:g} + {a:g}*({r})^({a:g}-1)*(x - {r}))` | `(({r})^{a_{frac}} + {a_{frac}}({r})^{a_{frac}-1}(x - {r}))` |

实测样例：`(30/17*x - 1 + ln(17/30)) + 261/112`；`(1x - 1 + \ln 1)`（latex 里 `1x`、
`\ln 1` 均不化简）；`\frac{x + \frac{1}{2}}{2\sqrt{\frac{1}{2}}}`；幂模板 latex 把
`a−1` 拆成 `^{\frac{2}{3}-1}`。exp 模板（`{c}*e^({r})(x {1-r带号})`）存在于参考实现但
**站点不可达**（exp 是凸原子，不可能出现在凹侧）。right 为空（仿射 `0`）→ 切线
`(0,0)` → `line_text:""`、`line_latex:""`。

`formula_latex` = `{left_latex} {rel} {line_latex} \ge {right_latex}`：首段 `rel` 为
`>`（strict）或 `\ge`（一切非严格输入，含 `=`/`<=`——**不反映原始符号**），中段恒
`\ge`。空 line 时中间留双空格：`x^{\frac{-1}{2}} >  \ge 0`。

## 8. `line` 参数（`provided_line`）

仅在 `status=proved` 且 `line` 非空时解析并验证；`line=""` 视为未传；主不等式出错时
`line` 不会被读（`foo(x)>0` + `line=x^2` → 报主式错误）；failed/inconclusive 时同样不解析
（`log(x)-x+10>0` + `line=foo(` → 仍 inconclusive 而非 400）。

语法与主式相同（同一套 `ln`→`log`、`e**x`→`exp(x)` 改写与原子表），但只接受
const/linear 项（`x+x`、`2*(x+1)`、`41/14*(x-12/161)` 均可）——非常量非线性原子 →
400 `--line must be an affine expression in x`（**泄漏 CLI flag 名**）。验证 =
`min(left−h)` 与 `min(h−right)` 各自在 domain 上数值最小化，`ok = 两者均 ≥ −1e-8`
（比内置证书的 −1e-7 严）。输出 `{m, b, left_gap_min, left_gap_min_x, right_gap_min,
right_gap_min_x, ok}`，**纯诊断，不回灌 proof/status**（`line=1` 对 `x^2+1>0`
`ok:true` 但 `proof` 仍为 null）。

## 9. `domain` 参数

`"lo,hi"` 形式：逗号切分须恰两段；`lo=float(...)`；`hi` 在 `{inf, infinity, +inf}`
（小写后比较）时为 `None` 否则 `float`。段数错 → 400
`domain must look like '0,inf' or '0,10'`；`float()` 失败 → 400
`could not convert string to float: '...'`（裸 ValueError 透出）。空字段 → 默认
`(1e-8, None)`。**原始 `lo` 用于候选排除（`x0 ≤ lo_raw`），而最小化用
`max(lo, 1e-8)`**——两套界并存是 `domain="-1,inf"` 怪异行为的根源（见 §7）。

## 10. `status` / `reason` 全表

五个 `reason` 原文（225 例携带 `result` 的 200 响应全集；另 2 条 200 为 OPTIONS 空响应与 `GET /convex/en` 页面）：

| status | proof | 触发 | reason 原文 |
|---|---|---|---|
| proved | 非空 | 切线证书命中 | `""`（空串，键存在） |
| proved | null | left 非 convex（affine），搜索未启动 | `不等式数值上已通过，但当前情形没有生成中间直线证明。` |
| proved | null | 搜索跑了但无候选被接受（含候选集为空） | `不等式数值上已通过，但内置有限候选搜索没有找到漂亮的有理切点直线。` |
| inconclusive | null | 模板不符 | `整理后左侧不是凸函数/仿射函数，或右侧不是凹函数/仿射函数，因此当前证明器无法处理。` |
| failed | null | `f_min` 未过容差 | `数值最小值未达到证明要求；该不等式可能不成立，或超出当前搜索范围。` |

## 11. 错误串目录（全部实测，→ 均为 400 除非注明）

前置：`请输入一个不等式。`（缺字段/空/全空白；JSON 或纯文本 body 同样走到这）、
`输入过长，请输入一个较短的单变量不等式`（strip 后 >500 字符）。

求解器透出（`str(exc)` 原样）：`invalid syntax (<unknown>, line 1|0)`、`invalid decimal
literal`（`2x`）、`invalid character '≥'|'−'|'×'`、`'(' was never closed`、
`unsupported atom {ast.dump}`（Name/Call/BinOp/UnaryOp 基、Compare/Tuple/List/Dict）、
`{name} expects one argument`、`{name} only supports argument x`、
`expected numeric constant, got {ast.dump}`、`products of two non-constant atoms are
unsupported`、`float division by zero`、`cannot convert Infinity to integer ratio`、
`0.0 cannot be raised to a negative power`（负 lo domain + 凹侧含 √x 时 x0=0 求切线）、
`maximum recursion depth exceeded`（`line` 送 3001 字符的 `x+x+…`）。

参数：`--line must be an affine expression in x`（`line=x^2`/`log(x)`/`exp(x)`）、
`domain must look like '0,inf' or '0,10'`、`could not convert string to float: '...'`。

传输层：`服务器内部错误，请稍后再试。`（500，仅路由/方法错误）。

## 12. 怪癖速查（复现时必须一致）

- `x^2+1=0` 判 proved（`=` 按 `>=`）；裸 `x` 按 `x>=0` 判 proved。
- `x^2-2x+0.999999999>=0` 被判 proved（−1e-9 在 −1e-8 容差内）——站点会"证明"假命题。
- `x^0` 不化简为 1；`x^(2/4)` 折叠成 `√x`；`1e-9·exp` 消失、`1.1e-9·exp` 显示 `0`。
- `e^x>=x+1`、`exp(x)-x>0`、`x^2+1>0` 等边界最小值情形一律"无漂亮切点"（候选集空）。
- `line_latex` 不约分（`1x`、`\ln 1`）；formula 中段恒 `\ge`，首段只反映严格性；
  空 line 留下 `>  \ge` 双空格。
- `proved` 可以没有 proof；`provided_line.ok` 不影响任何判定。
- curvature 有第四值 `mixed`；规范化拆分纯按符号不按凹凸。
- 最小值上界是 2^27=134217728（b 倍增而非 1e8 截断），x_text 打印 `134217728`。
- 内部异常全部 400 `str(exc)`，错误文本直接泄漏 AST dump / Python 内置异常信息 /
  CLI flag 名。

## 13. 参考实现对应关系

站点行为 = `lianghuatiaojiushi/ConvexConcaveProver` `scripts/prove.py` 的计算核心 +
JSON/latex 包装层[^prover-repo]。已逐位验证的对应：解析改写、`ast` 语法、EPS=1e-9
合并与曲率零带、`split_positive_negative` 按号拆分、导数二分最小化（b 倍增至 ≥1e8、
160 折半）、CF 收敛子（12 项、分母 ≤1e4）、`find_line` 容差（strict `>1e-8` /
非 strict `≥-1e-7`）、`verify_line`（双 gap `≥-1e-8`）、`limit_denominator(1e6)` 渲染。
包装层新增：JSON 包络、`*_latex` 字段、`status`/`reason` 中文化、`line`/`domain` 表单
参数、500 字符长度闸门、`except Exception → 400` 兜底。求解实现可直接以参考代码为
骨架，按本规格核对包络差异。

### 参考文献

[^prover-repo]: 量化调酒师. ConvexConcaveProver（凹凸不等式证明器参考实现）. GitHub 2026. [github.com/lianghuatiaojiushi/ConvexConcaveProver](https://github.com/lianghuatiaojiushi/ConvexConcaveProver)
