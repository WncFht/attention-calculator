# log 族实测笔记（ln_q / ln_q_square / artanh_q / arcoth_q）

实现：`src/attention_calculator/kernels/log_family.py`。数据源：`/tmp/probe_ln*.jsonl` 126 条探测 + `bench/data/golden.jsonl` 434 条（356 ln_q + 78 ln_q_square；artanh/arcoth 仅探测 16 条成功样本）。

## 核与矩空间

| type | 被积函数 | 矩空间 | P |
|---|---|---|---|
| ln_q | x^m(1-x)^n (a+bx) / (1+cx)^s | {ln q, 1} | 2 系数 |
| ln_q_square | x^m(1-x)^n (a+bx+cx²)·ln(1+cx) / (1+cx)^s | {ln²q, ln q, 1} | 3 系数 |
| artanh_q | 归约：artanh q = ½·ln q~, q~=(1+q)/(1−q) | 同 ln_q | 2 系数 |
| arcoth_q | 归约：arcoth q = ½·ln q~, q~=(q+1)/(q−1) | 同 ln_q | 2 系数 |

其中 c = q~−1（ln 型 q~=q）。矩闭式见 `log_family.iota`/`kappa` docstring。

- artanh/arcoth 的求解目标 = {ln: ±1/2, 1: ∓bound}，即**回报的 a,b 直接带 1/2**（u 不变）。c_val 复用为 q~ 字符串（"3"、"11/9"、"21"），供渲染端还原分母；ln_q 的 c_val 恒 "0"，ln_q_square 的 c_val 是 P 的二次系数。
- a_val/b_val/c_val 是 au/u 的**约分**形式（如 u=8155、au=22008 → a_val="3144/1165"）。

## s 规则

`s = max(m, n, 1)`——固定取值、不参与搜索（协调者 92/92 golden 50dps 数值拟合确认，本族全部 257 条参数记录复验一致）。

## (m,n) 遍历顺序：奇偶规则

m+n 升序 → |m−n| 升序 → **同 |m−n| 时 m+n 为奇数则 m 小者先、为偶数则 m 大者先**。即 sum=2 的序是 (1,1),(2,0),(0,2)；sum=3 是 (1,2),(2,1),(0,3),(3,0)。

- ln_q_square `3/2 < 17/100` 定案：(2,0) 可行且 (0,2) 也可行，站点取 (2,0)——engine.mn_order 的纯 m-升序在此选错。
- 该规则在全部 257 条唯一参数记录上一致（ln_q 192、ln_q_square 49、artanh+arcoth 16）。ln_q 的 46 处镜像冲突全部是奇数和，两种序不可区分；偶数和的判定证据全部来自 ln_q_square。疑似站点全局行为，其它族建议复核。

## 方向判定

- 命中"被积函数恒 ≤0"的解 → 立即 WrongDirection（引擎内）。
- 搜索耗尽（m,n ≤ 10）后由 `solve.prove` 用 float64 比较常数与界：命题为假 → WrongDirection（`要证明的式子不等号方向反了`），为真 → NoSolution（`在指数不超过10的范围内未找到{<|>}方向的解`）。
- golden 中 174+39 条 wrongdir 全部命中该路径（57 条原本误报 nosol，加 float 兜底后全部吻合）。

## 站点异常（双败计）

`ln_q_square` 在 q=5（另有 q=7）上**双向皆 500** `服务器内部错误，请稍后再试`——含 `5>2` 这种立刻可判伪的命题，属站点自身崩溃。我们的实现对真命题正常出解、对 1e-8 级紧致界给 NoSolution；parity 按双败处理（协调者已定）。

## 域检查

按 trig 族惯例，文案检查内置于 kernel `check_input`（ValueError，parity 直接比对错误文本）；顺序 = 右侧格式 → 左侧格式 → 域：

- bound<0 → `右侧有理数格式无效`；power<0 → `左侧系数格式无效`。
- ln_q / ln_q_square：q ≤ 1 → `请在ln后输入一个大于1的数`（server.domain_error 同文案 404 拦截在前）。
- artanh_q：整数或不属 (0,1) → `请在输入一个在(0,1)内的分数，本情况不支持整数`（q=1 与 q>1 同文案，含站点"在"字笔误）。
- arcoth_q：q ≤ 1 → `请在输入一个大于1的数`（同笔误）。
- **`server.domain_error` 缺 artanh_q/arcoth_q 两条规则**——kernel 的 ValueError 会被 server 映射成 400 而站点是 404，需协调者补。

## 渲染规则（257 条方程提取，byte 级验证 167/167 golden + 85/90 probes）

- 方程形：`{lhs} = \int_0^1 \frac{NUM}{DEN} \mathrm{d} x > 0`；lhs `>` → `{C} - {bound}`，`<` → `{bound} - {C}`。
- **bound 原样回显输入串**（`\dfrac{448}{250}` 不约分）。render_equation 收到的是已约分 Fraction，无法复现——parity.py 传的也是 Fraction。全族共性问题，待协调者定（需把原始数位透传）。
- C：`\ln2` / `\ln\dfrac{3}{2}`；`\ln^22` / `\ln^2\dfrac{3}{2}`（注意是 `\ln^2` 非 `\ln^{2}`）；`\mathrm{artanh}\dfrac{1}{2}`；`\mathrm{arcoth}2` / `\mathrm{arcoth}\dfrac{11}{10}`。
- DEN：c=d/e 约分，g=gcd(u,e^s)，t=e^s/g，u'=u/g。s≥2 → `{u'} \left(d x + e\right)^{s}`（u'=1 省略）；s=1 → 展开 `{u'd} x + {u'e}`（系数 1 省）。
- NUM 因子序 [x] [(1−x)] [P] [log]：
  - t 并入**第一个存在的**因子：x 部 → `2 x`；m=0,n=1 → `\left(2 - 2 x\right)`；m=n=0 → P 系数整体乘 t（`(4x²+10x+2)`）。
  - 单项 P 并入 x 部：`4 x^{2}`、`x`、`2 x^{2}`（指数 m+deg，系数 t·coef）。
  - P ≡ 1−x（au=1,bu=−1,cu=0）并入 (1−x)^{n+1}；`9−9x` 之类**不**并入。
  - 常数 P → 裸标量；值为 1 且有其它因子时**省略**（`2<1/2` 的 `(1-x)² log`）；唯一因子时为 `1`/`2`。
  - 多因子时 Add 因子加 `\left(\right)`；唯一因子裸排（`x + 50`、`101 x + 100`、`2 - x`）。
  - ` \cdot ` 规则：左因子以 `}` 结尾 且 右因子是 `\left(` 起、首字符为数字、以 `\right)` 结尾（即无外层指数）→ 加 cdot。`(1-x)^{n}` 作右因子恒不加（有 `^{n}`）；`\left(1-x\right)`、`\left(2-2x\right)`、P 括号加；P 若以 `-` 开头（`\left(- 270x²...\right)`）不加。
  - log 因子：`\log{\left(c x + 1\right)}`，sympy 风格（`\frac{x}{2} + 1`）。
- solution 串：`a = 1/15, b = 4/9`；square 加 `, c= 0`（等号后无空格——站点原样）。
- unified_form 恒 {}。

## 站点不做的事（已排除）

- 大 q **不**拆 ln（ln257、ln64 单核直解到 (7,9)/(10,10)）。
- **无 Padé 形态**——所有 ln 输出都是 (1+cx)^s 核。
- 0<q<1 的 ln_q 域拒绝（不是核内变换）。

## 待办/已标记

1. `engine.mn_order` 纯 m-升序与站点奇偶规则冲突——log 族已用本地 `mn_order_log`；其它族若出现偶数和镜像分歧应切同规则。
2. `server.domain_error` 缺 artanh_q/arcoth_q 域规则（上文文案）。
3. bound 原样渲染需要未约分数位——`render_equation` 的 bound 参数与 parity.py 链路都要透传。
4. ln_q_square q=5/7 站点 500——无法对齐，双败计。
