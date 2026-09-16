# decompose-notes — /decompose_inequality 机制规格

组合不等式（如 `pi^2+8*pi>35`、`e*pi+phi+sin(1)<11`）拆成若干基本型子证明。 端点见 api-spec.md；实现 `src/attention_calculator/decompose.py`。 判官 `bench/parity_decompose.py`：combo 87/87 字节级 + decompose 182/182 JSON 级 （transport_error 记录只需 400+error 键）。全部机制由线上探针钉死。

本文只钉站端冻结行为；exact 侧两篇见 `docs/2026-09-16-decompose-math.md`（可证构造的界分配）与 `docs/2026-09-16-decompose-exact-kinds.md`（原子型全表覆盖）。

## 输入与解析

- 输入先做 `replace(" ","").replace("\n","")` 去空白再解析。
- 文法：`e^pi` 合法；`ln2` 不带括号 → 400；含非基础常数乘积 → 400 `当前乘积证明只支持基础常数的乘积`；空/不可解析 → 400 `请输入一个只包含一个 > 或 < 的不等式`。
- `SympifyError` 是 `ValueError` 子类——须重抛 `RuntimeError` 让路由 `except Exception` 兜到 500 `组合证明生成失败，请稍后再试`；其余 ValueError → 400 原文回显。

## sympy 表示

- 函数全部用 undefined `sp.Function`，创建序 ln,sin,cos,tan,arctan,sinh,tanh,exp—— 未定义应用函数按类创建序排 canon，逐字节复刻站端项序。
- 常数全是 Symbol（e/pi/gamma/golden/catalan/gauss/varpi/zeta3），`π` 字面量→`PI`。

## 记录链 k_min（惰性）

- up 链从第一个满足 `ceil(vk)/k ≤ 6v/5` 的 k 起（ln2/ln3→4、ln5→3、pi→2、 gamma/golden→3、cos/arctan→5 全由此出）。
- lo 链从第一个 `floor(vk) ≥ 1` 的 k 起（e 在 k=1 即有界 2，`pi+e+phi>6` 实测 钉死，否定 5v/6 对称律）。
- **单原子 |coef|=1 项走原子自己的链**，其余走 |项值| 的 val_chain（k_min=1）—— flip 后落到 −1 系数的项曾误走 val_chain，是最后的批量 diff 根因。

## '<' 无乘积项 → 翻转到 '>' 域分配

全部系数取负、R→−R 按 '>' 规则分配，resid 位置仍按原 '<' 规则选 （最后负系数位，否则位 0），done 界最后翻回。

## 倒数/因子拆序

= sympy `Mul.make_args` canon 序：`pi/e`→(pi,1/e) 取分子侧、`pi^2/e`→(1/e,pi²) 取分母侧、`e^2/pi`→(e²,1/pi) 取分子侧，`reciprocal_split(first_j=0)` 直通。
