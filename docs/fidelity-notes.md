# 方向判定的 float64 保真机制（已钉死）

`solve.prove` 在进核前复现站端的方向预检：站端用 float64 求值命题常数 `c`，
再把有理界 `bound`（精确 Fraction）与 `c` 做 **精确比较**——Python 的
`Fraction > float` 语义把 bound 与 float 的精确 dyadic 值比较，bound 侧不舍入。

- `'>'` 命题：`bound > c` → 404 `要证明的式子不等号方向反了`，不进搜索。
- `'<'` 命题：`bound < c` → 同上；`bound == c` 越过预检进扫描。
- 扫描途中命中非正解：`'>'` 仍报方向反了；`'<'` 报 `未找到<方向的解`
  （站端 '<' 的方向判定只在预检发生，非正 P 在扫描里等同搜不到）。

这统一解释了两类保真分歧：

1. **超紧真命题误判方向**（cap:*-tight 系列）：c 低于真值的型
   （pi/e/sin_q/ln_q/e_q/ln_q_square——`fl(C)<C`），真命题界 `bound<C` 但
   `bound>c` 时站端直接判反。float 判定点距真值仅 ~1e-32 的界 float 上
   与 c 等值，故只有精确比较能解释。
2. **float 等值界的 '<' 假命题**（equal:pi-float-lt、e-float-lt）：
   `bound == c` 不触发 `<` 预检，扫描途中非正解按站端语义报未找到。

## 各型常数 c

取站端字面表达式的 float64 结果，不一定是正确舍入值：

| kind | c 表达式 | 备注 |
|---|---|---|
| pi / e / catalan / gamma / golden / varpi / gauss / e_pi / zeta3 | `float(q)·K` | K 见下 |
| e_q | `math.exp(float(q))` | |
| pi_n | `math.pi ** float(q)` | |
| arctan_q / arccot_q | `math.atan(f)` / `math.atan(1/f)` | |
| sinh/cosh/tanh/coth_q | `math.sinh/cosh/tanh(f)` / `1/math.tanh(f)` | |
| sin_q / cos_q / tan_q / cot_q | 同名 `math.*`；`cot_q` 为 `1/tan(f)` | 值域外不预检 |
| ln_q | `math.log(f)`，仅 `q>1` | 域校验先于方向 |
| ln_q_square | `math.log(f)**2`（二次舍入），仅 `q>1` 且 `q∉{5,7}` | q∈{5,7} 走核内特判 |

`K`：`math.pi`、`math.e`、正确舍入的 `fl(Catalan)`、`fl(ζ(3))`、`fl(e^π)`
（注意 **不是** `math.exp(math.pi)`，那个低 1 ulp）、Euler 常数字面量
`0.5772156649015329`（与 kernels.gamma 一致）、`(1+√5)/2`、以及 varpi/gauss
的 mpmath 值。`math.e**3` 之类的幂表达式会再差 1 ulp，必须用 `math.exp(q)`。

## 例外与边界

- **zeta3 `'>'` 是唯一不走精确预检的方向**：站端实为 `float(bound) > c`。
  证据：fid4:zeta-gt-fe* 系列界在 `(Cf, Cf+1ulp)` 内站端报未找到；
  `Cf+1ulp` 起报方向反了。由 NoSolution 后的 float 差兜底覆盖。
- `power=0` 系数型：`c=0`，`'>'` 仅 `bound>0` 触发（0>1 → 方向反了），
  `0>0`/`0<0` 进核给退化 ∫0dx 证明（200）。
- 负界、格式非法、分子分母 ≥10^16 在 server 表层就被 400 挡掉，到不了预检；
  prove() 内部对 `bound<0` 也不做预检，交给核内 check_input 报右侧格式错。
- artanh_q/arcoth_q/trig_pi 四型的方向判定路径未探测区分，保持核内现状。

## 证据锚点（bench/data/{edge,fidelity}-probes.jsonl）

- `'>'` 预检判反：cap:pi-gt-tight、cap:e-gt-tight、cap:e_q-tight、
  cap:sin_q-tight、cap:ln_q-tight、cap:ln_qsq-tight（站端全报方向反了）。
- `'>'` 反向对照（c 上溢，真紧界放行）：cap:cos_q-tight、cap:atan3-tight
  （站端报未找到>方向的解）。
- `'<'` 等值放行：equal:pi-float-lt、equal:e-float-lt、fid3:cat-lt-+0、
  fid3:zeta-lt-+0/+1/+2（站端全报未找到）；fid4:*-<-1 系列（Cf−1ulp）判反。
- zeta3 `'>'` float 判定：fid3:zeta-gt-+0（Cf，未找到）vs fid3:zeta-gt-+1/+2
  （方向反了）vs fid4:zeta-gt-fe*（(Cf,Cf+1ulp) 内，未找到）。
