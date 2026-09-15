# golden 扩容采集笔记（2026-09-15）

`bench/cases.py` PARAM_SETS 扩充了 7 个最薄类型，devbox 端 `bench/harvest.py
--only golden` 续跑采集，golden.jsonl 由 2969 → **3453** 条（+484，无重复键）。

## 扩容内容

7 个常数类的 `power` 语义已实测为**系数 k**（证 k·C vs r，见
behavior-notes §4），格式层只收 `^\d+(/\d+)?$` 非负有理数——负系数一律
400，没有采集价值。每型统一扩为 `{0, 1, 2, 3, 1/2}`：

| type | 原记录 | 现记录 | 新采集成功 |
|---|---|---|---|
| gamma | 20 | 96 | 23/76 |
| golden | 22 | 90 | 35/68 |
| catalan | 20 | 90 | 36/70 |
| zeta3 | 20 | 84 | 33/64 |
| e_pi | 22 | 90 | 32/68 |
| varpi | 22 | 96 | 13/74 |
| gauss | 22 | 86 | 24/64 |

（"新采集成功" = 新增 case 中站点返回 200 的比例；其余为方向反/未找到/500，
同样是有效样本。）`true_value` 相应改为常数乘系数 x；power=0 时 v=0 走
`bounds_for` 的固定探测界 {0, 1/2, 1}。

## parity 结果（`bench/out/parity-golden-3453.jsonl`）

旧 2969 条依旧**全项逐字节一致**（site_ok_ours_fail=0、body 2969/2969、
err_match 1578/1578）。新增 484 条暴露的分歧全部集中在常数类新系数上：

- `total=484 site_ok_ours_fail=5 site_fail_ours_ok=0 param_exact=190
  body=476/484 eq=158/189 both_fail=288 err_match=286 crashes=0`
- **err_diff ×2**：`gamma 0 > 1/2`、`gamma 0 > 1`——站端严格方向预检
  404 方向反了，本方核内除零先崩 500（edge-notes §7 已记的已知分歧，
  现落入 golden）。
- **site_ok_ours_fail ×5**：`varpi 1/2<2`、`varpi 0<1/2`、`varpi 0<1`
  （本方报未找到）、`gauss 0<1/2`、`gauss 0<1`（本方 500）——见下节
  "有理余项模板"。
- **param_diff ×1**：`gauss 1/2<1`——站端走 (0,0) 余项模板，本方正常搜到
  (5,1) 真恒等式，双方各自成立但参数不同。
- **eq_diff ×31**：gamma 系数折叠 22 条 + 退化渲染 9 条，详见下。

## 新发现的站点行为

### varpi/gauss 的"初等核 + 有理余项"模板

界足够松时站端不用 `x^{4m+r}(1-x)/√(1-x⁴)` 主族，而是返回 m=n=0 的
**分子侧 √ 核加界外常数余项**：

```
varpi 1/2 < 2  →  2 - ϖ/2 = ∫₀¹ (21x³(1-x)/2)·√(1-x⁴) dx + 1/4
gauss 1/2 < 1  →  G⁻¹ - 1/2 = ∫₀¹ 5x²(1-x)·√(1-x⁴) dx + 1/3
gauss 0 < 1    →  G⁻¹ - 0   = ∫₀¹ 5x²(1-x)·√(1-x⁴) dx + 5/6
gauss 1/2 > 0  →  G/2 - 0   = ∫₀¹ (3/2-3x/2)·√(1-x⁴)/π dx + 3/16
```

参数特征：`m=n=0, bu_val=0`，`au_val/u_val` 是多项式系数，`b_val` 是
积分号外的有理余项，x 幂次仍按 r 残数（varpi `<`→x³、gauss `<`→x²、
gauss `>`→x⁰，与主族 `x^{4m+r}` 一致），`cu_val∈{1,2}` 含义待查。
 LHS 形如 `r·G⁻¹ − k`（gauss `<`）或 `k − r·ϖ` 的系数归一写法
（`gauss 1/2<1` 渲染成 `G^{-1}-1/2`，系数被搬进界里）。本方对该模板
无实现：loose 界下或耗尽搜索（未找到）或走正常 (m,n) 搜出参数不同的
另一组真解。

### power=0 的退化 200

catalan/e_pi/golden/varpi/gauss/zeta3 六型 `0 < b`（b≥0）与 `0 > 0`
均返回 200 退化恒等式，**包括 `0>0` 这类假命题**（零被积函数
`∫0·K dx > 0`）。`0 > 正数` 一律 404 方向反了。gamma 例外：
`0 < 任意` 与 `0>0` 全 500、`0 > 正数` 404（方向预检先于核）。
edge-notes §3 的探测结论现已被 golden 覆盖。

### 渲染层差异（仅方程文本，参数一致）

- **gamma 系数折叠**：站端把 k 折进第一段核——`2γ-1 = ∫2(1/(1-x)+
  1/ln x−1/2)dx`；本方渲染缺这个前置系数（22 条 eq_diff 全因此）。
- **gauss `>` 的 P 分布写法**：站端 `(3/2 − 3x/2)`（/u 逐项分配），
  本方 `(3 − 3x)`。
- **退化零多项式的化简不对称**：golden `0<0` 站端把整个被积函数化简成
  字面 `0`（`∫0 dx`），本方保留 `√(x+4)·(0)/1`；zeta3 相反——站端保留
  未化简的 `x(2x²+2)/(x²+1)·ln²x`，本方化简成 `2x·ln²x`。
- **catalan 0 vs 0**：/calculate 站端 200（`0−0C = ∫0·ln(1/x)dx`），
  但本方 /get_integral_image 对该参数组 **500**——image 端退化路径缺失。

### varpi 成功率随系数塌缩

varpi 新系数成功率 1/2→4/24、2→3/22、3→2/22，绝大多数界（含真方向）
耗尽 cap-10 报"未找到"——与 behavior-notes §3"varpi/gauss 在 d≈1e-6
即失败"的结构性约束一致：系数放大目标值后，主族核需要更高指数才够紧。

## 复跑方式

```
rsync -az --exclude .venv --exclude __pycache__ ./ devbox:~/src/attention-calculator/
ssh devbox "cd ~/src/attention-calculator && .venv/bin/python bench/harvest.py --only golden"
rsync -az devbox:~/src/attention-calculator/bench/data/golden.jsonl bench/data/
.venv/bin/python bench/parity.py bench/data/golden.jsonl
```

harvest 按键 `(type,power,comparison,rational)` 断点续跑；本轮 484 条约
11 分钟（限速 0.7s/req，成功 case 双请求）。
