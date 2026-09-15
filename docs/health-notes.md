# /health 健康计算器克隆笔记

线上行为由 `bench/data/health-probes.jsonl` 钉死：五批共 567 次 POST/GET
（serial ≥1 req/s，429 自动重试），420 个唯一 tag 全部收到有效响应。
`bench/parity_health.py` 逐条回放并做**字节级比对**：当前 420/420 一致
（`record_id` 归一化后）。本文记录规则与残存歧义。

## 路由面

| 路由 | 方法 | 行为 |
|---|---|---|
| `/health`、`/health/` | GET | 同一份中文模板（两路并存，无 308 跳转） |
| `/health/en` | GET | 独立英文模板（与 /convex/en 共用中文页的做法不同） |
| `/health/calculate` | POST | 唯一计算端点；GET/POST-with-slash → 500 |
| `/health/calculate` | OPTIONS | 200 空 body（Flask 自动应答，天然一致） |
| `/health/<其它>` | 任意 | 500 `{"error":"服务器内部错误，请稍后再试。","ok":false}` |
| `/healthxyz` | 任意 | 主站域：404 `{"error":"请求的页面不存在"}`（裸包络） |

姊妹应用域内（`/health`、`/health/*`、`/convex`、`/convex/*`）未匹配路径与
错误方法**一律 500 + `ok:false` 包络**，错误文案带句号——与主站
`服务器内部错误，请稍后再试`（无句号、`{"error"}` 裸包络）不同。
server.py 的 `in_sibling()` 按精确前缀判定，`/healthxyz` 不误伤。

## 请求校验链（成对探针钉死的顺序）

1. 传输层：非 JSON CT / 坏 JSON / 非 dict 顶层 → 400 `提交内容格式不正确。`
   （`application/json; charset=utf-8` 正常放行）。
2. `website` 蜜罐：任何真值 → `提交未通过校验。`（`""`/null 放行）。
3. `nickname`：`str()` 强转 + strip；空 → `请填写昵称。`；**先长度**
   `昵称不能超过24个字符。` 后字符集 `昵称只能包含中文、字母、数字、下划线或短横线。`
   （`a`*24+`@` 实测报长度错）。字符集为 `\w`+CJK+`-`，日文/希腊文/café 均收。
4. `sex`：`str().strip().lower()` 后须 ∈ {male,female}（`MALE`/` male ` 均收；
   `["male"]` 这类非标量经 `str()` 强转后照样枚举拒，s5:sex-list 实测 400）。
5. 数值字段统一三段式解析（label 各字段不同）：
   bool → `{label}格式不正确。`；float() 失败（含 list/dict）→ `{label}必须是数字。`；
   NaN/±Inf → `{label}必须是有限数字。`；再范围 `{label}应在 {lo}～{hi} 之间。`；
   整数字段最后查 `{label}必须是整数。`（**范围先于整数**：17.9 报范围错）。
   - 必填：`age`[18,120,int] `height_cm`[1,250] `weight_kg`[1,500]
     `pal`[1.4,2.4] 且须 ∈ {1.55,1.85,2.2}（范围内非三档 → `请选择页面提供的活动水平。`）
   - 选填（null/缺席/`""` 同义）：`waist_cm`/`hip_cm`[1,300] `resting_hr`[20,250,int]
   - `exercise_*`：**类型在场先查枚举**（`bogus` 单发也报枚举错），再成对
     `运动类型和运动时长需要一起填写，或都留空。`，再 `exercise_minutes`[1,1440,int]
   - `bp_*`：**场景在场先查枚举**，再成组 `血压测量场景、收缩压和舒张压需要一起填写，或都留空。`，
     再 `systolic_bp`[50,300,int] → `diastolic_bp`[30,200,int] →
     `sys<=dia` → `收缩压通常应高于舒张压，请检查输入。`（等号也拒；此检查在范围之后）
   - `sleep_hours`[0,24]（0 是合法值，label 为 `平均每晚睡眠时长`）

## 公式与阈值（≥8 个探针点逐一回归验证）

- 数值展示 `round()`（Python 银行家舍入）；**状态判断一律用 raw 值**
  （bmi 18.475 → 显示 18.5 但状态 `体重过低`；whtr 87.5/175=0.5 → `达到或超过`）。
- bmi = w/(h/100)²：<18.5 过低(down) / [18.5,24) 正常 / [24,28) 超重(up) / ≥28 肥胖(up)。
- target_weight = [18.5·h², 23.9·h²]（h 米，1dp）；weight_adjustment 仅
  raw<18.5（gain: 18.5h²−w）或 raw≥24（lose: w−23.9h²），1dp。
- CUN-BAE 体脂率（b=raw bmi）：男 `−44.988+0.503a+3.172b−0.026b²−0.02ba+0.00021b²a`，
  女再 `+10.689+0.181b−0.005b²`。界值 男20/25 女25/30（<warn 未达 / [warn,alarm)
  偏高区间 up / ≥alarm 报警 up）。
- **不宜旗标**：CUN-BAE 与 Deurenberg 共用 **[0,70]** 闭区间，越界 →
  `结果超出公式的合理解释范围`（无方向）。低侧 bf=−0.1 触发 0.3/0.5 不触发；
  高侧 bf=69.5 不触发、70.18 触发（s4:cb-a120-w220 vs s5:cb-hi-flag 级联实测），
  deur=69.5 不触发、70.3 触发（s5:deur-69.5 vs s4:deur-w189）。
- CUN-BAE 旗标**级联停算** fat_mass/ffm/cunningham_ree（字段 null、卡片缺席、
  专属提示替换默认对照提示）；Deurenberg 旗标不级联。
- 超龄（age>80）：状态改 `超出CUN-BAE原始18～80岁验证范围`（无方向），
  **不级联**（a81 实测 fat_mass 照出）；且与不宜旗标可叠加（a120 旗标实测：
  状态显示超龄、派生仍停、两条提示都出）。
- Deurenberg：男 `1.2b+0.23a−10.8−5.4`，女 `1.2b+0.23a−5.4`；同套 20/25·25/30 界值。
- fat_mass = bf·w/100（1dp）；ffm = w−fat_raw（1dp）；
  janma 男 `9270w/(6680+216b)` 女 `9270w/(8780+244b)`（1dp）；
  bsa = √(h·w/3600)（2dp）。
- 腰围：男 85/90、女 80/85 两档（前期/中心型）；whtr=waist/h 界值 0.5；
  whr 界值 男0.9/女0.85；rfm `64−20h/w(+12 女)`；bai `hip/h_m^1.5−18`；
  absi `waist_m/(bmi^(2/3)·√h_m)`；bri `364.2−365.5√(1−(w_m/2π)²/(h_m/2)²)`；
  conicity `w_m/(0.109√(w/h_m))`。
- Watson 体水分：男 `2.447−0.09516a+0.1074h+0.3362w`（注意实测系数 0.09516
  而非论文 0.09156），女 `−2.097+0.1069h+0.2466w`；**tbw≤0 或 tbw>w（pct>100）
  → 两字段 null + 卡片缺席 + 提示**（w30@107% 触发、w33@99.3% 不触发）。
- Nadler 血容量：男 `0.3669h_m³+0.03219w+0.6041` 女 `0.3561h_m³+0.03308w+0.1833`（2dp）。
- mifflin 男 `10w+6.25h−5a+5` 女 `−161`；harris 男 `88.362+13.397w+4.799h−5.677a`
  女 `447.593+9.247w+3.098h−4.330a`；cunningham `500+22·ffm_raw`；
  **tdee = round(mifflin_raw·pal, 0)**（1592.5×1.85→2946 证明用 raw 不用展示值）。
- 心率：max_hr raw `208−0.7a`（下游全用 raw：karvonen rhr250 出现 `217.0～203.0`
  倒挂不交换）；moderate/vigorous = raw·[0.5,0.7]/[0.7,0.85]；
  hrr = raw−rhr；karvonen = rhr+hrr_raw·[0.5,0.7]；vo2max = `15.3·raw_maxhr/rhr`。
- 蛋白 `1.4w～2.0w` 0dp；运动 kcal = `met·3.5·w/200·min` 0dp。
- 血压：sys <90/`低于90参考值`down、[90,120) `理想范围`、≥120 `高于理想范围`up；
  dia 同理 60/80。类别**高侧优先**（150/50 → 1级 而非偏低）：clinic
  3级(≥180|≥110) → 2级(≥160|≥100) → 1级(≥140|≥90) → 正常高值(≥120|≥80) →
  偏低(<90|<60) → 正常；home 界值(≥135|≥85) → 正常高值 → 偏低 → 正常。
- 睡眠：[7,8]（a<65）/ [6,7]（a≥65，含端点）；越界方向 down/up。
- 展示精度：0dp→float（mifflin/harris/cunn/tdee/心率五元组/karvonen/hrr/kcal/
  protein）；1dp（bmi/bf/deur/fat/ffm/janma/tbw/pct/rfm/bai/vo2/waist/wadj/
  target）；2dp（bri/bsa/blood）；3dp（whr/whtr/conicity）；4dp（absi）；
  int 回显（rhr/bp/exmin/age）；sleep 原样 float。

## 提示语顺序（tips / tips_en 平行数组）

身高提示（h<100）→ 体重调整 → BMI 三态 → 腰围（仅 waist_status 前期/中心型，
whtr/whr 单独越界不触发）→ CUN-BAE 块（不宜替换默认对照；超龄提示紧随其后，
二者可并存）→ Watson 不宜 → 静息心率异常 → 血压（up/down 各一条）→
睡眠（low/high 各一条）。无 exercise/whr/whtr/bf 报警类提示。

## 服务端状态与限流

- `record_id`：全站递增计数器，仅成功提交占号（失败/429 实测不消耗）。
  本站是共享状态（探针观测值 ~500），克隆用 JSON 文件计数器
  （`HEALTH_DB` 环境变量可覆盖，默认 `health-records.json`）。
- **限流未复现**：线上约 1 r/s 即 429 `{"error":"提交过于频繁，请稍后再试。","ok":false}`。
  属站点基础设施行为；复现会让 parity 回放自相残杀，故只记录不实现。
- 探针覆盖盲点（实测不可得，按相邻规则推断）：
  - CUN-BAE/Deurenberg 上界钉到公共区间 (69.5,70.18] → 取 70
    （`>`/`>=` 的端点歧义在浮点上打不到）。
  - a>80 且 bf<0 时状态是 `超龄` 还是 `不宜` 未直接命中（实现取超龄优先、
    旗标照算——a120+旗标实测派生已停，状态超龄）。
  - 非 age 字段的 `格式不正确` 已由 pal/身高/腰围/运动时长四处确认统一模式；
    `臀围/静息心率/收缩压/舒张压/睡眠` 的 bool 用例未单独探，按同一解析器实现。
