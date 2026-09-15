# HANDOFF — 暂停点上下文（2026-09-15）

在 devbox 上继续开发时读本文件。**Canonical 仓库：`devbox:~/src/attention-calculator`**；
本地 `~/src/attention-calculator` 只是镜像。Python 用 `.venv/bin/python`（两边都有）。
devbox sudo 密码 REDACTED（仅装系统包用；优先 uv/pip 用户级安装）。

## 项目是什么

字节级复现 https://zhuyidao.net/ —— 「注意力计算器」：输入 `常数 ⋚ 有理数`，
返回构造好的定积分恒等式证明。复刻目标是**站点行为**而非数学正确性——站端的
bug 也要原样复现（见下「站点 bug 清单」）。作者源码不公开（GitHub 上只有 skill
薄客户端，`~/src/reference/AttentionCalculator` 有镜像）。

**范围**（用户已确认扩到全站）：主站 29 型 + 姊妹应用 `/convex` + `/health`。

## Git 状态

- `master` = 最后一个全绿提交 `1a3f283`（golden 原 2969 条全字节 parity）。
- `wip/pause-2026-09-15`（当前 checkout）= 暂停时刻 6 个 agent 的在制品快照
  `882b717`（含 handoff 文档）—— kernels/render/tests/scipy 依赖/模板/探测脚本都在里面。
  续作直接在 wip 分支上继续；完成并验证后合回 master。
- **wip 快照的 parity 实测（devbox，3454 条全量）：status_match=3454/3454、
  body_exact=3454/3454、err_match=1866/1866、crashes=0、param/solution 1588/1588，
  仅剩 eq_match 1571/1588 = 17 条方程文本 diff**——kernel-edge 停前已把
  「elementary+remainder」新模板和大部分边缘修复做到接近完工，wip 比看起来更接近收敛。
- `bench/data/*.jsonl` 是 gitignored 数据资产，不进 git 但 rsync 会带；
  golden.jsonl 现 **3454 条**。

## 已验证的基线（全部实测）

| 评测 | 数字 | 入口 |
|---|---|---|
| golden parity（/calculate + /get_integral_image，含响应体逐字节） | master：原 2969/2969 全绿；**wip 快照：3454/3454 status+body，仅 17 条 eq 文本 diff** | `bench/parity.py bench/data/golden.jsonl` |
| verify（恒等式数学真值，50dps） | 1391 成功记录中 1358 真 + **33 假=站点 bug 复现**（30 trig-bias + 3 gauss-window），0 未知 | `bench/verify.py` |
| edge 重放（305 条边缘探针离线重放） | 273 match / 32 mismatch（30 条在修 + health/convex 页面） | `bench/replay_edge.py` |
| fuzz 重放（922 条随机探针离线重放） | 906 match / 16 mismatch（全部 kernel-edge 范围） | `bench/replay_fuzz.py` |
| decompose parity | 2/87——decompose.py 还没写 | `bench/parity_decompose.py` |
| 页面字节级 | `/`、`/en`、`/attention` 与站端逐字节一致 | test_client vs `bench/data/site-*.html` |
| pytest | 460 绿 + 21 skip | `.venv/bin/python -m pytest tests/` |

全量评测一条命令：`bench/run_all.sh`（同步 devbox → parity → decompose parity →
两份重放 → verify → 聚合报告）。

## 站点行为模型（已钉死的关键点）

### 协议与字节语义
- 全部响应：`json.dumps(sort_keys=True, ensure_ascii=True, separators=(",",":")) + "\n"`；
  错误体裸 `{"error": ...}`（主站）或 `{"error","ok":false}`（姊妹应用新约定）。
- `/calculate`：POST 表单；type 缺省→pi、comparison 缺省→'>'；校验序 =
  type→comparison→**右侧整组**（格式 `^\d+(/\d+)?$`→分母 0→≥10^16）→**左侧整组**→
  类型值域 404→求解。数值字段先全局删 `" "`、`"\n"` 再卡格式（`"2 2/7"`→22/7，tab 拒）。
- `/get_integral_image`：GET only；逐字段 400——int 字段 strip 后 `^-?\d+$`
  （`+3` 拒、两端空白含 tab 可），m,n∈[0,30]，u_val 下限按型（gamma=0 其余=1）；
  frac 字段同 /calculate 语法但允许负号；**coef/rational 不校验、原文回显**。
- `/decompose_inequality`：POST `problem`；空/不可解析→400 `请输入一个只包含一个 > 或 < 的不等式`；
  含非基础常数乘积→400 `当前乘积证明只支持基础常数的乘积`。
- 非 2xx 全 JSON：未知路径 404 `请求的页面不存在`；错方法→500 `服务器内部错误，请稍后再试`；
  `/en/`→404（无斜杠跳转）；`/favicon.ico`→204 空。
- 求解语义：方向假（float64 判定）报 `方向反了` 先于搜索耗尽；真但搜不到→`未找到{方向}方向的解`
  （cap：pi/e=30 其余=10，按单指数计 m,n 各自 ≤cap）；精确有理相等→404 `二者相等`（仅 Niven 点可达）。

### 站点 bug 清单（必须原样复现，verify-report.md 有全案）
- gauss `<`：m≤4 正常解，m=5 档跑**转置错误**的方程组（rows=基矩向量转置），发出数学为假的恒等式。
- varpi `<` 从不触发转置（转置解恒 a<0）。
- `>` 方向 varpi/gauss 走 sqrt 路径（cu_val=1, t=3·power gauss / t=5·bound varpi）。
- trig_pi 四型**不积分**：α 代入 Mathematica 预存闭式；(1,8) j=0 存式带
  δ(α)·(C−1) 伪项（`kernels/trig_pi.py` BIAS18_NUM）；50dps 方向预检 + defer 非正判负。
- ln_q_square q∈{5,7}：命题真→500，假→404 方向反了。
- 渲染端逐字回显（不约分、`\frac` 宏原样）；type 回显内核归一名（degree→`_pi_q`）。
- 休眠类型确认已死：`ln_pi/arcsin_q/arccos_q/Gamma_1_3_2_3/psi/zeta/erf` 服务端同样 400。

## 暂停时的在制品（6 个任务，按文件归属续作）

### 1. decompose.py（组合拆解）— 最接近完成
- 现状：规则模型在 `bench/decompose_model.py` 基本收敛——固定处理序 + 贪心链分配
  （resid 项规则、wrap-around 处理序、share=R−Σcontrib、乘积项吸附）；停前在探
  `ln5<7/4` 的判别变体（站端选 5/3 的方式还没钉死）。
- 剩余：写 `src/attention_calculator/decompose.py` + 过 `bench/parity_decompose.py`
  （combo 87 条逐字节 + decompose 182 条 JSON 级）。
- 文法情报（已验证）：`e^pi` 合法（e_pi 型）；`pi^e`、`sin(30°)`、`sin(pi/5)`、
  `ln(2)^2`、`phi^2`、`sin(1)*pi`、`ln(10)/pi` 全 400；`ln(2)` 带括号合法、`ln2` 不带 400。

### 2. kernel 边缘分歧 — ✅ 已完成（wip `b6af6b9`+`b3bb0ab`）

**验收全绿**：parity 3454/3454（eq_match 1588/1588）、replay_fuzz 922/922、
replay_edge 剩余 10 条全属其他 agent（7 float + 1 decompose + 2 姊妹页）。
钉死的机制（详见 git 历史与 test_edge.py 84 例）：'<' 退化模板在普通搜索耗尽后、
转置尝试前发射（t≥0∧b>0 严格门）；崩溃型核先做 float64 方向预检
（CONST_F/PRE_F 覆盖双曲/arctan/pi_n/e_q/gamma——q=0 不对称即由此产生）；
pi_n `wire_pair` 不约分原对；gamma coef 折叠 + u=0 `\tilde\infty` + cf/u 约分显示；
golden/zeta3 0-vs-0 sympy Mul 坍缩；选择序 `>`→sqrt_bound→`<`→plain→lt_bound→transposed。
- **新 proof family（最重要）**：varpi/gauss「elementary + remainder」模板——
  宽松界走 m=n=0、`∫poly·(1-x)·√(1-x⁴)dx + 有理余项`（余项在积分号**外**）、
  cu_val∈{1,2}、x 幂按 4k+r 余数。触发条件（vs 主搜索谁先谁后）要探测钉死。
  golden 新数据里 5 条 site_ok_ours_fail + 1 条 param_diff 属此。
- 退化解：gauss/varpi power=0 的 200 退化参数组（含 cu_val=2 变体）、
  `"00"` power、10^15 溢出垃圾、可解时退化仍触发；power=0 渲染不对称
  （golden 0<0 站端输出字面 `∫0dx`、zeta3 保留未化简式、catalan 退化参数让
  我方 image 500）。
- 双曲 q=0 崩溃不对称：cosh 0 `>`→站 404 我 500；sinh/tanh 0 `<`→站 500 我 404；
  cosh 10^15 `<`→站 500 我 404。
- gamma 渲染：系数折叠进首段 kernel（`∫2(1/(1-x)+1/lnx−1/2)dx`，golden 新数据
  22 条 eq_diff）；u_val=0 → 字面 `\tilde{\infty}`；gamma 0>正数 → 方向预检 404。
- pi_n：power∉[1,10]→500（非 400）；负系数渲染 `(\pi^{-9/8})^{8}`。
- 六条双 200 方程 diff：cos 内层 `\dfrac` 剥壳、e_q 负系数丢负号、ln_q 负参分母、
  tan_q 负因子布局、degree 族缺 `\cdot`。
- coef 回显归一：仅 `"1"`/缺省→裸常数；`"1/1"`→`1\pi`；`"01"`→`01\pi`；
  `" 1 "`→` 1 \pi`；`"2/4"`→`\dfrac{2}{4}\pi` 不约分——共享 `coef_tex` 助手进 render.py。
- 自查工具：`bench/replay_edge.py`（305 条）+ `bench/replay_fuzz.py`（922 条）
  离线重放，剩余 mismatch 应清零。

### 3. float-fidelity 保真边界
文件：`engine.py`、`solve.py`、`integrand.py`、`tests/test_fidelity.py`、
`docs/fidelity-notes.md`。
- 假说：站端 P 定号/方向判定在 float64 下进行（我们是精确有理）。
- 已知分裂：假命题+float 等值界（`pi < float64(π)`）→站 `未找到`（float 判不了
  1e-16 的非正 P）我 `方向反了`；~1e-32 超紧真命题→站误报 `方向反了` 我 `未找到`/200。
- 停前进展：zeta3 打破了单阈值模型——`>` 阈值 ≥1.46e-31 而 `<` 在 ζ3_f 即 NS，
  两个方向边界不同；正在区分「扫描中途 float 误判非正 P」vs「预检」两种机制。
- 工具：`bench/fidelity_probe*.py`（三波探针脚本）、`bench/sim_*.py`（5 种机制仿真）、
  `bench/data/fidelity-probes.jsonl`（已采集数据）。
- 注意 `solve.py` 已有 NoSolution 后的 float64 兜底（`float(C)−float(r)` 严格判号）——
  要改的是**扫描途中**的 P 定号路径。

### 4. /health 克隆 — ✅ 已完成（wip `349646a`）

**417/417 字节级 parity**（`bench/parity_health.py`，~560 次探测去重后 417 个 tag）；
568 测试绿。`health.py` = 校验链 + ~20 公式 + 36-key items + tips/tips_en +
JSON record_id 计数器（`HEALTH_DB` env）。细节全部钉死，见 `docs/health-notes.md`。
残留 4 处不可判定/假设点已记录（70 边界 ≥vs>、超龄优先级、bool 报错一致性、429 限流不复现）。
注意：`/convex/prove` 的 catch-all 目前用主站 500 文案（无句号），若 convex 探针
发现站端 500 带「。」需换 `SIBLING_INTERNAL_ERROR`（已转告 convex-probe agent）。

### 5. /convex 行为探测（进行中，数据是求解器的前提）
- 契约：`docs/sibling-apps.md` §1。语法 = Python ast（错误串直接泄漏 AST repr）。
- 待测绘：原子清单与报错全集、归一化拆边规则、凹凸分类边界、minimum 的数值方法
  （x_text ~12 位有效数字，疑似 scipy minimize_scalar——对 scipy 输出版本敏感）、
  **有理切点候选集与排除规则**（核心谜题：`e^x>=x+1` 的 x0=0 为何不被认可）、
  provided_line 语义、全部 reason 原文。
- 产出目标：`bench/data/convex-probes.jsonl`（~200-400 条）、
  `bench/data/convex-golden.jsonl`、`docs/convex-behavior.md`。
- 停前在找原作者关于凸性证明器的专栏文章（文章归档在
  `~/Desktop/obsidian/output/zhihu-mathematical/`，那里没有的话可抓知乎原文）。

### 6. /convex 求解器 — 骨架完成，等探测数据收尾
- `src/attention_calculator/convex.py` 已写：ast 解析（错误串=`ast.dump` 逐字）、
  归一化、凹凸分类、状态机+全部 reason 原文、scipy brentq 数值引擎、
  provided_line 处理、响应 dict 全字段；路由已挂（/convex/ + /convex/en +
  POST /convex/prove + /convex/static/*）。
- 样例精度：`minimum.value` 逐位一致，`x` 差 1 ulp——等 convex-probe 钉
  bracket/xtol/候选集后对齐。全部待定点收在 `ProbeConfig` +
  `convergent_candidates`（候选集假说：CF 渐近分数，样例 17/30 吻合）。
- 姊妹应用共享 `in_sibling()` 错误包络（{"error","ok":false}）。

## 协作约定（照此执行过的）

- 每个 agent 只动自己名下的文件；server.py 是共享区，改动保持最小且成块。
- 探测站端限速 ≤1 req/s；数据文件一律存 `bench/data/*.jsonl` 带原始响应体。
- 判分工具独立于我写的实现（parity_decompose/replay_edge/replay_fuzz 都是
  leader 写的离线重放器，防止 agent 自己给自己判分）。
- 文档只记核实过的事实；站点行为一律以实测为准，文章/文档只做线索。
- 每完成一块：跑该域 parity → pytest → commit（Conventional Commits）→
  rsync 回本地。
