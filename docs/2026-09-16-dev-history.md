# attention-calculator 开发史（2026-09-15 → 2026-09-16）

本文是对本项目两天开发过程的事后复盘，由 attention-calculator-ba 会话基于 devbox 上的全部会话 transcript、git 历史、mba 存档 bundle 与仓库文档重建。所有时间均为 UTC+8。

## 一句话摘要

两天、两台机器、218 个提交、8 个 Claude 会话、59 个 subagent：9-15 全天在 Mac（mba）从零复刻 zhuyidao.net「注意力计算器」，当晚迁移到 devbox 通宵收口；9-16 上午字节级复刻冻结为 `v1.0.0-site-parity`；9-16 中午用户一句话把目标从「复刻站点行为」切换到「数学正确性」，当天下午到晚间以 leader + 判官 + 实现 agent 的编队把 exact 系统扩到 26 个新类型、四条第二证明器路径、全量判官 12898 条记录零 BUG。

## 证据边界

- devbox 侧 8 个会话的 transcript 完整保留在 `~/.claude/projects/-home-wncfht-src-attention-calculator/`，包括 59 个 subagent 的 `.meta.json`，本文的指令清单、对话与过程细节主要来自它们。
- mba 侧的会话 transcript 已随仓库删除不可考。Mac 侧工作的权威证据是 `~/archive/mba-src-2026-09-16/` 存档（`attn-calc-2026-09-16.bundle` 全量 git 历史 62 提交 + 脏工作区 patch + untracked tar，制作于 9-16 08:37）以及 HANDOFF/ROADMAP 等过程文档。9-15 白天的逐会话细节（Mac 上有几个会话、各自分工）只能从提交流与文档推断。
- 时间线数字来自 git 提交时间戳（UTC+8）与 transcript 内 UTC 时间戳换算。

## 全景时间线

| 时段 | 机器 | 事件 |
|---|---|---|
| 9-15 09:29 | Mac | 首个提交 `41d8c73`「project skeleton with kernel/api/bench specs」 |
| 9-15 白天 | Mac | 引擎/8 族 kernel/探测/benchmark/webapp 全线推进，约 48 提交 |
| 9-15 ~21:00 | Mac | 语料扩到 3454 条、范围扩到姊妹应用；`cd73188`「wip: in-flight agent work paused for devbox move」等三个暂停/交接提交（21:00–21:09），HANDOFF.md 落成 |
| 9-15 21:08 | devbox | leader 会话 5e（f004866f）在 devbox 启动：「HANDOFF.md 现在请你继续做，多用 subagents」 |
| 9-15 21:25/21:27 | devbox | toolchain 会话 12（916ce5e4）、side-api 排障会话 09（b79fce11）先后启动 |
| 9-15 21:45–22:29 | devbox | leader 侧 kernel-edge 修复、/health 克隆、/convex 探测规格落盘；约 22:29 仓库回同步到 Mac |
| 9-15 22:30 | Mac | `b34c3c2`「docs: mark convex probing complete in handoff」——Mac 上最后一个提交，未回传 devbox，是全部历史中唯一未进 master 的提交 |
| 9-15 23:39–9-16 05:02 | devbox | 通宵收口：12 装 format/lint 工具链（23:39–23:48）并在 04:38–05:02 逐文件刷完全仓 ruff；5e 落地 convex 求解器移植（00:09 `4675217`）、decompose（04:32 `6e9b850`）等 |
| 9-16 07:50–10:50 | devbox | 用户连环验收：进度→subagent 全面审计→pipeline 讲解→演示→前端修复→优化审计→合并 wip→master（10:30 `4527520`）→markdown 一段一行（10:41 `f64c230`）→规则写入全局 AGENTS.md |
| 9-16 08:37 | devbox | 用户把 Mac 仓库存档（bundle+dirty patch+untracked tar）后删除 Mac 侧仓库，devbox 成为唯一 canonical |
| 9-16 10:53 | devbox | **转型指令**（全文见下）：打 tag、写报告、转向数学正确性、疯狂扩展 |
| 9-16 11:21 | devbox | annotated tag `v1.0.0-site-parity` 打在 `f64c230`（156 提交）；11:29 状态报告+W0–W5 方案提交（`592ec3d`） |
| 9-16 11:51–12:01 | devbox | W0 骨架：prove_exact + certified_cmp + exact_check 契约 + exp 参考实现（`6742443` 等） |
| 9-16 12:37 | devbox | a9 会话（8b04178a）被用户派来「帮他们一起推进进度」，成为独立判官/调研线 |
| 9-16 13:13–16:55 | devbox | 五波 agent 交付：W1 失真修复 → a9 判官首扫抓 62 BUG → W3 新类型梯队 → W4 Padé → W5 证书 → W7 第二证明器（AGM/Euler–Maclaurin/ln4/decomp-ext） |
| 9-16 15:13 | devbox | 71 会话（179eeab9）做全仓文档整理与工具沉淀 |
| 9-16 17:47–18:2x | devbox | 7e 会话（51fc4ded）起 web 给用户体验；实测发现站端 float64 方向预检误判窗 (fl(C), C)，双模式钉子入档（`09eb457`） |
| 9-16 17:58–19:10 | devbox | f0 会话（ca107b99）审计知乎文章 → 讨论项目意义 → lean4 接入调研（3 个 lean-* agent），结论「先不加，写文档」（`d97bea5`，19:10，当日最后提交） |
| 9-16 18:40–18:52 | devbox | a9 提案的统一 q=0 闸门落地（`7d240c9`）；a9 终报 12898 条零 BUG/FLAG 收进状态文档（`f8e7832`） |
| 9-16 ~19:2x | devbox | 本会话（ba）受用户委托撰写本文 |

## 人类指令全录

以下是两天内人类用户发出的全部消息（按会话分组，逐字摘录；时间 UTC+8）。Mac 侧 transcript 已删，9-15 白天 Mac 会话中的指令不可考——但从提交内容看主要是复刻工作的推进指令，21:00 前后的暂停/交接指令产生了 HANDOFF.md。

### 主线 leader 会话 5e（f004866f，devbox，16 次上下文压缩）

| 时间 | 消息 |
|---|---|
| 9-15 21:08 | `/home/wncfht/src/attention-calculator/docs/HANDOFF.md 现在请你继续做,多用 subagents` |
| 9-16 00:05 | `continue` |
| 9-16 07:50 | `现在进度怎么样` |
| 9-16 07:52 | `我希望你不只看 HANDOFF,我希望你用多个 subagents 全面看看整个仓库的目标是否达成了` |
| 9-16 08:00 | `continue` |
| 9-16 08:19 | `搞完给我讲一下现在整个项目的情况,整个 pipeline 是怎么样的,原理是怎么样的,有什么功能` |
| 9-16 08:22 | `你能不能给我演示一下看看` |
| 9-16 08:34 | `你把 mathjax 啥的配好,我看现在前端还有很多问题` |
| 9-16 09:06 | `还不错,我希望你现在给我全面看看,使用多个 subagents,看看整个项目还有什么可以优化改进的地方,以及你整理整理我们这个项目` |
| 9-16 10:34 | `1. 你合并一下,然后讲讲现在这个项目是不是已经彻底做完了?要不要重构之类的,适合发布吗,我们的代码优雅程度怎么样 2. 把 markdown 文档全部修复一下,一段应该在一行内` |
| 9-16 10:43 | `请你把这个一行一段的 markdown 要求给我补充到全局的 AGENTS.md 以及 mba 上面的` |
| 9-16 10:48 | `为啥 docs 下面的没有过 format ?` |
| 9-16 10:50 | `讲讲我们现在这个整个仓库在做什么,然后架构是怎么样的,所有的情况都能覆盖吗?` |
| 9-16 10:53 | `笑死我了,我觉得你应该要数学正确性而不是站点行为本身,现在请你给我们现在这个版本打个 tag,然后写个详细的报告记录现在的情况,后面我们转换目标到数学正确性上面,然后尝试疯狂扩展我们整个系统,你可以用 subagents,请你做完以后给我详细讲讲,我们应该怎么扩展整个系统` |
| 9-16 17:19 | `讲讲现在是什么进度,是不是已经做的差不多了?` |
| 9-16 17:34 | `一定要等他吗?你把整个仓库的内容看看,总结一下今天做了什么,以及后面的发展计划是什么` |

### toolchain 会话 12（916ce5e4，devbox）

| 时间 | 消息 |
|---|---|
| 9-15 ~22:20 | `/setup-toolchain 请你帮我把这个仓库的 toolchain 配置好` |

### side-api 排障会话 09（b79fce11，devbox；与本仓代码无关但同机同时段）

| 时间 | 消息 |
|---|---|
| 9-15 23:45 | `你帮我去详细看看 mba 上面的 side-api,为啥在 mba 上面的 claude code 好像都连不上他们本地的那个 side-api 啊?` |
| 后续 | `感觉你可以把 A 和 B 都做一下` / `continue` / `还是有问题 ⏺ API Error: 502 ... [::1]:3003,请你做好验证,同时去读一下他们的日志等等`（最终定位 VS Code autoForward 抢占 127.0.0.1:3003，改走 [::1]） |

### 判官/调研会话 a9（8b04178a，devbox）

| 时间 | 消息 |
|---|---|
| 9-16 12:37 | `f004866f-59d0-4cad-b012-ee9f31771979 话说你能不能帮他们一起来推进进度,你可以用多个 subagents` |

### 文档整理会话 71（179eeab9，devbox）

| 时间 | 消息 |
|---|---|
| 9-16 14:03 | `看看这个项目,给我整理更新所有文档,然后我希望你仔细看看相关的所有聊天记录…看看我们能不能总结出有用的脚本或者小工具等等,你可以用多个 subagents 来做这个事情` |

### 演示会话 7e（51fc4ded，devbox）

| 时间 | 消息 |
|---|---|
| 9-16 17:47 | `把这个项目起个 web 让我体验一下,然后把你让我体验的输入都给我,我来看看情况` |
| 9-16 18:01 | `你能不能给我一个非常非常紧的不等式,然后不要太复杂` |
| 9-16 18:12 | `你去跟这个项目的负责人 f004866f-... 讨论一下这个事情` |

### 审计/调研会话 f0（ca107b99，devbox）

| 时间 | 消息 |
|---|---|
| 9-16 17:58 | 知乎文章审计请求（逐条核对站点文章内容是否都被复现） |
| 9-16 18:09 | `我们在做的这个实际的用处是什么?能不能搞出什么更加高维的认知或者直觉之类的?` |
| 9-16 18:11 | `给我讲讲 lean4 的原理是啥,我们这个要接入吗?有啥用?` |
| 9-16 18:14 | `详细调研一下, lean4 具体是怎么接入的?` |
| 9-16 18:53 | `讲一下 mathlib lrat_proof 63MB 证书 OOM 到 95.7GB 这个是什么意思` |
| 9-16 19:06 | `所以个这个 lean4 对我们有啥用吗` |
| 9-16 19:09 | `算了,你先别加了,在文档里面写一下` |

### 本会话 ba（3427f38b）

| 时间 | 消息 |
|---|---|
| 9-16 ~19:2x | 委托撰写本文（全面阅读两台机器上的开发历史与会话历史，整理人类指令、开发流程、结果、时间线、subagent 使用情况，可用 subagents 调研） |

## 会话与 subagent 名册

devbox 侧共 8 个会话，显示名 ↔ sid 对应关系由 transcript 自报与消息内容交叉确认：

| 显示名 | sid 前缀 | 角色 | subagent 数 | 备注 |
|---|---|---|---|---|
| attention-calculator-5e | f004866f | leader/主线：HANDOFF 收口 → 上午验收响应 → exact 阶段总指挥 | 43 | 26.7MB transcript，16 次压缩；所有 git 提交与注册表收口都归它 |
| attention-calculator-12 | 916ce5e4 | 工具链：pre-commit/format/lint 全套、CI、AGENTS.md 布局 | 1 | 凌晨 04:38–05:02 的全仓 ruff 逐文件提交出自它手 |
| attention-calculator-09 | b79fce11 | side-api 排障（Mac 侧服务，与本仓无关） | 0 | 产出进了用户全局记忆而非本仓 |
| attention-calculator-71 | 179eeab9 | 文档整理/工具沉淀 | 5 | 接手 README/AGENTS/kernel-spec/api-spec 刷新 |
| attention-calculator-a9 | 8b04178a | 独立判官 + W3 调研机队 | 7 | judge_correct 语料/判官全归它；w3 调研笔记（erf/gamma-quarter/li2/sicin/psi 等）是它的子代理产出 |
| attention-calculator-7e | 51fc4ded | web 体验/demo | 0 | 实测贡献：(fl(C), C) 方向预检误判窗的实证 |
| attention-calculator-f0 | ca107b99 | 知乎文审计 + lean4 调研 | 3 | lean-api/lean-eng/lean-nonneg 三个调研 agent；产出 `2026-09-16-lean4-integration.md` |
| attention-calculator-ba | 3427f38b | 本文作者：开发史重建 | 0 | — |

leader 的 43 个 subagent 按波次：site 阶段（decompose/decompose2/decompose3、convex-skel、convex-probe、convex2、health、fidelity、kernels、replay-capture-gap、src-refactor）→ 上午审计波（audit-src-quality、audit-packaging、audit-judge-suite、audit-golden-parity、audit-docs、audit-divergences、audit-coverage、audit-conformance、audit-bench-tests）→ exact 阶段（exact-quadlog/log/trig/hyp/beta 五型失真修复、w3-zeta/lnpow/arcsin/invhyp/gauss-erf/pisqrt2/dixon/li2/trigamma/sicin/decomp、w4-pade、w5-cert、w7-agm/euler/ln4/decomp-ext/beta）。a9 的 7 个是调研 agent（w3-erf/w3-gamma-quarter/w3-li2/w3-psi/w3-sicin 等），产出 impl-notes 供 leader 派实现 agent。

## 协作机制（实测成型的规则）

- **leader 制**：5e 是唯一提交者——agent 只交文件和报告，wiring（TYPES/EXACT_TYPES/FAMILY/certified_cmp/渲染/证书 TeX/域检六处注册点）、测试、提交、关 agent 全归 leader。高峰期它保持 5–6 个实现 agent 在编 + a9 的判官/调研线并行。
- **判官/实现分离**：判官脚本（parity/replay/judge_correct）与实现者不同人，a9 的 judge_correct 用 `constant_mpf` 区间比较独立裁决 ground truth，不信仰站端标签也不信仰实现。
- **交付即关**：交付物落盘确认后立即关 agent，roster 不留僵尸；leader 多次主动清 roster。
- **文件级归属**：每个 agent 只动自己名下文件；共享区（server.py、solve.py、注册表）归 leader。跨会话消息用 SendMessage，git 写操作前互相通气。
- **先规格后码**：kernel-spec.md 是唯一事实源；agent 报告必须附 spec 草稿节，leader 合并时统一口径。
- **质检文化**：每条 emitted proof 先过 exact_check 自检（ℚ 字典相等，全称验证非采样），失败直接 InternalError——假证明结构性发不出去。

## 事故与坑（值得记住的）

- **pre-commit stash 窗口竞态**（9-15 深夜）：12 的格式化提交触发 pre-commit 的 stash 机制，反复吞掉 5e 侧 agent 的未暂存工作；5e 的应对是「在对方下一次 wipe 前几秒内抢提交」。教训沉淀为：共享工作树不跑 `pre-commit run -a`，中间产物放 /tmp。
- **vscode autocorrect formatOnSave**：编辑器侧 autocorrect 独立于 pre-commit 链生效，会改 docs/ 里逐字引用站端字符串的 inline code span（内容即证据，改了就失真）——docs/ 整目录划出全部 formatter + 编辑器关 formatOnSave。
- **压缩幻觉**：health agent 的压缩摘要声称 s5 探测已完成，磁盘上实际没落盘——agent 自己发现后补跑。教训：压缩摘要不可信，以磁盘为准。
- **Python 版本差**：站端跑旧 CPython——float 求和用朴素 sum（非 3.12+ 的 Neumaier），`ast.dump` 需 `show_empty=True` 才逐字节对齐；两处都是为 parity 专门钉的。
- **mpmath 精度陷阱**：`import mpmath as mp; mp.dps = N` 是模块属性静默无操作（真实精度仍 15），须 `mp.mp.dps` 或 workdps；`Fraction * mpf` 会塌缩成 float64——decompose_exact 的 slack 曾因此失真（`e00e5cd`）。
- **判官抓真 bug**：a9 首扫 8830 条抓出 62 个 BUG（双曲 q<0 未奇偶归约 41 条 WD-on-true、power=0 崩溃 21 条），同日修复转阴（`460c5b1`）；后续又把「checker 零值键过滤」「nonneg 恒零守卫」打成全族契约。
- **跨会话撞车**：5e 与 a9 在 kernel-spec、bench/cases.py 等文件上有明确的分区通气；7e 的演示实例（:8081）与 leader 的测试并存无冲突。
- **测试也犯数学错**：至少两次 agent/leader 写的测试用错了界（把 WrongDirection 的正确答案当成 bug 修），靠高精度数值核对纠正——判官是机器，判题人也得被机器判。

## 最终结果

### 阶段一：站点复刻（冻结于 `v1.0.0-site-parity`，156 提交）

全判官零分歧：golden 3454/3454（eq 1588/1588、响应体逐字节）、decompose 87+182、health 420、convex 333、edge 376、fuzz 922、fidelity 228 + probes 160、页面字节级、619 测试绿。复现含全部已知站端 bug：gauss/varpi '<' 转置解、trig_pi (1,8) 损坏预存公式、ln_q_square q∈{5,7} 500 分裂、渲染逐字回显。权威失真面：1588 条成功记录中 46 条假恒等式 + 10 条零被积函数空洞证明，全部按簇归因入档。

### 阶段二：数学正确性（9-16 一天，62 提交）

- **26 个 EXACT_TYPES**：zeta5/7/9/11、beta4/6/8/10、ln³/ln⁴、arcsin/arsinh、gaussint/dawson/erfiint、pi_sqrt2、pi3/pi3_u/pi3_a（Dixon 格点）、li2_q、psi1_q、si/cin、gamma14/34/12（composite DAG）。
- **四条第二证明器路径**：Padé（ln/arctan）、AGM 区间包络（gauss ~1e-2000 间隙、varpi ~1e-16）、Euler–Maclaurin（gamma，实测地板 1e-200 双向）、composite DAG（Γ 特殊值 + varpi÷AGM）。三弱型覆盖率被直接抬穿：varpi 26→108、gauss 28→105、gamma 32→105，36 条站端假证明簇（真命题）全数证出。
- **证书体系**：classic/pade/agm/composite 四种形状在 verify_cert CLI、verify_response、/calculate、decompose_exact 五处闭环，均校验证书内嵌命题与请求一致。
- **判官终报（a9，12898 条最终码全量重验）**：零 BUG/FLAG，4547 条 emitted 全过 exact_check 零容忍线，145 条 composite 端到端复验全过。覆盖率强（≥80%）14 型（gamma12 96% 领衔），残余 unsolved-true 全是 '<' 紧界的诚实搜索缺口（pi_sqrt2 24%、pi3 26%、gamma14 27%、gamma34 31%、psi1_q 42% 为候选第二证法面）。
- **测试套件**：约 1724 绿；site 路径全程零回归。

### 明确否决/遗留

- 否决：erf 本体（√π 障碍）、Γ(1/4)/Γ(1/3) 字面幂（格点奇偶锁）、Glaisher/lnA（矩空间符号数超 P 系数）、arcosh_q（span 无 "1" 方向）、Lean4 接入（f0 调研后用户拍板「先别加了，写文档」）。
- 遗留：ln_q_square q=13 预言崩溃未探（行列式根因已定位 q∈{5,7} 同型）；β(12)+/ζ(13)+/ln⁵ 是一行表项的同构扩展，刻意停在对称截面；五弱型的第二证明器是下一波候选。

## 观察

1. **这个项目本质上是在用工程方法做数学**：站端原作者用 Mathematica 预计算查表；本项目用闭式矩现算。leader 把它总结为「矩机器」——唯一需要数学的一步是找核使矩落在有限维 ℚ-空间，其余全是机械工程。这个抽象是 9-16 下午能安全「疯狂扩展」的根因：每型正确性由 exact_check 逐条全称保证，扩张不会把错误混进来。
2. **人类指令密度低但方向感强**：两天约 34 条消息，模式是「委派—验收—转向」。关键转向两次：9-15 21:00 前后决定迁机（产生 HANDOFF 制），9-16 10:53 一句话换目标（「笑死我了…要数学正确性而不是站点行为本身」）。每条指令都改变了资源分配，没有微管理。
3. **「全面审计」是用户的标准验收动作**：上午 07:52 的「不只看 HANDOFF，用多个 subagents 全面看看」直接孵化出 9 个 audit-* agent；09:06 的「还有什么可以优化改进」又孵化一波。每次用户验收都伴随一次全仓并发审计，这与「求全不妥协」的全局偏好一致。
4. **判官/实现分离产出真实价值**：a9 不是附属——它首轮就抓 62 个真 BUG，后面每波扩张都靠它的语料收尾；它还能反过来纠正 leader 的归因（转置簇的命题真假之辨）。两个会话之间出现了真实的「同行评审」动态。
5. **转型的执行速度快得反常**：从 10:53 转型指令到 11:21 打 tag、11:29 报告+方案、11:51 W0 骨架落地——不到一小时完成阶段切换，因为 Moment dict 结构让精确验证几乎免费。下午 26 型落地靠的是「调研 agent 出 impl-note → leader 验证数学 → 实现 agent → leader 收口」的流水线。
6. **多会话协作的真实成本**：stash 竞态、压缩幻觉、口径漂移（q=0 分类前后两版）、僵尸 roster、交付回声消息——都被实战踩过并沉淀成 collab-notes 的规则。这套规则本身是项目的重要产出。
7. **证据链完整是刻意的**：每个数字都落到可重跑的判官上（bench/README 声明「benchmark 是一切结论的依据」），每篇调研笔记头部有状态行，每个「不可能」都有结构性论证。这份文档能写成，本身就受益于这个纪律。
8. **一些有趣的边角**：唯一丢失的提交 `b34c3c2` 只是个 HANDOFF 状态行，实际信息零损失；leader 曾因 commit message 里的反引号被 shell 吞掉一段而返工（collab-notes 从此禁用反引号）；7e 演示会话顺手实测出站端 float64 预检的 (fl(C), C) 误判窗——这是「用户体验」指令意外产出的站点行为新发现。
