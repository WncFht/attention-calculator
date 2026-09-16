# 多 session 协作 playbook（2026-09-16 实测沉淀）

本仓常态是多个 Claude session + 每 session 多个 subagent 并行改同一棵树。以下规则都是踩坑换来的；AGENTS.md 只列了摘要，这里记全。

## git 纪律

- **写操作前先通气**：commit/branch/rebase 前 `ListAgents` 看谁在、`SendMessage` 报自己要动的文件。曾有不通气撞车的先例。
- **`git add` 只按文件名逐个加，永不 `-A`**——并行 session 的未暂存文件会被一起卷走。
- **交付即提交**：agent 交付物验收后 leader 立刻 commit——"防再被抹"（pre-commit 的 stash/restore 窗口、其他 session 的 git 操作都吃过未提交工作）。文档类交付用一次整文件 Write 而不是多次 Edit，缩短脆弱窗口。
- **绝不在共享树上跑 `pre-commit run -a` 类全量批处理**——stash 机制会竞态吞掉他人未提交编辑。
- 丢 commit 别慌：`git log --oneline --all -- <file>` + `git reflog`，并行 session 下悬空 commit 进 lost-found 是正常产物。

## 分工模式

- **文件级归属**：每个 agent 名下写死一组文件；server.py 等共享区只准 leader 动或最小成块改。
- **判官/实现分离**：replay 判官、对抗语料（judge_correct + cases_correct）由 leader 或非实现者 session 写——打分器不能出自被评分者。
- **/tmp 当中转站**：agent 交付物先落 `/tmp/`（免疫 stash 竞态），leader `cp` 进仓再提交。
- **racer 模式**：agent 停滞 >1h 无产出时，共享同一份书面 brief（如 `/tmp/decompose-brief.md`）派并行实现者，先过判官者胜。decompose 三 racer、convex 两 racer；convex 赢家是**上游源码移植**（找到原作者 `ConvexConcaveProver` 的 prove.py → 333/333）——探测推断不出的行为谜题，找上游源码比逆向快。
- **模块注入评测**：不改仓就能给候选实现打分——`sys.path` 插 `src`+`/tmp`，`sys.modules['attention_calculator.X'] = 候选模块`，然后原样跑官方判官。
- **及时关 agent**：交付物落盘确认后 TaskStop；roster 不留僵尸（同名再派会路由歧义）。追问时同名可恢复。

## 站端探测治理

- 全网only指定一个 agent 碰网络，串行 POST ≤1 req/s（授权突发 1.6s），逐字响应体追加进 `bench/data/*.jsonl`。
- 其他 agent 写探测愿望清单，leader 按名批准批次。
- 跨 session 撞探测（曾有 ~295 条来自邻 session）走 SendMessage 协调，不打速率战。

## 等待与收尾

- **不要前台 sleep 轮询**——Bash 工具 2 分钟超时，`sleep 240` 直接被杀；等长任务用 ScheduleWakeup（racer 期间 25min 档实测好用）或任务完成通知。
- commit message 里**别用反引号**——双引号 `-m "…\`exact\`…"` 会触发 shell 命令替换把消息吃掉，救回来要 --amend。
- **`pytest … | tail` 会吃掉退出码**——曾因此带着失败测试过了 pre-commit。要看尾部输出用 `pytest -x --tb=short` 或分开跑。
- 杀 dev server 用 `scripts/dev-server.sh`（按端口 `ss` 找 pid）；`pkill -f` 会匹配到调用者自己的 wrapper 进程，历史上两次把 shell 杀掉。

## 数学/工具坑位

- `mp.workdps = N` 是无效赋值（workdps 是方法不是属性，静默不改精度）——写 `mp.dps = N`（限 `from mpmath import mp` 拿到 context 对象的形态；`import mpmath as mp` 拿到模块时此句同样静默无效，须 `mpmath.mp.dps = N` 或 `mpmath.workdps(N)`）。
- `mp.quad` 精度上限 ~1e-16；要高精度参考值用收敛幂级数算，别用 quad 对 quad。
- Moment dict 省略零系数键——手搭 target 字面量要 `{k:v for k,v in t.items() if v}`，否则字典相等误判。
- golden 记录的 `parameters` 是站端参数编码，`raw_calculate` 才是逐字响应体——两个 schema 别混用。
- exact_check 复核要含 `bool(integrand)` 非空守卫：站端真有 `0 = ∫0 dx > 0` 式的"成功"记录，identity_ok=True 但命题为假。
