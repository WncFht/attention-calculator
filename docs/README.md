# docs/ 索引

项目分两阶段：site 字节级复刻（冻结于 `v1.0.0-site-parity`）与 mode=exact 数学正确性（当前主线）。文档按用途分组；日期前缀文件是时间点文档，其余是持续维护的活文档。

## 活规格（改动实现前先读/先改）

- `kernel-spec.md` — 29 个站端类型 + 26 个 exact-only 型的核函数/矩空间/搜索顺序/渲染怪癖，数学层唯一事实源
- `api-spec.md` — 线上 API 协议实测（含 mode=exact 本地扩展节）
- `behavior-notes.md` — 站端行为规格（power 语义、校验顺序、错误分层）
- `fidelity-notes.md` — 方向判定的 float64 保真机制
- `decompose-notes.md` — /decompose_inequality 拆解机制
- `convex-behavior.md` + `convex-doc-audit.md` — /convex 行为规范与独立复核
- `sibling-apps.md` + `health-notes.md` — 姊妹应用侦察与 /health 克隆笔记
- `log-notes.md` — log 族实测（s=max(m,n,1) 等）
- `verify-notes.md` — verify.py 方法学（mpmath/acb/Sturm 三层；含重建伪影重归因补记）
- `edge-notes.md` / `fuzz-notes.md` — 边缘输入探测（376 探针）与随机模糊（922 条）记录
- `article-audit.md` / `surface-diff.md` — 作者专栏文章语料审计、站点表面审计

## exact 阶段（W0–W5）

- `2026-09-16-math-correctness-plan.md` — 阶段方案与工作流划分（头部有当日落地速记）
- `2026-09-16-exact-status.md` — 当前状态汇总：各 W 落地情况、判官结果、在途项
- `2026-09-16-certificate-spec.md` — 证书 schema（含 Padé 变体）
- `2026-09-16-pade-notes.md` — Padé 第二证明器推导与覆盖
- `2026-09-16-ln-cube-derivation.md` — ln_q_cube 矩空间推导
- `2026-09-16-ln-quad-impl-notes.md` — ln_q_quad 落地：ln³ 核矩空间、四次 P 的 Sturm 非负判据、实测深度
- `2026-09-16-agm-impl-notes.md` — AGM 区间第二证明器（gauss/varpi）：isqrt 包络、cert schema、实测地板
- `2026-09-16-euler-gamma-impl-notes.md` — Euler–Maclaurin γ 第二证明器：EM 界证明、ln2 子证 DAG、Bernstein 符号引理重证
- `2026-09-16-beta-even-impl-notes.md` — beta4/beta6 落地：Dirichlet β 偶数点、zeta_odd 同构镜像
- `2026-09-16-dixon-impl-notes.md` — pi3/pi3_u/pi3_a 落地：Γ(1/3) 格三常数（Dixon 椭圆周期 π₃）、六配置双向覆盖
- `2026-09-16-li2-impl-notes.md` — li2_q 落地：三符号矩空间、q<0 凹性引理更正、交错级数替代 lifted 族
- `2026-09-16-pi-sqrt2-impl-notes.md` — pi_sqrt2 落地：lemniscate 格点末个单核型（S=π√2）、两方向各一核
- `2026-09-16-sicin-impl-notes.md` — si_q/cin_q 落地：一族两型同机、四维 span 的三次非负判据
- `2026-09-16-trigamma-impl-notes.md` — psi1_q 落地：对称望远镜核双向二次收敛、N 不进参数、一维 m 扫描
- `2026-09-16-decompose-exact-kinds.md` — decompose_exact 原子表全型覆盖 + arg<0 预存 bug 修复记录
- `2026-09-16-decompose-math.md` — exact 模式 decompose（可证构造的界分配）
- `2026-09-16-lean4-integration.md` — Lean4 接入调研：verified-checker 路线选型、逐族矩引理工作量、W8 分阶段建议
- `2026-09-16-w3-research-*.md` — 七篇新型可行性调研（erf、Γ(1/4)、Γ(1/3)、Glaisher、Li₂、Si/Cin、trigamma）；**各篇头部有落地状态行**，可行未注册与已否决一眼可查

## 里程碑与过程记录

- `2026-09-16-site-parity-status.md` — site 阶段终点：判官全绿基线 + 失真面权威盘点（46+6）
- `2026-09-16-dev-history.md` — 两日开发史复盘：双机时间线、人类指令全录、会话/subagent 名册、事故与观察
- `HANDOFF.md` — site 阶段完成记录（架构/分工/遗留）
- `ROADMAP.md` — site 阶段路线图与实测行为记录（含 tie-break 更正）
- `verify-report.md` — golden 全量数值验证报告（含重建伪影重归因）
- `harvest-notes.md` — golden 2969→3454 扩容采集记录
- `collab-notes.md` — 多 session 协作 playbook（文件归属、判官/实现分离、探针治理、踩坑清单）
