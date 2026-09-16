# 数学正确性扩展方案（2026-09-16）

基线：tag `v1.0.0-site-parity`（f64c230），站点字节级复刻全绿。新目标：**产出数学上为真的恒等式证明**——恒等式等号在 ℚ 上精确成立、被积函数定号有精确证书、方向判定有认证依据；并在此基础上大规模扩展类型与命题覆盖。站点复刻行为整体冻结为 compat 层。

## 总体判断

现有架构对这个转向几乎是准备好的：`Moment` 就是 ℚ-向量（`{常数符号: Fraction}`），每个核已有真矩生成器，`solve_moment` 在 ℚ 上求解。数学验证 = `combine(coeffs, true_basis(m,n)) == target` 的字典相等，零新数学依赖。真正的工作在三层：(1) 每族暴露真 basis/目标构造供校验复用；(2) 四个失真簇的正确路径；(3) 新类型的矩空间推导。瓶颈只有第 (3) 层要研究，其余是工程。

## 双模式架构

`solve.prove(kind, power, comp, rational, exact=False)` 增加模式旗标；`mode=site` 行为逐字节冻结（全部现有判官钉死它），`mode=exact` 走正确路径。/calculate 加 `mode` 参数，默认本期保持 site（判官不动），正确性 benchmark 达标后再考虑翻默认。exact 模式响应附 `certificate` 字段（见 W7）。

## 工作流划分

### W0 精确验证层（主干，最先做）

新增 `exact_check/` 包：harness + 每族一个模块。统一接口

```python
def check(kind, power, comp, bound, params) -> bool
```

重建 (m,n) 的**真** basis 矩（复用核内已有生成器，绝不复用 trig_pi.site_basis 等有偏路径）→ 从 params 还原 P 系数（`au_val/u_val` 等）→ `combine` 与核的目标构造器（如 `lemniscate_target`）比字典相等 → 再查 `poly_nonneg`。等于把 verify.py 的数值审计搬进 ℚ——更强：字典相等是全称证明，不是 50dps 采样。

同时它也是 mode=exact 的行内断言（emit 前自检）与 benchmark 判官的裁决器。特殊结构各自处理：gamma 复合积分（主核 + ln 子证明两段矩相加）、beta 兜底模板（t·poly·√(1−x⁴)[/π] + b 的矩已知）、tan/cot/tanh/coth 的整体除因子。

### W1 失真面修复（exact 路径）

- **trig_pi (1,8)**：exact 模式直接用 `basis_moment`（真矩，核内已存在），不走 `site_basis`。顺带把"预存闭式"升级为真矩递推，解除 (m,n) 存储表范围限制。
- **beta '<' 转置**：删 `transposed_lt_proof`，正确解继续向上搜（LT_M_LIMIT 是站端任意上限）；保留两个 sqrt 兜底模板（其恒等式本身为真，W0 会核证）。落地决定：exact '<' 扫描上限 EXACT_LT_LIMIT=256——更深的真证明（如 gauss 1<1398/1675 需 m~10⁶）诚实 NoSolution，预算上限作为可调参数留档；bench/verify.py 的 gamma 重建伪影不单独修，exact_check 已取代它做正确性裁决（verify.py 保留为站端审计工具）。
- **系数缩放**：站端把命题规约到 `C ⋚ bound/power` 求解却印 `power·C − bound` 与未缩放被积函数（claimed = power×actual）。exact 路径令目标向量 = 印刷 LHS 向量：解 `{C: s·power, 1: −s·bound}` 或解 `C ⋚ bound/power` 后被积函数整体乘 |power|——两种印法任选，关键是 ∫印刷物 == 印刷 LHS。
- **ln_q_square q∈{5,7}**：诊断站端崩溃根因（疑似矩系统奇异/实现除零），exact 路径按真系统正常求解。
- **方向判定认证**：float64 预检/兜底替换为 `mp.iv` 区间比较（精度递增至符号确定）；代数常数（golden=√5 型）走 ℚ 精确平方比较。等值命题报 `二者相等` 保持。

### W2 正确性 benchmark

`bench/judge_correct.py`：语料 = 全 golden 输入（ground truth 用区间比较认证，不信仰站端标签）+ 随机/对抗生成（紧界 1e-30、巨大分子分母、power∈{0,±1/2,±2,…}、边界等值点）。指标：每条 emitted proof 过 W0（**零容忍假恒等式**）、真命题证明覆盖率、方向判定正确率、与 site 模式的输出分歧清单（分歧应恰好落在失真簇上）。运行仍全本地，不打站端。

### W3 类型扩展（主体扩展面，下一阶段并行铺开）

作者文章目录 ~45 型，站端活 29。每型 = 推导矩空间 + 一个核文件 + 进注册表 + benchmark 用例，天然一型一 agent。按可达性排序：

1. **矩机器现成的**：ζ(5)、ζ(7)…（quadlog 的 ln^{k} 核直接给出 η/β 矩，kernel-spec 关键数学事实节已载）；(ln q)ⁿ（log_family 升幂同构）；π^k 更宽幂次。
2. **Beta 积分族**：arcsin q、arccos q（1/√(1−x²) 核或变量替换进 lemniscate 机器）；Γ(1/4)、Γ(1/3)（B(a,b) 表，ϖ/G 族相邻）。
3. **ln π、e^π 已在外需复核边界**；ψ′(q)（级数核）、Glaisher、Li₂(q)、Si/Cin、Ein/Ei、erf、arsinh/arcosh —— 各需一次矩空间推导（调研任务，先出推导笔记再实现）。
4. decompose 数学模式：子界分配用连分数最佳逼近保证可证；项文法扩到商/幂/嵌套。

### W4 证明器强度

Padé 插值证法（kernel-spec 已有完整规格：误差函数恒正有理函数，构造性严格）作 ln/arctan 族的独立第二证法，与 (m,n) 搜索互为补充——能证搜索预算外或高次 P 才够的界；P 升次（≥3 系数，矩空间维数随升）；预算策略改自适应（先小预算快答，未命中递增）。

### W5 证书输出与接口

exact 响应附 `certificate`：`{integrand, domain, nonneg: {rule, coeffs}, moments: 目标向量, basis: 符号表}`——一个独立脚本（或 sympy/mpmath 小段）可离线复核。CLI `attention-calculator --exact`；docs 出证书格式 spec。远期可出 Lean/sympy 可检形式（非承诺）。

## 推进顺序与分工

W0/W1/W2 本期落地：leader 先做 mode 管线 + exact_check 骨架 + solve.py 公共修复（power 归一、区间方向）作收敛点；然后 5 个 agent 按族文件归属并行（quadlog / exp+hyperbolic / log / trig_q+trig_pi / beta+gamma），bench agent 独立文件并行。全部判官保持绿是硬约束——site 路径任何 agent 不许动。W3–W5 待 W2 数字出来后按上面优先级展开，每型一个 agent。
