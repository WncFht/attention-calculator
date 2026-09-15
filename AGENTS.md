# 注意力计算器复现项目

目标：完整复现 zhuyidao.net（注意力计算器）——输入"常数 ⋚ 有理数"，输出积分恒等式证明。

## 必读文档

- `docs/kernel-spec.md` — 全部 29 种类型的核函数/矩空间/搜索顺序（唯一事实源，改实现先改它）
- `docs/api-spec.md` — 线上 API 协议（实测）
- `bench/README.md` — benchmark 数据集与评测口径

## 代码风格（继承自用户全局约定）

- 精确算术用 `fractions.Fraction`，禁止 float 进入证明路径；数值只用于校验。
- 不防御性编程：不为不可达场景加 guard/try-except；函数都要 docstring。
- minimal code，线性控制流优先，不为一次性代码建抽象；十行相似代码好过过早抽象。
- 项目自有标识符不用前导下划线；不发明 dunder。
- 矩表等预计算产物是数据文件（JSON/pickle），生成脚本放 `tools/`，产物放 `src/attention_calculator/tables/` 或 `bench/data/`。

## 工作流约定

- 重计算（矩表生成、大规模验证、benchmark 抓取）在远程机 `devbox` 上跑（ssh devbox，12 核，Python 3.14）；本机只做开发。
- 远端工作目录 `~/attention-calculator`；用 rsync 同步仓库，不要手工改远端文件。
- 所有结论以 benchmark 为准：golden 数据集在 `bench/data/golden.jsonl`，评测脚本 `bench/`。
