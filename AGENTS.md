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

## 格式化工具链

`git commit` 会走 pre-commit：formatter 经 git-format-staged 改写暂存内容并同步工作区（两侧一致，commit 不被格式化阻断，未暂存编辑不受污染）；check 类 hook 失败才拦。前置条件：`npm install`、系统工具（本机 Arch 走 `pacman -S shfmt shellcheck gitleaks taplo-cli ruff actionlint`，CI ubuntu-latest 走 brew 同源）、`pre-commit install`。npm 工具版本以 `package.json` 为唯一事实源，编辑器（`.vscode/settings.json` 的 prettierPath）与 hook 同源。

- `*.sh`：`shfmt -i 2`（gfs）+ `shellcheck -S warning`。zsh 脚本不在链内——两者都不支持 zsh，`scripts/fmt-shell.sh` 对 zsh shebang 原样透传。
- `*.py`：`ruff format`（gfs）+ `ruff check`（规则集在 `pyproject.toml`）。
- `*.md`：`markdownlint-cli2 --fix` 原地改写（改写会 fail 一次，重新 `git add` 再提交）→ `autocorrect --stdin | prettier`（gfs）。
- `*.yaml`/`*.yml`：`prettier`（gfs）。缩进规则：yaml 2 空格、md 4 空格，见 `.prettierrc` overrides。
- `*.toml`：`taplo fmt`（gfs）；`uv.lock` 是生成物，hook 已排除。
- `.github/workflows/*`：`actionlint` check。
- `src/attention_calculator/templates/*.html`、`bench/probes/*.json` 是 parity 夹具/数据——不在链内，formatter 不许碰（已在 `.prettierignore` 划出）。
- `docs/` 整目录划出 markdownlint/autocorrect/prettier——inline code span 逐字引用站端字符串（全角标点、带空白回显），内容即数据，格式化会篡改证据（hook `exclude:` + `.prettierignore` + cli2 ignores + 编辑器 `markdownlint.ignore` 四层划出）。
- gitleaks 拦 secret；自定义规则与精确值放行写法见 `.gitleaks.toml` 注释。

CI（`.github/workflows/ci.yml`）与本地同源，本地不过 CI 必挂。

## Agent 入口约定

`AGENTS.md` 与 `.agents/` 是入库的唯一事实源（skill 本体在 `.agents/skills/`，各带 `agents/openai.yaml` 元数据）。`CLAUDE.md` 与 `.claude/skills/` 是指向它们的软链，属本机便利层、不入库；clone 后跑 `scripts/agent-links.sh` 重建。
