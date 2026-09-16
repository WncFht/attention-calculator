# 证明证书格式（W5，2026-09-16）

mode=exact 每条 emitted proof 附一份机器可检证书：`src/attention_calculator/certificate.py` 的 `build()` 产出、`verify_cert()` 独立复核、`tools/verify_cert.py` 做离线 CLI。证书的目标是"持有者可不信我方管线、只信这段短代码即可复核命题"——重算走 `exact_check.verify`（各族 checker 从核内真矩生成器重建 basis，与被积函数矩做 ℚ 字典相等比较），证书本身只是断言记录。

## 字段

```json
{
  "kind": "pi",
  "power": "1",
  "comparison": ">",
  "bound": "3",
  "parameters": {"m": 0, "n": 4, "a_val": "0", "b_val": "0", "c_val": "0", "au_val": "47", "bu_val": "-13", "cu_val": "0", "u_val": "120", "unified_form": {}},
  "check": {
    "integrand": {"1": "-3", "pi": "1"},
    "target": {"1": "-3", "pi": "1"},
    "identity_ok": true,
    "nonneg": true
  }
}
```

- `kind`：命题类型，`solve.FAMILY` 注册名（site 29 型 + `EXACT_TYPES`）。
- `power` / `bound`：请求的常数系数与有理界，`str(Fraction)` 序列化——分母为 1 时是 `"n"`，否则 `"n/d"`；不带空白、不接受 float。
- `comparison`：`">"` 或 `"<"`，原命题方向（`power·C ⋚ bound`）。
- `parameters`：/calculate 响应的 parameters 原样透传（`m`、`n` 为 int，`*_val` 为有理数字符串，`unified_form` 为对象；`ln_pow` 等核可扩展 `du_val` 等额外键）。证书的数学内容全靠它：checker 从中还原 P 系数与 (m,n)。
- `check`：`exact_check.verify` 的结果记录。`integrand` / `target` 是 Moment 序列化——`{常数符号: "n/d"}`，符号排序、缺省即零；`integrand` 是被积函数矩的精确值，`target` 是印刷声称向量 `s·(power·C − bound)`（`s` 由 comparison 定）。`identity_ok` 记 `integrand == target` 的 ℚ 字典相等判定；`nonneg` 记 P 的定号规则（`poly_nonneg`、beta 模板的 `t,b ≥ 0`、gamma 的结构性引理等）加上被积函数非恒零。

## 验证语义

`verify_cert(cert)` 为真当且仅当三件事同时成立：证书通过 schema 校验（字段齐全、类型正确、有理数可解析、`kind` 已注册）；用证书内 `parameters` 重跑 `exact_check.verify` 得到的 `identity_ok`、`nonneg` 皆为真；且重算结果与记录的 `check` 块逐字一致。第三条是关键——记录的矩与判定不是证据而是声称，篡改 `bound`、系数、记录矩或判定位都会造成声称与重算不符。重算不依赖证书里任何数学断言，basis 矩完全从核内生成器重建，且 `verify` 先做参数域校验（`m,n ≥ 0`、`u_val > 0`、`power` 在 kind 域内）并把可由 kind/power 派生的字段（`c_val` 等）重算比对——不是只查格式。

因此证书证明的是：`integrand == target` 是 ℚ 上的全称恒等式（不是数值采样），且 `nonneg` 给出被积函数在积分域上定号的精确证书；两者合起来蕴含 `power·C ⋚ bound` 为真。

## 序列化与 CLI

规范化序列化用 `json.dumps(cert, sort_keys=True, ensure_ascii=False)`；Moment 内部键已排序。`cert_tex(cert)` 给一行人类可读摘要（target 向量渲染为 `22/7 - \pi = \int f dx > 0` 形）。离线复核：`tools/verify_cert.py cert.json`（或 stdin），打印重算的 `identity_ok`/`nonneg`/命题行，退出码 0/1/2 分别对应核验通过、核验失败、输入或格式错误。

## Padé 证书变体（W4 第二证明器）

`prove_exact` 在 (m,n) 搜索 NoSolution 后回落 `pade.py`（ln_q/arctan_q 的 Padé 插值证法）；其响应无 `parameters`，`prover: "pade"`，证书是另一套 schema：

```json
{"kind": "ln_q", "comp": ">", "q": "...", "p": "...",
 "n": 1, "m": 3, "a": "...", "b": "...", "resid": "0", "serr": "..."}
```

含 `serr` 键即 Padé 证书，`verify_cert`/`recheck`/`verify_response` 见到该键自动分发到 `pade.verify_cert`：重推 [L/M] 逼近与误差单项形（正系数单项 + 分母无正极点）、重核 `resid` 恒等式，全程 ℚ。配对规则与退化（resid≠0）情形见 `docs/2026-09-16-pade-notes.md`。
