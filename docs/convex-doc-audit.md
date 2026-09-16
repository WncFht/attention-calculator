# convex-behavior.md 独立复核报告

日期：2026-09-16（时间点文档，一次性复核快照）。

复核对象：`docs/convex-behavior.md`（303 行）vs `bench/data/convex-probes.jsonl`（333 条）。 方法：全量 digest + 逐节断言→记录对照 + 程序化交叉扫描。全程离线，未探站。

## 结论

**零硬矛盾。** 文档每条可检验断言都找到至少一条探针记录佐证，且全部一致。 以下为逐节证据定位 + 7 处"未被数据钉死但亦无矛盾"的弱支撑点 + 2 处措辞瑕疵。

注：本复核完成后 convex2 已移植上游 `prove.py` 源码实现，§二所列 7 处弱支撑常数与 站端同源（同一源码出处），其取值不再需要探针钉死。

## 一、关键断言 ↔ 记录定位（逐节）

### §1 传输层
- POST-only / 500s：`http:get` `http:put` `http:trailing`(/convex/prove/) `http:case`(/convex/Prove) `http:root-post`(/convex/) `http:get-params` — 均 500 `服务器内部错误，请稍后再试。`（带句号、带 ok:false）。
- OPTIONS 200 空体 `http:options`；HEAD 500 空体 `http:head`；GET /convex/en 200 `http:get-en`。
- JSON/text-plain 不解析：`http:json` `http:textplain` → 400 请输入一个不等式。
- query string 不参与：`http:query-post` → 400。
- 多余字段忽略：`http:extra`(foo=bar) `http:lang` `http:no_line` `param:no_line` `param:no-line` `param:noline`。
- multipart + urlencoded 均接受：`http:multipart` `http:urlenc-line` 均 200。
- 包络字节级核对（raw 字段）：sort_keys、紧凑分隔符、`\uXXXX`、尾部 `\n` 全部实测吻合（atom:xx / err:empty raw）。
- `result.ok ≡ status=="proved"`：225/225 程序化验证，0 违反。
- 500 字符闸门：`len:cap500`(500 通过→solver 错) `len:cap501b`/`len:501`(501→输入过长) `len:raw501`(raw 501/strip 500→通过) `len:blank500/501`(全空白→空错)。
- `line`/`domain` 无长度限制：`line:long`(3001 字符→进入解析→maximum recursion depth exceeded)、`dom:long`(3000 字符→进入 domain 校验)。
- 内部异常一律 400 str(exc)：`err:xd0`(float division by zero)、`err:infexp`(cannot convert Infinity)、`line:long`(recursion)。合法请求零 500 ✓。

### §2 输入语法
- 改写三规则：`atom:ln`(ln→log) `atom:estar`(e**x) `atom:ecaret`(e^x→exp(x)，证明 ^→** 在 e**x 规则之前) `err:eparen`(e^(x) 不匹配→BinOp unsupported)。
- 大小写：`err:E`(E^x) `err:EXP` `err:Log` → unsupported atom。
- 缺括号：`err:ln_noparen` `err:exp_noparen` → invalid syntax。
- 不等号序与切分：`err:eqeq`/`err:neq`/`err:xgtgt`(==/!=/>> →invalid syntax) `err:chained`(x>0>1→unsupported atom Compare)。
- 无不等号→右 0 按 >=：`dir:noop`(裸 x→proved)。
- `=` 按非严格：`dir:eqconvex`(x^2-2x+1=0→proved+proof，formula `x^{2} + 1 \ge 2x \ge 2x` 首段 \ge ✓) `dir:eq`/`dir:eqx`。
- `<`/`<=` 翻转保持/取消严格性：`dir:lt`/`dir:le`(x^2+1<2x→翻转后 right=x^2+1 convex→inconclusive ✓) `dir:flip2`(0<x→x>0→failed strict)。
- Call 检查序（个数→裸 x→名字分派）：`err:max`(max(x,1)→expects one argument，证明个数先于名字) `err:log2arg`/`err:log0arg` `err:log2`/`err:exp2x`/`err:loglog`/`err:sqrtx2`/`err:expX`/`err:expneg`/`err:expx1`/`err:expexp`(→only supports argument x) `err:foo`/`err:sin`/`err:tan`/`err:abs`(→unsupported atom Call)。
- 常量形态：`atom:hex`(0x10→16) `atom:underscore`(1_000) `atom:1e3` `atom:sqrt4`(sqrt(4)→2) `atom:int`/`atom:frac`/`atom:dec` `atom:starstar`(x**2 原生 ** 合法)。
- `expected numeric constant`：`err:xx_pow`(x^x) `err:1dx`(1/x) `err:xdx`(x/x)。
- `products of two non-constant atoms`：`err:xmulx`(x*x) `err:explogmul`。
- 其它 unsupported atom：`err:y`/`err:yalone`/`err:nan`(Name) `err:tuple`/`err:list`/`err:dict` `err:2dx`(2^x) `err:negpow`((-x)^2) `err:parenpow` `err:exppow` `err:xgrp`(x*(x+1))。
- 系数分配：`atom:grp3`((x+1)*3→3*x+3) `norm:biggroup`(2*(x+exp(x))+3→2*x + 2*e^x + 3) `atom:xneg1d2`(x^-1/2→½x^-1，证明除法折入系数) `atom:x2r`(x*2→2*x)。
- 优先级 `x^-1/2` = x**(-1)/2：`atom:xneg1d2` left `\frac{1}{2}x^{-1}` ✓。

### §3 normalized
- 首现序保留：`norm:order4`(x+x^2→x + x^2) vs `norm:order5`(x^2+x→x^2 + x)；`norm:order1/2/3`。
- 按号拆分不按凹凸：`curv:negsq`/`norm:negsq`(-x^2+10→left `10` right `x^2`，文档原例) `curv:x2sqrt`(x^2+√x 正系数凹项留左→mixed)。
- |c|≤1e-9 丢弃含等号：`curv:coeff-eps`(1e-9·exp 消失→left `1`) `curv:coeff-eps2`(1.1e-9 保留显示 `0`，曲率仍 convex、f_min=1.1e-9 内部生效 ✓)。
- 合并同原子：`norm:dup`(exp+exp+x+x→2*e^x + 2*x) `norm:cancel`(x^2-x^2 消去) `atom:xx`(x+x→2*x)。
- 渲染：`atom:longcoef`(0.123456789→10/81) `disp:1e-6coef`(1e-6→1/1000000*x) `atom:01coef`(0.1→1/10) `norm:coeffrac`(1/3*x) `atom:x15dec`/`atom:x32`(x^1.5→x^{\frac{3}{2}}) `atom:xsqrt2`(→665857/470832) `atom:x0`(x^0 不化简) `atom:x24`(x^(2/4)→√x) `atom:xhalfdec`/`atom:xhalf`(x^0.5/x^(1/2)/sqrt(x) 全→√x) `atom:xneg12`(x^{\frac{-1}{2}} 负号进分子) `atom:x20`(x^2.0→x^2) `atom:x1`(x^1→x^1 不化简) `err:bigint`(x^1e+20)。
- latex 系数并置：`atom:coefexp`(2e^x) `tan:*` 各记录（`\frac{1}{10}x` 等）。

### §4 曲率
- EPS 乘积判零边界对：`curv:eps2`(x^1.0000000001→affine) `curv:eps1`(x^1.000000001→convex)。
- `mixed` 第四值：`curv:x2sqrt`/`curv:x2log`/`curv:explog`/`curv:x2x3`/`curv:sqrtlog`。
- 模板规则：程序化扫描 225 条 result——inconclusive ⟺ left∉{convex,affine} ∨ right∉{concave,affine}，零违反。模板通过时 difference 恒 convex/affine（第二道检查不可达 ✓）。
- 幂指数曲率符号：`curv:099`(concave)/`curv:101`(convex)/`curv:neg2`/`curv:neg01`/`curv:01`/`curv:third`。
- 多凹原子合成：`curv:twoconc`(right `√x + x^0.6` concave→模板过→failed 于数值)。

### §5 最小值
- lo_eff=max(lo,1e-8)：`min:edge`/`min:edge2`(x>0 系最小值落 1e-08) `dir:noop` `atom:x0`。
- 上界分支：`dom:bind`(100>x dom=0,10→d(hi)≤0→(10,90)) `dom:lobound`(dom=5,10→d(lo)≥0→(5,20)) `dom:hi-xstar`(dom=0.5,inf→x*=0.5 边界) `dom:inverted`(10,0→(10,90)，文档原例)。
- 内部二分：`min:interior`(x*=2) `min:golden`(x*=0.5,f=-1.25) `min:quarter`(0.25,-0.25) `dir:explog261`(0.56714329041)。
- b 倍增至 2^27=134217728（非 1e8 截断）：`min:1e8`(x^-1→x_text `134217728`,f=7.45e-9) `min:decr`/`min:decrconv`/`atom:hex`/`atom:underscore`/`curv:neg2`/`curv:neg01`/`curv:twoconc`。
- `%.12g` 文本 + 原生 float repr：全部 min 记录吻合（含 `-1.11022302463e-16`、`9.99999981579e-09` 等）。
- 下溢：`min:underflow`(x^1000→f=0)。

### §6 判定容差
- 容差表 8 行逐行核对：`tol:eq0gt`(0,failed) `tol:eq0ge`(0,proved) `tol:pos5e9`(4.9999999696e-9,failed) `tol:pos1e8`(9.9999999392e-9,failed) `tol:pos2e8`(1.9999999878e-8,proved) `tol:neg1e9ge`(-1.0000000827e-9,proved) `tol:neg1e8ge`(-9.9999999392e-9,proved) `tol:neg11e8ge`(-1.1000000021e-8,failed) + `tol:neg1e9`/`tol:pos1e8ge`。
- status 映射：`ok`→proved 恒成立（程序化）。

### §7 切线候选搜索
- 触发条件（left 严格 convex 才搜）：`norm:negsqrt2`/`tan:left-affine`/`norm:neglog`(left affine→proved 无 proof) vs `min:interior`/`dir:explog261`(convex→有 proof)。
- CF 渐近分数 + 首中序悬崖：`tan:explog:c=1.8`→1、`c=2.1`→1/2、`c=2.31`→4/7、`tan:cliff1`(2.3303)→17/30、`tan:cliff2`(2.33036)→17/30、`tan:cliff3`(2.330366)→38/67、`tan:cliff4`(2.330366124762=c* 零最小值)→failed。
- 分母 ≤10000 含等号：`cf:denom1e4`(1/10000 接受) `cf:denom-over`(1/10002 被滤→无证书) → 上界钉在 {10000,10001}。
- `x0 ≤ lo_raw` / `x0 ≥ hi_raw` 跳过：`tan:zero-gap`(候选仅 0/1，0≤0 跳过→集空) `cf:denom-over` `dom:neglo-ta0`(lo_raw=-1→x0=0 不跳→tan='0') `dom:neglo-crash`(x0=0 处 √x 求切→0.0 cannot be raised) `dom:hibound`(x*=0.3968,hi=0.4→1/2=0.5 被跳→落到 1/3)。
- `>=`/`>` 不同切点：`tan:ge-variant`(>=0→tan=1，gap -1.1e-16 被接受) vs `param:no_line`/`tan:explog:c=2`(>0→x0=1 拒→1/2)。
- 搜索无果 reason：`tan:zero-gap`/`min:explinx`/`tan:noright`/`cf:denom-over` 等 33 条同一原文。
- proof 字段与模板：`dir:explog261` line_text=`(30/17*x - 1 + ln(17/30)) + 261/112`（文档引用样例逐字）；`tan:pow23`/`tan:pow13` 幂模板（latex `^{\frac{2}{3}-1}` 拆指数 ✓）；`tan:multiright` 多项 ` + ` 连接；`tan:ge-variant`/`tan:plain` lltex `(1x - 1 + \ln 1)` 不化简；`tan:ge-affine`/`tan:ge-tangent2`/`tan:ge-affine2`/`tan:gt-affine` 仿射 right→h=right。
- 空 right→空 line→双空格：`atom:xneg12` formula `x^{\frac{-1}{2}} >  \ge 0` 逐字；`curv:neg01` tan=134217728。
- 首段 rel 只反映严格性：`dir:eqconvex`(\ge) `tan:ge-*`(\ge) `tan:gt-affine`(>) `tan:plain`(>)。

### §8 line 参数
- 门控（仅 proved 且非空才解析）：`line:empty`(空→pl=null) `line:onfail`/`line:oninconcl`(不解析) `err:precedence`(主式 foo(x) 错→line=x^2 未读)。
- 仅仿射：`line:nonaffine`(x^2) `line:log` `line:exp` → `--line must be an affine expression in x` 逐字。
- 同套语法：`line:ln`/`line:out-syntax`(ln(2)→log only supports argument x) `line:badparse`(foo(→'(' was never closed) `line:badatom`(foo(x)→unsupported atom) `line:group`(2*(x+1) 分配) `line:article`(41/14*(x-12/161) 接受) `line:frac`/`line:space`/`line:xx`/`line:neg`/`line:float`/`line:const`/`line:zero`。
- 双 gap ≥ -1e-8：`line:x`/`line:const`/`line:frac`/`line:zero`/`line:xx`(ok:true) vs `line:group`/`line:neg`/`line:space`/`line:bad`/`line:simple`/`line:float`(ok:false)，与各自 gap 值全部一致。
- 纯诊断不回灌：`line:x`(ok:true 但 proof 仍 null)。

### §9 domain
- 段数校验：`http:domain-bad`(abc→domain must look like) `dom:halfbad`('0,'→could not convert '') `dom:badlo`('a,b'→could not convert 'a')。
- lo_raw vs lo_eff 双界：`dom:neglo`/`dom:neglo-ta0`/`dom:neglo-crash`（§7 已列）。

### §10 reason 全表
- 数据全集恰 5 串（程序化去重）：空串×77、`没有生成中间直线证明`×25、`内置有限候选搜索`×33、`模板不符`×25、`数值最小值未达到`×65 —— 与文档五串逐字一致。

### §11 错误目录
- 文档列出的每一类错误串在数据中均有 ≥1 条记录（含 unicode_ge/unicode_minus/times 三字符、`invalid syntax` line 0/1 两变体、`(' was never closed`、bigint、sqrtneg、pipe、semi、import、dot、blank 等）；反向扫描数据中全部错误串均落在文档目录内，无漏收。

### §12 怪癖
- 全部有记录支撑（上文已逐条定位）：`=` 判 proved、假命题宽容、x^0 不化简、x^(2/4)→√x、1e-9 消失/1.1e-9 显示 0、边界最小值无切点、line_latex 不约分、proved 可无 proof、mixed 第四值、2^27、异常全部 400。

## 二、弱支撑点（无矛盾，但数据未完全钉死——值来自参考实现）

1. **非严格候选接受阈值 `≥ -1e-7`**（§7）：唯一区分性证据是 `tan:ge-variant` x0=1 gap≈-1.1e-16 被接受——与 `≥-1e-7` 和 `≥-1e-8` 都兼容；缺 gap 落在 (-1e-7,-1e-8) 区间的探针。
2. **非严格整体判定 "含等号" `≥ -1e-8`**（§6）：最近探针 `tol:neg1e8ge` f_min=-9.99999993922529e-9 严格高于 -1e-8，`>` 与 `≥` 不可区分。
3. **`provided_line` 双 gap `≥ -1e-8`**（§8）：ok:true 的最小 gap 为 0（line:zero）与 +1e-8（line:x），ok:false 的最大 gap 为 -0.25——(-1e-7,0) 判定带无探针。
4. **CF `max_terms=12`**（§7）：数据只证明 ≥5 项（38/67 是第 5 渐近分数仍被生成）；12 的具体值无探针区分。
5. **bisection "160 次折半"**（§5）：float64 在 ~53 次后饱和，迭代次数从输出不可区分；数值结果一致。
6. **`domain` hi 集合 `{inf, infinity, +inf}`**（§9）：仅 `inf` 实测（dom:neglo* 等）；`infinity`/`+inf`/大写 `INF` 无探针。`hi=''` → float 错（dom:halfbad）间接支持"不在集合内走 float()"分支。
7. **常量 `a**b`（如 `2**3`）**（§2）：无直接探针（`atom:starstar` 是非常量侧 x**2）；分配 `(x+1)/2` 亦无探针（仅 `atom:xneg1d2` 证明除法折入系数）。

另注：候选严格侧 `>1e-8`（§7）——拒绝事件存在（ge-variant 的 `>` 版本 x0=1 被拒）但被拒 gap 值未打印，精确界未隔离；与整体判定的 1e-8 一致性合理。

## 三、措辞瑕疵（非错误）

- §10 "227 例 200 响应的全集"：227 是 HTTP 200 计数（含 OPTIONS 空响应 + GET /convex/en HTML），实际携带 result/reason 的为 225 条。五串确为 result-bearing 全集，措辞略松。
- §1 "≤1.05 req/s"：实测 333 条平均 ~0.125 req/s（ts 跨度 2663s），不等式成立但数值来源不明。

## 四、总评

文档可放心作为实现 spec 使用：所有转录（reason 原文、line_text/formula 模板、normalized 项序、minimum 数值、错误串目录、容差阶梯、CF 悬崖序、包络字节格式）均与 333 条语料逐字/逐值吻合。§二列出的 7 处弱支撑点全部朝"参考实现给定常数"方向，数据无一反例；如需把它们也钉死，需追加探针（gap 落在 -1e-7~-1e-8 的 `>=` 用例、`domain=0,infinity`、`2**3>0`、CF 第 6-12 项用例等）。
