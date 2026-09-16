# mode=exact decompose: provable-by-construction bound allocation

The site allocator (`decompose.py`) splits a composite inequality `Σ c_i·Π atoms ⋚ R` into per-term rational bounds chosen for byte parity with zhuyidao.net — record-chain-style approximants with no provability argument. This note documents the exact-mode allocator in `src/attention_calculator/decompose_exact.py`, which instead allocates bounds so that every emitted sub-claim is *provable by construction* inside the (m,n) search budget, then verifies that claim by actually running `solve.prove(..., exact=True)` on each step.

Related docs: `docs/decompose-notes.md` (frozen site behavior spec) and `docs/2026-09-16-decompose-exact-kinds.md` (atom kind coverage).

## The provability frontier

Each kernel proves `C ⋚ β` by scanning `(m,n)` in depth-major order (`engine.mn_order`); a bound is reachable iff the moment system has a non-negative solution at some `(m,n)` within the per-type limit. Empirically (and verified exactly for quadlog by solving the affine-in-`r` system at `r∈{0,1}`), the reachable set per `(kind,power,comp)` is a **one-sided interval** `(-∞, r_max]` — or `[r_min, ∞)` for `'>'` — with no gaps below the frontier. So a bound proves iff its margin `|U − β|` clears a per-unit **frontier floor** `ε(kind,power,comp)`.

`bench/decompose_math_probe.py` maps the frontier: for each unit reachable from the decompose grammar it scans a decimal ladder of margins `10^-k` and records the minimal depth `m+n` of the first proof (`NoSolution` marks the rung below the budget frontier).

| kind | power | comp | 1e-0 | 1e-1 | 1e-2 | 1e-3 | 1e-4 | 1e-5 | 1e-6 | 1e-7 | 1e-8 | 1e-9 | 1e-10 | 1e-11 | 1e-12 | 1e-13 | 1e-14 | 1e-15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pi | 1 | > | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 |
| pi | 1 | < | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 |
| pi | 8 | > | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 | 32 |
| pi | 8 | < | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 | 30 | 32 |
| e | 1 | > | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 8 | 9 | 10 | 11 | 11 | 12 | 13 | 13 |
| e | 1 | < | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 8 | 9 | 10 | 11 | 11 | 12 | 13 | 13 |
| e_q | 2 | > | 3 | 4 | 5 | 7 | 8 | 9 | 10 | 11 | 12 | 12 | 13 | 14 | 15 | 16 | 17 | 17 |
| e_q | 2 | < | 3 | 4 | 6 | 7 | 8 | 9 | 10 | 11 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 17 |
| e_q | 1/2 | > | 0 | 1 | 2 | 3 | 3 | 4 | 5 | 6 | 6 | 7 | 8 | 8 | 9 | 9 | 10 | 11 |
| e_q | 1/2 | < | 0 | 1 | 2 | 3 | 3 | 4 | 5 | 6 | 6 | 7 | 8 | 8 | 9 | 9 | 10 | 11 |
| e_pi | 1 | > | 5 | 8 | 10 | 12 | 14 | 16 | 18 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| e_pi | 1 | < | 6 | 8 | 10 | 12 | 14 | 16 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| pi_n | 2 | > | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 20 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| pi_n | 2 | < | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| pi_n | 3 | > | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| pi_n | 3 | < | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| pi_n | 5/2 | > | 1 | 1 | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | NoSolution | NoSolution | NoSolution |
| pi_n | 5/2 | < | 0 | 0 | 2 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | 20 | NoSolution | NoSolution | NoSolution |
| ln_q | 2 | > | ValueError | 1 | 2 | 4 | 5 | 6 | 8 | 9 | 10 | 11 | 13 | 14 | 15 | 17 | 18 | 19 |
| ln_q | 2 | < | 0 | 1 | 2 | 4 | 5 | 6 | 7 | 9 | 10 | 11 | 13 | 14 | 15 | 17 | 18 | 19 |
| ln_q | 3 | > | 0 | 2 | 4 | 5 | 7 | 9 | 10 | 12 | 14 | 16 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution |
| ln_q | 3 | < | 0 | 2 | 4 | 5 | 7 | 9 | 10 | 12 | 14 | 16 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution |
| ln_q | 3/2 | > | ValueError | 0 | 2 | 2 | 4 | 4 | 6 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
| ln_q | 3/2 | < | 0 | 0 | 2 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
| sin_q | 1 | > | ValueError | 1 | 1 | 3 | 4 | 5 | 6 | 7 | 7 | 8 | 9 | 10 | 11 | 11 | 12 | 13 |
| sin_q | 1 | < | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 7 | 8 | 9 | 10 | 11 | 11 | 12 | 13 |
| cos_q | 1 | > | ValueError | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 8 | 9 | 10 | 10 | 11 | 12 | 13 |
| cos_q | 1 | < | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 8 | 9 | 10 | 10 | 11 | 12 | 13 | 13 |
| tan_q | 1 | > | 0 | 1 | 3 | 4 | 5 | 5 | 7 | 7 | 8 | 9 | 9 | 10 | 11 | 12 | 12 | 13 |
| tan_q | 1 | < | 0 | 2 | 3 | 4 | 5 | 5 | 6 | 7 | 7 | 8 | 9 | 9 | 10 | 11 | 12 | 13 |
| arctan_q | 1 | > | 0 | 2 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 20 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| arctan_q | 1 | < | 0 | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| arctan_q | 2 | > | 0 | 4 | 9 | 14 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| arctan_q | 2 | < | 0 | 4 | 9 | 15 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| sinh_q | 1 | > | 0 | 1 | 1 | 3 | 4 | 5 | 6 | 7 | 7 | 9 | 9 | 10 | 11 | 11 | 12 | 13 |
| sinh_q | 1 | < | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 7 | 9 | 9 | 10 | 11 | 11 | 12 | 13 |
| tanh_q | 1 | > | 0 | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 7 | 8 | 9 | 10 | 11 | 11 | 12 | 13 |
| tanh_q | 1 | < | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 7 | 8 | 9 | 10 | 11 | 11 | 12 | 13 | 13 |
| gamma | 1 | > | 0 | 0 | 6 | 15 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| gamma | 1 | < | 0 | 2 | 9 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| golden | 1 | > | 0 | 0 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
| golden | 1 | < | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
| catalan | 1 | > | 0 | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| catalan | 1 | < | 0 | 0 | 3 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| zeta3 | 1 | > | 0 | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| zeta3 | 1 | < | 0 | 1 | 3 | 5 | 7 | 9 | 11 | 13 | 15 | 17 | 19 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| varpi | 1 | > | 0 | 2 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| varpi | 1 | < | 0 | 3 | 32 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| gauss | 1 | > | 0 | 0 | 10 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |
| gauss | 1 | < | 0 | 1 | 10 | 104 | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution | NoSolution |

(`ValueError` at 1e-0 is the log-family domain check rejecting bound < 0 — expected, not a frontier datum.)

### Depth-vs-margin relation

For every unit the minimal depth grows **linearly in digits of margin**: `depth ≈ a + b·log10(1/ε)`, i.e. the frontier margin decays geometrically `ε ~ A·ρ^d`. Rates observed: pi/pi_n/quadlog units ~0.5 decades per depth level; e/e_q/trig/hyperbolic ~1; ln_q ~0.5–0.75; catalan/zeta3 ~0.55 until their budget (limit 10 → depth ≤ 20) floors near 1e-9–1e-10; gamma ~1e-3–1e-4; varpi '>' only reaches ~1e-1; gauss '>' ~1e-2 while gauss '<' reaches 1e-3 only at depth 104. `e_pi` saturates ~1e-7 (limit-reached); `arctan_q(2)` stalls ~1e-3.

Two provisos on reading this table. First, `ln_q`/`arctan_q` have a second prover: after (m,n) exhaustion, exact mode falls back to the interpolation-Padé certificates (W4), which reach *any* margin — their effective floor is 0 (confirmed by `margin_floor` measurement: `ln_q(2)` and `arctan_q(1)` prove below 1e-20·|U|). Second, several "floors" are budget artifacts, not kernel limits: pi's row still climbs at 1e-15 (depth 30 of 60 available), so its true frontier is far below 1e-15.

Measured floors used by the allocator on the golden problems (log-bisection, ~8 oracle calls each): `gamma '<'` 1.8e-3, `varpi '<'` 2.0e-3, `gauss '<'` 6.3e-4, `pi_n(2) '>'` 9.9e-10, `golden '<'` 2.9e-19; `pi`, `e`, `sin_q`, `ln_q` measure 0 (provable below the 1e-20 probe floor within budget).

## Is CF-tightness desirable? No

A continued-fraction convergent `p/q` of `U` achieves margin `~1/q²` — the *smallest* margin attainable at denominator `q`, hence the *deepest* proof at that denominator. Tightness is pure cost: the site fixture proves `π² > 227/23` (margin 4e-5, depth ~20) where the exact allocator's `π² > 986731/100000` (margin 2.3e-3) needs only depth 6. The right objective is the opposite: margins as *large* as the slack allows while respecting each unit's floor — which also yields small denominators (snap quantum ~`1/margin`).

## Allocation algorithm

Claim `Σ_i t_i ⋚ R` with slack `S = |Σ v_i − R|` (certified by escalating-dps evaluation with the same `2^(30−dps)` guard convention as `certified_cmp`; 99% of the certified slack enters the budget).

1. **Sub-claim model.** Each factor contributes `U_j = constant_mpf(kind_j, |arg_j|)`; the term is `outer · Π U_j^{σ_j}`, `σ_j = sign(arg_j)`. Coefficient kinds (`pi`, `e`, `gamma`, `golden`, `catalan`, `zeta3`, `varpi`, `gauss`, `e_pi`) fold `|c|` into the power slot so `U = |c|·C`; argument kinds keep the coefficient in `outer`. A reciprocal factor (`arg<0`) is proved on its base — `e^{-2}` becomes `e² ⋚ β` contributing `1/β`.
2. **Frontier floors.** For each distinct `(kind,power,comp)` measure `ε` by log-bisection on the prove oracle (provable margin monotone ⇒ bisection valid; result within a small constant factor, always overestimating). Cached process-wide.
3. **Slack split ∝ frontier need.** Term `i`'s need is `need_i = |v_i|·Σ_j ε_j/|U_j|` — for a `k`-factor product the relative margin `η` must cover `Σ ρ_j` with `ρ_j` each factor's relative margin, and `ρ_j ≥ ε_j/|U_j|` is provability (for `σ=−1` factors `dt/t ≈ dU/U`, so relative margins carry over). Splitting `S ∝ need_i` gives every factor realized margin `≈ K·ε_j` with the *same* safety ratio `K = S/Σ need` — uniform proof headroom across units. A 2%-of-mean damping floor keeps zero-floor terms off the boundary (margin 0 would snap to `U` itself → EqualClaim).
4. **Factor split.** Within a product, `ρ_j ∝ ε_j/|U_j|` blended 25% toward uniform — every factor gets a nonzero share regardless of measured floor.
5. **Bound placement.** Each factor bound snaps to a decimal grid `Q ≈ 4/margin` *inside* its margin window (realized margin ∈ [share/2, share]): lower bounds `⌈(t−μ)Q⌉/Q`, upper bounds `⌊(t+μ)Q⌋/Q`, with a one-quantum bump if rounding lands on the wrong side of `U`. All emitted quantities are `Fraction`; mpf only evaluates constants.
6. **Direction rules.** A `'>'` claim needs `d_i ≥ R`-side bounds: `need_lower_j = dir·sign(outer_i) > 0` decides whether factor `j` needs a lower bound on `t_j`; the proved direction on `U_j` is `'>'` iff `need_lower_j == (σ_j > 0)` (reciprocals flip).
7. **Certificate check + escalation.** `Σ term.bounds ⋚ R` is checked exactly over ℚ. Sub-claims returning `NoSolution` get their term's weight ×4 and re-allocate (≤5 rounds); `WrongDirection`/`ValueError` are hard failures — no margin rescues them.
8. **Direct paths.** Single-atom single-term claims prove against `R/outer` itself (the loosest valid bound — no allocation needed); `atom ⋚ atom` claims chain two proofs through a rational midpoint.

## End-to-end results

| problem | slack | steps | depths | result |
|---|---|---|---|---|
| `pi^2+8*pi>35` (golden) | 2.3e-3 | π²>986731/100000, 8π>50265437/2000000 | 6, 11 | proved, sum 35.00003>35 |
| `e*pi+phi+sin(1)<11` (golden) | 7.6e-4 | e<2.718282611, π<3.141593559, φ<1618777/10⁶, sin1<8414759/10⁷ | 7, 12, 3, 5 | proved, sum 11−8e-7 |
| `pi^2+8*pi>35002/1000` | 3e-4 | π²>4934633/500000, 8π>628318447/25000000 | 8, 13 | proved |
| `e*pi<4443/520` | 2.6e-4 | e<271899/10⁵, π<3142411/10⁶ | 4, 6 | proved (bound within 1e-5 of R) |
| `pi^2/e>36/10` | 3e-2 | π²>1959/200, e<10⁶/367571 (reciprocal) | 3, 4 | proved |
| `varpi+gauss<35/10` | 4.3e-2 | ϖ<26543/10⁴, G<8451/10⁴ | 10, 10 | proved — both units near their floors |
| `gamma+pi<19/5` | 8.1e-2 | γ<821/1250, π<157119/50000 | 2, 6 | proved |
| `e^pi>2314069/100000` | 2.6e-6 | e^π>R direct | 18 | proved |
| `e^pi>2314069263/10⁸` | 2.8e-9 | e^π>R | — | honest `NoSolution` (slack < e_pi floor ~1e-7) |
| `pi^2+8*pi<35` | — | — | — | certified false, `claim_certified_false` |

`tests/test_decompose_exact.py` pins these: 14 cases green in ~100 s (floor measurements dominate; cached process-wide).

## Remaining gaps

- **Honest frontier, not magic.** Claims whose slack falls below `Σ need` (e.g. `varpi '>'` almost anything, `gamma` with slack < ~2e-3, `e_pi` below ~1e-7) report `NoSolution` steps rather than fabricating proofs — that's the correct answer given the kernels' budgets, but it means the prover budget is the ceiling on composite claims.
- **Grammar inherited.** Products may contain only non-func atoms (site rule); `ln(1/2)`-type negative units are domain-rejected by the kernels and surface as failures.
- **Floor model.** `margin_floor` overestimates ε by a small factor (conservative — wastes slack, never under-allocates); escalation rescues residual `NoSolution`s. A per-`(m,n)` reachability certificate would replace measurement, at the cost of one affine solve per unit.
- **Not wired to the server.** `decompose_exact` is a library entry point; `server.py` keeps only the frozen site path.
