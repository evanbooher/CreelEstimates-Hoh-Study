# `bss_b_T2_variability.csv` — data dictionary

One row per **basin × fishery series × likelihood type**. Produced by
`fit_series()` in `06_variability_analysis.R`, which runs a random-effects
meta-analysis over the years of one series.

Worked values throughout are **Snohomish fall salmon, vehicle (`b[1]`)**.

Everything is modelled on the **log scale**, because `b` is strictly positive
with a `lognormal(0, σ)` prior. Values reported as `b` are exponentiated back;
values reported as variances or SDs (`tau2`, `tau`, `sd_log_b_raw`) stay in log
units, where they read as proportional change — `tau = 0.05` means roughly ±5%.

---

## Identity

### `basin`
**What it is** — Skagit, Snohomish or Stillaguamish.
**How** — joined from `bss_b_comparability.csv`, derived from the fishery name.
**Tells you** — the grouping level for the T3 variance decomposition.
**Example** — `Snohomish`

### `fishery_type`
**What it is** — the fishery series: the fishery name with the year removed, so
all years of one fishery share a label.
**How** — `fishery_type_from_name()` in `common.R`; `"Snohomish fall salmon 2023"` → `"Snohomish fall salmon"`.
**Tells you** — the unit a `b` series belongs to. One row of T2 = one of these × one likelihood type.
**Example** — `Snohomish fall salmon`

### `bias_type`
**What it is** — which count type the bias term belongs to. `vehicle` is `b[1]`, `trailer` is `b[2]`.
**How** — set in `get_bss_bias()` when the two elements of the Stan `b` vector are extracted.
**Tells you** — which likelihood you are reading. **These behave very differently and should never be quoted interchangeably.**
**Example** — `vehicle`

---

## How much data went in

### `n_years_available`
**What it is** — how many years of this series have a `b` estimate at all.
**How** — `n_distinct(year_start)` **before** the usability filter.
**Tells you** — the size of the fitted record. Compare against `n_years_used`.
**Example** — `5`

### `n_years_used`
**What it is** — how many of those years actually entered the model.
**How** — `nrow(df)` passed to `metafor::rma()`.
**Tells you** — the real sample size behind every statistic to its right. A τ² on two years is a very different claim from one on five.
**Example** — `5`

### `n_years_unusable`
**What it is** — the gap between the two.
**How** — `n_years_available − n_years_used`.
**Tells you** — that estimates were dropped, and how many. **Any value above 0 means read `bss_b_T2_excluded.csv`**, which names each dropped fishery-year and why (no posterior median, no posterior SD, non-finite `log(b)`, or zero/non-finite sampling variance).
**Example** — `0` — Snohomish is the only basin where this is 0 across both likelihood types.

### `n_informed`
**What it is** — how many of the used years carry real data rather than mostly reproducing the prior.
**How** — `prior_contraction ≥ 0.10` **and** not flagged unconverged. Prior contraction measures how far the posterior moved off the prior.
**Tells you** — whether the series is describing the river or the prior. A fishery-year with no trailer index counts still produces a `b[2]` posterior — it just returns the prior. **`n_informed` well below `n_years_used` invalidates everything to the right.**
**Example** — `5` of 5

---

## The raw spread, before any modelling

### `mean_b`
**What it is** — the plain average of the yearly estimates, on the log scale and converted back.
**How** — `exp(mean(log_b))` — a geometric mean, **unweighted**.
**Tells you** — a reference point only. Differs from `pooled_b` because that one weights precisely-measured years more heavily. A large gap between them means the years differ a lot in precision.
**Example** — `1.366` vs `pooled_b` `1.352`

### `sd_log_b_raw`
**What it is** — the standard deviation of the yearly point estimates, in log units.
**How** — `sd(log_b)`.
**Tells you** — **the naive answer to "how variable is `b`", and it is an overstatement.** It contains real year-to-year change *and* the measurement error in each year, mixed together. Its whole purpose here is to be compared with `tau`.
**Example** — `0.088` (≈8.8%) against a real `tau` of `0.053` (≈5.3%) — the naive figure overstates movement by about 65%.

---

## The separation — the core of the table

### `tau2`
**What it is** — the variance of **real** year-to-year change in `b`, after each year's own measurement error has been removed.
**How** — `f$tau2` from `metafor::rma(yi = log_b, vi = vi, method = "REML")`. The model is
`log b_i = μ + u_i + e_i`, where `u_i ~ N(0, τ²)` is real change and `e_i ~ N(0, v_i)` is measurement error. **`v_i` is not estimated — it is supplied**, as the posterior SD from that year's Stan fit. Knowing the expected wobble is what makes the split possible: only the excess is attributed to τ².
**Tells you** — how much the underlying bias term genuinely moves.
**Example** — `0.00284`

### `tau`
**What it is** — the same quantity as a standard deviation instead of a variance, so it is readable.
**How** — `sqrt(tau2)`.
**Tells you** — **the single most useful number in the row.** In log units it is approximately the proportional year-to-year movement: `0.053` ≈ ±5%.
**Example** — `0.053` for vehicle; `0.187` (≈20%) for trailer.

### `I2`
**What it is** — the percentage of the observed spread between years that is real rather than measurement noise.
**How** — `f$I2`, which is `τ² / (τ² + v̄)` as a percentage.
**Tells you** — in plain terms: *of the differences you see between years, this share are real.* Near 0 means the years look different only because each was measured imprecisely. Near 100 means the bias term genuinely moves.
**Example** — `50.6%` for vehicle — half the visible wobble is real. `82.9%` for trailer.

### `Q_p`
**What it is** — the p-value from Cochran's Q test of "is there **any** real year-to-year variation?"
**How** — `f$QEp`. Q compares observed scatter against what measurement error alone predicts.
**Tells you** — whether real movement is *demonstrated*, as opposed to merely estimated. A small `tau` with a large `Q_p` means stability that is consistent with the data but not proven — which with 5 years is the expected result, not a failure. **A small `Q_p` is the stronger statement: variation is real.**
**Example** — `0.128` for vehicle (cannot demonstrate any movement) versus `0.00017` for trailer (movement demonstrated at better than 1 in 5000).

---

## The series mean and its interval

### `pooled_b`
**What it is** — the model's estimate of the series' underlying average bias.
**How** — `exp(f$b)`. An inverse-variance weighted mean with weights `w_i = 1/(v_i + τ²)`, so a precisely-measured year counts more. τ² sits in the denominator too, which stops any single year dominating once real variation exists.
**Tells you** — the best single number for the series as a whole. **Not the right number for predicting a new year on its own** — pair it with the prediction interval.
**Example** — `1.352`

### `ci_lb` / `ci_ub`
**What it is** — the 95% confidence interval on `pooled_b`.
**How** — `exp(predict(f)$ci.lb)` / `$ci.ub`; on the log scale, `μ ± 1.96 × SE(μ)`.
**Tells you** — **where the *average* `b` of the observed years lies.** It answers a question about the past, and it shrinks toward zero width as years accumulate.
**Example** — `1.261 – 1.450`

---

## The prediction — what a year without census gets

### `pi_lb` / `pi_ub`
**What it is** — the 95% prediction interval: where a **new, unobserved** year's `b` is expected to fall.
**How** — `exp(predict(f)$pi.lb)` / `$pi.ub`; on the log scale, `μ ± 1.96 × √(τ² + SE(μ)²)`. The only difference from the CI is `τ²` under the root — a new year is its own draw off the series mean, not another look at it.
**Tells you** — **this is the answer to "what do we use for a fishery-year we cannot measure, and how uncertain is it?"** Unlike the CI it never shrinks below `±1.96 τ`, however many years you add, because τ is a property of the river rather than of the sample size.
**Example** — `1.193 – 1.533` for vehicle; `0.568 – 1.297` for trailer.

### `pi_width_ratio`
**What it is** — how many times wider the prediction interval is than the confidence interval.
**How** — `(pi.ub − pi.lb) / (ci.ub − ci.lb)`, computed on the log scale. Algebraically this is exactly `√(τ² + SE²) / SE`.
**Tells you** — **the factor by which you would understate this year's uncertainty if you quoted the CI where the PI belongs.** That substitution is the most likely misreading of this table. It also indicates what more years can buy: a ratio near 1 means the interval is already dominated by uncertainty in the mean (more years help); a large ratio means it is dominated by τ (more years barely help).
**Example** — `1.80` for vehicle; `2.16` for trailer.

---

## The sensitivity check

### `tau2_informed_only`
**What it is** — `tau2` recomputed using only the years that carry real data (`prior_contraction ≥ 0.10`).
**How** — the same `rma()` refit on the informed subset; `NA` when fewer than 2 informed years remain.
**Tells you** — whether apparent stability is genuine or manufactured. **If this is much larger than `tau2`, the small τ² came from prior-dominated years all sitting near `b = 1` and pulling the series together — the prior talking, not a stable bias term.** Mostly a trailer concern, since trailer counts are often absent. Where it equals `tau2` exactly, every year was informed and no such artifact exists.
**Example** — `0.00284`, identical to `tau2` — all five Snohomish years are informed, in both likelihood types.

---

## Diagnostics

### `note`
**What it is** — why a row has no variance statistics, when it doesn't.
**How** — set by `fit_series()`. Two possible values:
- `single year — no interannual variance estimable` — only one usable year, so there is nothing to compare it against; `pooled_b` is just that year's estimate.
- `model did not converge` — `rma()` failed.
**Tells you** — that the blank columns are structural rather than a bug. `NA` here means the fit succeeded normally.
**Example** — `NA`
