# Thesis result figures — build directory

Seven publication figures for the UCL dissertation results chapter. Each figure
is one self-contained `fig-*.R` that sources `_common.R`, reads the canonical
analysis files, asserts the numbers it is about to draw, and writes a vector PDF
plus a 300 dpi PNG into this directory.

## Rebuild

```bash
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
cd /tmp/figures-build
for f in fig-primary fig-forest fig-complexity fig-engagement fig-strata fig-cas fig-actionsize; do
  /Users/wolpfor/miniconda3/envs/stats/bin/Rscript "$f.R"
done
```

Any `stopifnot` failure stops that figure's script before it writes output, so a
silent drift between the figures and the analysis pipeline is not possible.

## Conventions shared by all seven

| Item | Value |
|---|---|
| Palette | Okabe-Ito. reasoning `#0072B2`, control `#E69C00`, third level `#009E73` |
| Theme | `theme_minimal(base_size = 10)`, minor grid removed, no in-figure title |
| Intervals | **90% CI throughout**, stated on the axis label or in the figure's own note. No significance stars, no p values, no "significant" anywhere |
| Bootstrap | participant cluster bootstrap, `seed = 20260727`, `B = 2000`, percentile — the same estimator and seed as `explore_family.R` |
| Output | PDF (vector, width <= 5.5 in) + PNG at 300 dpi |

## Data sources

```
D   = /Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道
      $D/data/real-20260727/merged_cells_screened.csv          387 cells, 36 participants
      $D/data/real-20260727/raw/ReasoningEvent.csv             panel telemetry (dwell)
      $D/exploratory-注册族-20260729/participant_summary.csv    participant sufficient statistics
      $D/exploratory-注册族-20260729/cells_with_complexity.csv  387 cells + complexity tier
      $D/exploratory-注册族-20260729/reasoning_blocks.csv       273 panels, open flag
      $D/exploratory-注册族-20260729/results-exploratory.md     the reference numbers asserted against
```

Everything under `$D` is read-only. Nothing is written outside `/tmp/figures-build`.

---

## fig-primary

**Claim.** The raw arm gap (+8.1 pp) is small next to the between-participant
spread, and the two arms overlap heavily.

**Data.** `participant_summary.csv` (36 points, sized by `n_cells`),
`merged_cells_screened.csv` (cell-level rates).

**Assertions.**
- cell-level follow rate: reasoning `0.7740`, control `0.6927`
- participant-mean follow ratio: reasoning `0.7417`, control `0.6648`
- 18 participants per arm

Arm means carry a participant bootstrap 90% CI computed here:
control `[0.553, 0.764]`, reasoning `[0.645, 0.831]`.

## fig-forest

**Claim.** Every interval that adjusts for participant clustering contains OR = 1.

**Data.** Five **locked values** taken from the figure specification — the GLMM,
Firth and CR1 fits live in the confirmatory pipeline and are not refitted here.
They are written as an explicit constant table at the top of the script with the
source noted in the header comment. The descriptive twin of row 1 that *is*
reproducible from the exploratory family is `+0.081 [-0.043, 0.202]`
(results-exploratory.md, T7 ITT).

**Assertions.**
- five rows
- all four cluster-adjusted intervals straddle 1
- the fallback row (`1.51 [1.02, 2.24]`) does not, and is drawn with an open
  point and a dashed bar, with a note naming it as a different estimand with no
  cluster adjustment

## fig-complexity

**Claim.** The point-estimate ordering follows the design hinge (offsetting
largest), and every arm-difference interval contains 0.

**Data.** `cells_with_complexity.csv` for the rates; `participant_summary.csv`
for the bootstrap.

**Assertions.**
- rates: `.826/.789` single-factor, `.724/.582` two-factor-offsetting,
  `.750/.679` three-factor-mixed
- arm differences recomputed with the same estimator, seed and B as
  `explore_family.R` T2, and asserted equal to results-exploratory.md:
  `+0.037 [-0.109, 0.172]`, `+0.142 [-0.035, 0.331]`, `+0.071 [-0.093, 0.242]`

## fig-engagement

**Claim.** Voluntary re-opening was uncommon; dwell time was right-skewed, with
a median of 11.9 s and an IQR of 5.3--46.7 s.

**Data.** `reasoning_blocks.csv` (panel level); `raw/ReasoningEvent.csv`
`collapse` events carrying `dwell_ms`, restricted to the 36-participant set.

**Assertions.**
- 273 panels, 51 voluntarily reopened, overall `0.187`; the denominator (panel renders) is
  printed as `k/n` above every bar and again in the in-panel note
- dwell: 52 re-open--collapse pairs, median `11896 ms`, IQR `[5263, 46749]`

## fig-strata

**Claim.** As engagement rises, the perceived "it explained itself" reading
climbs while the follow rate does not.

**Data.** `participant_summary.csv`. Strata: control / reasoning never reopened /
reasoning reopened >= 1 panel. `MC1` is the H1 perceived-reasoning item, carried
into `participant_summary.csv` by `explore_family.R` (see deviation D1 in
报告.md).

**Assertions.**
- follow rate `0.693 / 0.795 / 0.758`, group n = `18 / 9 / 9` participants
- MC1 mean `2.444 / 3.875 / 6.222`, group n = `18 / 8 / 9`

Each group's 90% CI is bootstrapped on its own participants.

## fig-cas

**Claim.** Descriptively, the share of available capital committed on a follow
cell is higher in the reasoning arm.

**Data.** `merged_cells_screened.csv`, follow cells with a non-missing `cas`.

**Assertions.** n = `108 / 91` cells; mean `.505 / .315`; median `.441 / .229`.
The caption states the denominator as "CAS is defined only on cells adopted in the follow direction."

## fig-actionsize

**Claim.** Descriptively, positive capital moves are larger in the reasoning arm.

**Data.** `merged_cells_screened.csv`, cells with `invested > 0` (the primary
denominator).

**Assertions.** n = `161 / 124` cells; mean `460.66 / 333.08`;
median `300.00 / 277.33`. The caption states the denominator as "every executed capital move with amount > £0, including adopted step-back moves."
