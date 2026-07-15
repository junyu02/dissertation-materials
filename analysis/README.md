# Follow-Ratio analysis pipeline (B2)

Confirmatory stats for the 1×2 between-subjects study (`llm-reasoning` vs
`llm-control`). Fits the pre-registered covariate-adjusted model to the
**direction-aware** Follow Ratio and writes LaTeX `\newcommand` macros the
thesis pulls in with `\input`. The numbers in the write-up are never typed by
hand — they come from `results-values.tex`.

**Direction-aware DV (Plan X).** A cell counts as *followed* (`followed==1`)
when the participant takes the direction-compliant adoption: a **follow**
recommendation adopted by investing (`decision_type=='follow'`), OR an
**unfollow** recommendation adopted by divesting (`decision_type=='unfollow'`).
`direction` labels the *advice*; `decision_type` labels the *choice*, so a
declined unfollow keeps `direction=='unfollow'` with `decision_type=='decline'`.

The model is covariate-adjusted:
`followed ~ condition + literacy + is_HCI + factor(month) + factor(round) + (1|participant_id)`.

## Run

The **synthetic** generator produces data in the final merged schema (literacy /
is_HCI baked in), so it feeds `fit_glmm.R` directly. **Real** data arrives as two
exports (behavioural cells + questionnaire) and needs `prepare_cells.R` first to
derive and merge the questionnaire covariates.

```bash
# --- synthetic path (2 steps) ---
python3 synth_data.py                                   # -> data/synth_cells.csv (+ self-check)
/Users/wolpfor/miniconda3/envs/stats/bin/Rscript fit_glmm.R   # -> results-values.tex

# --- real-data path (3 steps) ---
Rscript prepare_cells.R study_cells.csv questionnaire_items.csv merged.csv   # derive+merge literacy/is_HCI
Rscript fit_glmm.R merged.csv results-values.tex                             # same step 2
#   then in the thesis preamble: \input{results-values.tex}  -> \FRORpoint, \FRORlow, ...
```

To exercise the fallback branch without a real non-convergence:
`FORCE_FALLBACK=1 Rscript fit_glmm.R`.

`prepare_cells.R` is **fail-loud**: a cells participant with no pre-task
questionnaire submission `stop()`s the run (never a silent drop). Orphan
questionnaire participants — a `participant_id` in the questionnaire but not in
the cells export (they carry an empty `condition`) — are the opposite case and
are correctly dropped by the cells-anchored join, not an error.

## CSV schema (the merged file `fit_glmm.R` consumes)

One row per decision cell. Synthetic: 40 participants × 20 cells = 800 rows.
Real: whatever the export contains (the test export is 62 cells). These are the
**required** columns; any extra columns in the real export (e.g.
`chatbot_change`, `newspost_change`, `cell_index`) are ignored by name.

| column | type | semantics |
|---|---|---|
| `participant_id` | string | unit of randomisation; each id is in exactly one `condition` |
| `condition` | `llm-reasoning` \| `llm-control` | between-subjects arm |
| `month` | int | study month; enters the model as `factor(month)` fixed dummies |
| `round` | int 1–4 | round within month; enters the model as `factor(round)` fixed dummies (§3.9 non-stationarity/learning clause, adopted 2026-07-15) |
| `round` | int | round within month |
| `decision_type` | `follow` \| `unfollow` \| `decline` \| `zero_optout` | what the participant did (see below) |
| `direction` | `follow` \| `unfollow` | direction of the **advice** (not the choice) |
| `followed` | 0/1 | **primary DV**; `1` iff `decision_type ∈ {follow, unfollow}` |
| `available_snapshot` | float (GBP) | budget available at decision time; **always present** |
| `invested` | float (GBP) or empty | amount committed/divested — see structural-NA rules |
| `cas` | float (0,1] or empty | conditional-on-follow share = `invested / available_snapshot`; **defined only when `decision_type=='follow'`** |
| `literacy` | int 0–3 | Lusardi Big-Three financial-literacy score; **participant-level** covariate |
| `is_HCI` | 0/1 | HCI / CS background (questionnaire F3); **participant-level** covariate |

On the real path `literacy` and `is_HCI` are not in the behavioural export — they
are derived from the questionnaire and merged in by `prepare_cells.R`.

### Three structural-NA types (design encoding, not missing data)

`cas` is the *conditional-on-follow* share, so it is deliberately NULL on every
row that is not a `follow`. The three non-follow-`cas` rows carry **different
`invested` encodings** so they stay separable in analysis:

| `decision_type` | `followed` | `invested` | `cas` | meaning |
|---|---|---|---|---|
| `unfollow` | 1 | **positive** (divestment amount; may exceed `available_snapshot`) | empty | divest advice **adopted** — an adoption, not a refusal |
| `decline` | 0 | **empty (NULL)** | empty | "Never mind": did not engage |
| `zero_optout` | 0 | **`0.00`** (a real recorded zero) | empty | engaged but committed **£0** |

`unfollow`'s `invested` is a withdrawal amount, so a *share of available* is
undefined — that is why `cas` is NULL even though it is an adoption. Its
`invested` can legitimately exceed `available_snapshot`; `fit_glmm.R` only checks
`cas == invested/available` on `follow` rows.

### True missing = a bug, and the pipeline aborts on it

A NULL `cas` on a `follow` row (or a non-NULL `cas` on any non-follow row, or
`followed` disagreeing with `decision_type`, or a `decline` with a non-NULL
`invested`, or a `zero_optout` without `invested == 0`, or an `unfollow` without
a positive `invested`, or an adoption whose `direction` contradicts its
`decision_type`) is a **data-integrity alarm**. `fit_glmm.R` counts these up
front and `stop()`s with the offending counts rather than silently modelling
corrupt rows. In clean behavioural data the true-missing count is `0` (reported
as `\FRnTrueMissing`).

## What `fit_glmm.R` does

1. **Audit** — separate structural NA from true missing; abort on any integrity alarm.
2. **Primary GLMM (confirmatory)** —
   `glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) + (1|participant_id), family=binomial, nAGQ=1)`
   (Laplace). Reports OR = `exp(β_condition)` with a **90% Wald CI** (primary)
   and a **90% profile-likelihood CI** (sensitivity). `literacy` and `is_HCI` are
   pre-registered participant-level covariates; `month` enters as **fixed**
   `factor(month)` dummies (not a second random effect) because it has < 8 levels
   — the pre-registered low-cardinality rule (outline §3.9).
3. **Coded convergence gate** — a joint criterion; all three required to accept
   the GLMM: (a) lme4 emits **no convergence message**, (b) the **relative
   gradient** `max(abs(solve(H, g)))` — the metric from lme4's troubleshooting
   doc — is **≤ 0.002** (threshold borrowed from `checkConv`'s default tol), and
   (c) the **optimizer return code is 0** (a channel not carried in `checkConv`'s
   messages). Joint strictness is therefore ≥ lme4's own default checks. Any
   failure — or a hard fit error — routes to the fallback. No post-hoc judgement
   call; the criterion is in the code.
4. **Pre-registered fallback (Firth)** — collapse to per-participant follow
   counts and fit a Firth-penalised logistic
   (`logistf(y ~ condition + literacy + is_HCI, weights=w)`, α=0.10). `literacy`
   and `is_HCI` survive the collapse (participant-level); the `factor(month)` / `factor(round)`
   dummies **cannot** (they vary within participant) — dropping them is a known,
   documented loss of this degraded path. Firth
   converges even under separation, so this is the convergence-guaranteed last
   resort. On this path the **primary CI (`\FRORlow`/`\FRORhigh`) is logistf's
   penalised-profile CI** (the recommended interval under Firth); Wald is demoted
   to `\FRORlowWald`/`\FRORhighWald`. **Locked commitment:** if the real data
   triggers this fallback, the write-up additionally reports a **Wilcoxon
   rank-sum on per-participant follow proportions** as a clustering sensitivity
   check, and states plainly that the pooled Firth CI is **not** corrected for
   within-subject correlation (anti-conservative). Fixed now to remove any
   writing-period degree of freedom. The taken path is recorded in `\FRfitPath`.
5. **Secondary descriptives (CAS)** — per-arm mean / median / IQR of the
   conditional-on-follow share. **No model is fit to CAS** (no Tobit, per spec).
6. **Emit** `results-values.tex`.

Both fit paths emit the **same macro set**, so the thesis compiles regardless
of which one ran. **The primary CI (`\FRORlow`/`\FRORhigh`) is path-dependent —
Wald on the GLMM path, penalised-profile on the Firth path — disambiguated by
`\FRfitPath`.** `\FRORlowWald`/`\FRORhighWald` always carry the Wald pair so a
Wald interval is available on either path.

> **Fallback caveat.** `logistf` on participant-collapsed counts reconstructs a
> pooled penalised likelihood; its CI is **not** adjusted for between-participant
> variance the way the GLMM random intercept is (expect a narrower CI). It is a
> stability escape hatch, not an equal substitute. If it ever fires on real
> data, say so in the write-up and treat its CI as approximate. Prefer fixing
> the convergence problem (more data, simpler RE structure) over shipping the
> fallback.

## Macros in `results-values.tex`

Names are all-letters (LaTeX-legal). Values may contain digits.

| macro | meaning |
|---|---|
| `\FRfitPath` | `GLMM (Laplace)` or `Firth (participant-collapsed)` |
| `\FRnTotal` | participants analysed (40) |
| `\FRnCells` | decision cells analysed (800) |
| `\FRORpoint` | OR, reasoning vs control |
| `\FRORlow` / `\FRORhigh` | 90% **primary** CI — Wald (GLMM) or penalised-profile (Firth); see `\FRfitPath` |
| `\FRORlowProfile` / `\FRORhighProfile` | 90% profile CI (sensitivity) — ordinary likelihood profile on GLMM, penalised profile on Firth |
| `\FRORlowWald` / `\FRORhighWald` | 90% Wald CI (both paths, for parity) |
| `\FRgradMax` | max scaled gradient (GLMM path; `NA` on fallback) |
| `\FRrateReasoning` / `\FRrateControl` | raw adoption rate per arm (`mean(followed)`, both directions) |
| `\CASmeanReasoning` / `\CASmeanControl` | CAS mean per arm |
| `\CASmedianReasoning` / `\CASmedianControl` | CAS median per arm |
| `\CASqOneReasoning` / `\CASqThreeReasoning` | CAS Q1 / Q3 (reasoning) |
| `\CASqOneControl` / `\CASqThreeControl` | CAS Q1 / Q3 (control) |
| `\FRnUnfollow` | count of `unfollow` adoptions (divest advice adopted) |
| `\FRnDecline` / `\FRnOptout` | count of `decline` / `zero_optout` rows |
| `\FRnStructuralNA` | rows with structural NULL cas = all non-`follow` rows (`decline` + `zero_optout` + `unfollow`) |
| `\FRnTrueMissing` | `follow` rows with NULL cas — must be 0 |

## Analysis recipes (secondary DVs — derived, not in `study_cells`)

These are **not** columns in the behavioural export; each is reconstructed from
the raw platform models. Where a join is heuristic, the heuristic and its limit
must be stated in the write-up and fixed in the pre-registration.

1. **Decision modality (typed vs clicked)** — join `Message.from_button` onto
   each decision cell by a **time-window** match (`Message.from_button` has no
   `round` column, so pair the button/typed message to the cell whose decision
   timestamp it is nearest). Heuristic; state the window and that ambiguous
   ties are dropped, not guessed. Feeds §5.3 / §5.7.
2. **Do-it-then-disconnect sensitivity** — rebuild the affected cells with a
   `Recommendation ⟕ UserAction` **outer** join: a `Recommendation` with a
   terminal decision but no paired invest `UserAction` and no following `Round`
   is a disconnect. Re-run the primary model with these cells excluded as a
   robustness check (§5.1).
3. **Questionnaire de-dup** — one submission per participant **per phase**: keep
   the **earliest** `created_at` submission (a re-take counts once). This is the
   rule `prepare_cells.R` already applies to `pretask`; apply the identical rule
   to `posttask` scales.
4. **Time-on-task proxy** — `Month` has no end timestamp, so use
   `updated_at` / `ended_reason` as the close signal, or difference consecutive
   `Round` / `UserAction` `created_at` within a month. It is a proxy; name which
   signal was used and its granularity limit.
5. **Offsetting vs non-offsetting classification** — not stored; reconstruct
   from `FactorWeight` + `FrozenChange` (whether the round's factor move opposes
   the standing position). Because it is derived, the exact rule must be
   **fixed explicitly in the pre-registration** before it touches any DV.

## Hardening (R2 cross-check, 2026-07-15)

Three independent AI reviewers + a supervisor pass added, and `test_alarms.R`
now proves, the following fail-loud guards (12 corrupt-input cases, each
asserted to stop the pipeline; a 13th control case asserts clean data still
passes — run `Rscript test_alarms.R` after `synth_data.py`):

- enum whitelists + NA guards on `decision_type` / `direction` (data.table
  treats NA conditions as FALSE, so an NA row previously slipped every alarm)
- `(participant_id, month, round)` uniqueness (a doubled export would silently
  double-count cells inside the GLMM)
- domain checks BEFORE integer coercion (`followed=1.7` no longer truncates to 1)
- direct `invested > 0` on follow rows (previously only implied via the cas
  identity, which an NA invested NA-skipped); `available_snapshot >= 0`
  (0 is legal — a fully-invested participant declining; negative/NA/Inf stop)
- `zero_optout` on an unfollow-direction cell = contract-undefined combination
  → alarm for manual adjudication
- prepare_cells.R: raw ISO timestamp preserved via `colClasses` (fread's
  auto-POSIXct would silently re-format and break the lexicographic de-dup),
  ISO-format + single-UTC-offset assertions, exactly-7-items-per-participant
  assertion, one-item-per-question assertion, `nchar(D1_CORRECT)==14` mojibake
  trap, and a D1 zero-pass-rate warning
- `cas` precision contract: the exporter writes cas at 10dp so the 1e-9
  consistency tolerance holds; if the export precision ever changes, revisit
- singular/boundary GLMM fits are deliberately NOT a fallback trigger (they
  degenerate to a pooled logistic; computed and printed for the record) —
  the coded gate is: no lme4 message, relative gradient <= 0.002, optimizer
  code 0
- profile-CI warnings are logged (not silently muffled)

## Environment

conda env `stats`, pinned in `environment.yml` (R 4.3.3, lme4 1.1.37,
logistf 1.26.1, data.table 1.17.8 on conda-forge; the freeze also had
Matrix 1.6-5, which lme4 1.1.37 is compatible with).

```bash
conda env create -n stats -f environment.yml
```

Python is stdlib-only (any Python 3.8+); `synth_data.py` has zero third-party
dependencies.
