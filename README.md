# Displayed Reasoning and Following Behaviour in AI Financial Advice — Analysis Materials

Analysis code and de-identified data for the MSc dissertation *Displayed
Reasoning and Following Behaviour in AI Financial Advice: A Controlled
Estimation Study* (Junyu Wei, UCL Interaction Centre, University College
London, 2026).

- **Experimental platform (source code):** https://github.com/junyu02/SocialTradingChatbot
- **Preregistration:** https://aspredicted.org/py622u.pdf

The study compared two arms of an LLM financial advisor — `llm-reasoning`
(advice plus an expandable reasoning block) and `llm-control` (identical
advice, no reasoning block) — on whether participants followed the advisor's
recommendation. The confirmatory analysis is a mixed-effects logistic
regression over decision cells nested within participants.

## Repository history

This repository has two commits, and the split is deliberate.

The **first commit** (`4b177e4`, dated 2026-07-15) is the registration-frozen
analysis pipeline, imported byte-identical from commit
`352792ce3a860a746a97c1443587e91e164eb37f` of the private project repository —
the hash recorded in the AsPredicted #301918 registration. It was frozen before
any outcome data was inspected. At that point `results-values.tex` had been
generated from the synthetic dataset, as its own header states.

It carries eight of the nine files in that frozen tree. The ninth,
`dev_smoke.csv`, is a development smoke-test fixture that no script in this
archive reads; it is withheld because it contains two platform-format
participant codes whose provenance could not be positively established.

The **second commit** is the final state: the post-window analysis, the
exploratory family, the de-identified data, and the figure programs.

So every change made to the analysis after the registration freeze is visible
in this repository's own history:

```
git diff 4b177e4..master -- analysis/
```

## Package structure

```
analysis/      Confirmatory pipeline (preregistered) and its integrity checks
exploratory/   Preregistered exploratory family: code, derived tables, results
posthoc/       Post-hoc robustness recomputation (optimiser sweep, bootstrap)
data/          De-identified analysis datasets and the LaTeX value macros
figures/       Figure programs for the seven result figures
LICENSE        MIT (code) + CC BY 4.0 (data)
```

### `analysis/`

| File | Role |
|---|---|
| `prepare_cells.R` | Step 1: joins pre-task questionnaire covariates (financial literacy, HCI background) onto the behavioural decision cells. |
| `fit_glmm.R` | Step 2: fits the confirmatory mixed-effects logistic model and emits the LaTeX value macros. |
| `test_alarms.R` | Fail-loud integrity checks on the step-2 alarm conditions. |
| `apply_screening.R` | Applies the preregistered attention/straight-lining screen; emits `screening_log.csv`. |
| `identify_disconnect_cells.R` | Post-window check for cells where a participant accepted a recommendation but disconnected before the capital move completed. |
| `apriori_precision_reconstruction.R` | Reconstructs the a-priori precision (CI half-width) expectation; self-contained, no data inputs. |
| `gender_sensitivity_20260730.R` | Post-hoc sensitivity adding participant gender to the registered covariate set. |
| `p2-computations-20260730.py` | Supporting computations over the exploratory derived tables. |
| `p2-computations-输出-20260730.txt` | Captured stdout of the above. |
| `synth_data.py` | Generates the synthetic dataset (`data/synth_cells.csv`) used to exercise the pipeline without real data. |

### `exploratory/`

`explore_family.R` computes the preregistered exploratory family;
`verify_explore_family.R` independently re-derives its numbers.
`results-exploratory.md` is the reported output. The four CSVs are the derived
tables the figure programs and `p2-computations` consume:

| File | Rows | Unit |
|---|---|---|
| `cells_with_complexity.csv` | 387 | decision cell, with advice-complexity tier |
| `participant_summary.csv` | 36 | participant (screened set) |
| `reasoning_blocks.csv` | 273 | reasoning block, with open/flash counts (18 participants in the reasoning arm) |
| `decision_latency.csv` | 411 | decision cell, unscreened, with latency in seconds |

### `posthoc/`

`recheck_allfit.R` refits the confirmatory model across every available
optimiser to confirm the estimate is not optimiser-dependent;
`recheck_allfit_stdout.txt` is its captured output.
`recheck_task3_variants.R` refits model variants.
`cluster_bootstrap_canonical.R` and `boot_seed.R` are the shared
participant-cluster bootstrap helpers (fixed seed, sourced by the callers).

### `data/`

| File | Rows | Participants | Contents |
|---|---|---|---|
| `study_cells.csv` | 411 | 38 | Decision cells before screening |
| `study_cells_screened.csv` | 387 | 36 | Decision cells after the preregistered screen |
| `merged_cells.csv` | 411 | 38 | As above, plus questionnaire covariates |
| `merged_cells_screened.csv` | 387 | 36 | **Primary analysis dataset** |
| `screening_log.csv` | 40 | 40 | Per-participant screening decision and its inputs |
| `disconnect_cells.csv` | 0 | 0 | Header only — no disconnected cells were found |
| `Round.csv` | 733 | 51 | Round-level event log (participant code, month, round index, timestamp) |
| `synth_cells.csv` | 800 | 40 (synthetic) | Synthetic dataset in the same schema, IDs `P01`–`P40`; contains no real participant data |
| `results-values.tex` | — | — | LaTeX value macros, unscreened |
| `results-values-screened.tex` | — | — | LaTeX value macros, screened (reported in the manuscript) |

`Round.csv` spans 51 participant codes because it logs everyone who began the
task, including those who did not complete it or were screened out; the
analysed sets are the 38 / 36 above.

### `figures/`

`_common.R` holds the shared theme, colour palette (Okabe–Ito, colour-blind
safe), data loaders and bootstrap helpers; each `fig-*.R` sources it and writes
one figure. `figures/README.md` documents the build. Rendered PDF/PNG outputs
are not included — the programs regenerate them.

## De-identification

Participants appear only as pseudonymous codes of the form `UCL` followed by a
random alphanumeric string (e.g. `UCL0ZGVZV11`). These codes were generated by
the platform and carry no information about the person.

No direct identifiers are present anywhere in this package: no names, email
addresses, IP addresses, recruitment-panel or crowdsourcing platform IDs,
dates of birth, or institutional affiliations of participants. The code-to-person mapping is not
part of this release and is not derivable from it.

All free-text material has been withheld. Specifically, the following source
tables were deliberately excluded and are **not** in this package:

- free-text chat transcripts between participant and advisor
- free-text questionnaire responses
- the prize-draw entry table and the platform's user-account table (both
  contained email addresses)
- participant and profile registry tables
- the platform database file

The `profile_name` column in `cells_with_complexity.csv` (e.g. `Aricka Lewis`)
names one of ten **fictional trader personas** authored as experimental
stimuli. These are not real people and not participants.

## Reproducibility scope

Withholding the free-text and raw source tables means the package reproduces
the pipeline from the analysis datasets onward, not from the raw platform
export. Concretely:

**Reproducible from this package:**

- the confirmatory model and its LaTeX value macros —
  `Rscript analysis/fit_glmm.R data/merged_cells_screened.csv out.tex`
- the same model on the bundled synthetic dataset, with no arguments —
  `Rscript analysis/fit_glmm.R` defaults to `data/synth_cells.csv`, so the
  pipeline runs out of the box
- the integrity-check suite — `Rscript analysis/test_alarms.R` runs against the
  bundled synthetic dataset, so the fail-loud alarm conditions are verifiable
  out of the box
- the post-hoc optimiser sweep and model variants (`posthoc/`)
- the a-priori precision reconstruction (`analysis/apriori_precision_reconstruction.R`)
- `analysis/p2-computations-20260730.py` over the `exploratory/` tables
- six of the seven figures (see caveat below)

**Not reproducible from this package**, because their inputs are withheld:

- `analysis/apply_screening.R` and `analysis/gender_sensitivity_20260730.R`
  read the free-text questionnaire table
- `analysis/prepare_cells.R` and `analysis/identify_disconnect_cells.R` read
  the raw platform export
- `exploratory/explore_family.R` and `verify_explore_family.R` read the raw
  export plus the platform's frozen stimulus files
- the dwell-time loader in `figures/_common.R` reads the raw reasoning-event
  table, so the engagement figure's dwell panel cannot be rebuilt

Their committed outputs are included so the reported numbers can be checked
even where the step itself cannot be rerun.

**Paths:** the scripts are released exactly as run, so several carry absolute
paths from the author's machine (`figures/_common.R`, `posthoc/recheck_*.R`,
`exploratory/*.R`, and the `/tmp/figures-build` output directory in
`figures/fig-*.R`). Adjust the `D` / `DIR` / `OUT` constants at the top of
those files before running. Nothing was rewritten, so what is published is the
code that produced the reported results.

**Environment:** R with `lme4`, `data.table`, `ggplot2`, `scales`; Python 3
with `pandas`. The captured `posthoc/recheck_allfit_stdout.txt` records the
exact R and package versions used.

## Licence

- **Code** (`.R`, `.py`): MIT
- **Data and outputs** (`.csv`, `.tex`, `.txt`, `.md`): CC BY 4.0

See `LICENSE` for both texts and a suggested citation.
