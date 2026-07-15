#!/usr/bin/env Rscript
# Step 1 of the REAL-data pipeline: merge pre-task questionnaire covariates
# (financial literacy, is_HCI) onto the behavioural decision cells, producing the
# same schema fit_glmm.R (step 2) consumes.
#
#   Rscript prepare_cells.R <study_cells.csv> <questionnaire_items.csv> <out_merged.csv>
#   Rscript fit_glmm.R      <out_merged.csv>  <results-values.tex>
#
# Synthetic data skips this step -- synth_data.py already emits the merged schema
# with literacy / is_HCI baked in.
#
# Style matches fit_glmm.R: data.table, fail-loud at every boundary. A cells
# participant with no questionnaire submission is a data problem and stops the
# run; it is never silently dropped. Orphan questionnaire rows (a participant in
# the questionnaire but not in the cells export) are the opposite case and are
# correctly dropped by the cells-anchored join.

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: Rscript prepare_cells.R <study_cells.csv> <questionnaire_items.csv> <out_merged.csv>")
cells_csv <- args[[1]]; quest_csv <- args[[2]]; out_csv <- args[[3]]

# encoding="UTF-8": the D1 correct answer contains a GBP sign (multibyte). We
# compare answers with `==`, which is byte-wise on UTF-8 strings and therefore
# locale-independent, but the file must still be READ as UTF-8 or the pound sign
# decodes wrong and silently zeroes every D1 score.
cells <- fread(cells_csv, na.strings = "", encoding = "UTF-8")
# created_at MUST stay the raw ISO string: fread's auto-detection would parse it
# to POSIXct and as.character() would re-format it WITHOUT the T/offset, breaking
# both the ISO assertion and (worse, silently) the lexicographic-min de-dup.
quest <- fread(quest_csv, na.strings = "", encoding = "UTF-8",
               colClasses = list(character = "created_at"))
stopifnot(nrow(cells) > 0, nrow(quest) > 0)
for (col in c("participant_id", "phase", "question", "answer", "created_at"))
  if (!col %in% names(quest)) stop(sprintf("questionnaire missing column: %s", col))
if (!"participant_id" %in% names(cells)) stop("cells missing column: participant_id")

## ---- 1. pretask only -------------------------------------------------------
pre <- quest[phase == "pretask"]
if (nrow(pre) == 0) stop("no pretask questionnaire rows found")

## ---- 2. de-dup: one submission per participant (pre-registered rule) --------
# A participant who re-took the pre-task has multiple submissions (distinct
# created_at). Pre-registered口径: keep the EARLIEST submission per participant.
# created_at is an ISO-8601 string with a fixed +00:00 offset, so lexicographic
# min == chronological min; no date parsing needed.
pre[, created_at := as.character(created_at)]
# Lexicographic min == chronological min ONLY for ISO-8601 strings sharing one
# UTC offset. Assert both before relying on it (R2 cross-check finding): a
# malformed or mixed-offset timestamp must stop the run, not silently reorder it.
stopifnot(all(grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}", pre$created_at)))
if (uniqueN(sub("^.*(\\+\\d{2}:\\d{2}|Z)$", "\\1", pre$created_at)) > 1)
  stop("mixed created_at UTC offsets: lexicographic de-dup would be wrong")
first_sub <- pre[, .(created_at = min(created_at)), by = participant_id]
pre <- pre[first_sub, on = .(participant_id, created_at)]
# The de-dup assumes one submission writes ALL its items with one shared
# timestamp. Assert the deployed shape (7 pretask items per participant) so a
# platform change to per-row timestamps, or a truncated first submission,
# fails loud instead of silently scoring missing items as 0.
item_counts <- pre[, .N, by = participant_id]
if (item_counts[N != 7, .N] > 0)
  stop(sprintf("pretask item count != 7 after de-dup for: %s",
               paste(item_counts[N != 7, participant_id], collapse = ", ")))

## ---- 3. score literacy (Lusardi Big-Three) + is_HCI ------------------------
# Correct answers hard-coded from the deployed questionnaire wording. Only an
# exact-correct answer scores 1; "Do not know" / "Prefer not to say" / blank all
# score 0 (a wrong-or-absent answer is not financial literacy).
#   D1 (compound interest)          -> "More than £102"
#   D2 (real return under inflation) -> "Less than today"
#   D3 (single stock vs fund, T/F)   -> "False"
#   F3 (HCI / CS background, yes/no) -> "Yes"  => is_HCI = 1
D1_CORRECT <- "More than £102"   # £ = literal GBP sign (U+00A3). Correctness
                                      # needs BOTH this source file and the data read as UTF-8
                                      # (fread encoding="UTF-8" above); == is byte-wise, so locale
                                      # cannot silently zero D1. The real-data run log below prints
                                      # per-participant literacy -- eyeball it against known values.
D2_CORRECT <- "Less than today"
D3_CORRECT <- "False"
F3_YES     <- "Yes"
# Mechanical mojibake trap (R2 cross-check finding): if this source file is ever
# re-saved in a non-UTF-8 encoding, the pound sign becomes two bytes ("Â£") and
# every D1 answer would silently score 0. nchar() catches that immediately.
stopifnot(nchar(D1_CORRECT) == 14L)

ans_for <- function(pid, prefix) {
  v <- pre[participant_id == pid & startsWith(question, prefix), answer]
  if (length(v) > 1)
    stop(sprintf("participant %s has %d '%s' items in one submission (expected 1)",
                 pid, length(v), prefix))
  if (length(v) == 0) NA_character_ else v[[1]]
}

pids <- unique(pre$participant_id)
lit_dt <- rbindlist(lapply(pids, function(pid) {
  d1 <- ans_for(pid, "D1."); d2 <- ans_for(pid, "D2."); d3 <- ans_for(pid, "D3.")
  f3 <- ans_for(pid, "F3.")
  data.table(
    participant_id = pid,
    literacy = sum(identical(d1, D1_CORRECT),
                   identical(d2, D2_CORRECT),
                   identical(d3, D3_CORRECT)),
    is_HCI   = as.integer(identical(f3, F3_YES))
  )
}))

## ---- 4. merge onto cells (cells-anchored) ----------------------------------
# A cells participant with no pretask submission is a hard error, not a silent
# drop. Orphan questionnaire participants (in lit_dt but not in cells) are simply
# not matched by the left join.
missing <- setdiff(unique(cells$participant_id), lit_dt$participant_id)
if (length(missing))
  stop(sprintf("cells participant(s) missing a pretask questionnaire submission: %s",
               paste(missing, collapse = ", ")))

merged <- merge(cells, lit_dt, by = "participant_id", all.x = TRUE, sort = FALSE)
stopifnot(nrow(merged) == nrow(cells),
          !anyNA(merged$literacy), all(merged$literacy %in% 0:3),
          !anyNA(merged$is_HCI),   all(merged$is_HCI %in% c(0L, 1L)))

fwrite(merged, out_csv)

## ---- 5. run log ------------------------------------------------------------
cells_lit <- unique(merged[, .(participant_id, literacy, is_HCI)])[order(participant_id)]
n_orphan  <- length(setdiff(lit_dt$participant_id, cells$participant_id))
cat(sprintf("[prepare_cells] %d cells rows, %d participants merged -> %s\n",
            nrow(merged), nrow(cells_lit), out_csv))
cat(sprintf("[prepare_cells] dropped %d orphan questionnaire participant(s) (no cells rows)\n", n_orphan))
cat("[prepare_cells] per-participant literacy / is_HCI:\n")
print(cells_lit)
# D1 is the easiest Lusardi item; a 0% pass rate across several participants is
# far more likely an encoding / wording-drift bug than real data. Warn, don't stop.
n_d1_correct <- sum(sapply(cells_lit$participant_id, function(pid)
  identical(ans_for(pid, "D1."), D1_CORRECT)))
if (nrow(cells_lit) >= 3 && n_d1_correct == 0)
  cat("[prepare_cells] WARNING: 0 participants scored D1 correct — check encoding / questionnaire wording drift\n")
