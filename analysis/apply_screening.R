#!/usr/bin/env Rscript
# apply_screening.R -- careless-responder exclusion (AsPredicted #301918, Q6(a)).
#
# WHAT THIS IS
#   Q6(a) of the registration commits to excluding careless responders identified
#   by a straight-lining screen on the post-task questionnaire battery. This
#   script executes that committed exclusion and emits the screened cell file
#   that the frozen pipeline (prepare_cells.R -> fit_glmm.R) then consumes.
#
# THRESHOLD PROVENANCE (read this before citing the numbers)
#   The registration commits to a straight-lining screen but DOES NOT FIX A
#   THRESHOLD. This script therefore adopts the STRICTEST defensible reading --
#   a participant is careless only if every numeric Likert response in the
#   post-task battery carries the same value (number of distinct values == 1),
#   i.e. a completely flat response vector. A relaxed reading (distinct values
#   <= 2) is also computed and PRINTED as a sensitivity reference, but is NOT
#   applied. That relaxed reading flags FOUR ADDITIONAL participants who are
#   present in the cell file (UCL3GNOX0FG, UCLGVXX8SZQ, UCLSPQP3TAQ,
#   UCLYFQMGG1H); applying it would take N from 36 to 32.
#
#   WHY THE STRICT READING IS THE DEFENSIBLE ONE -- not because it excludes
#   the fewest people (fewest exclusions is not the same thing as fewest
#   researcher degrees of freedom), but because the relaxed reading misfires on
#   genuine respondents: UCLGVXX8SZQ answers with {2, 7}, two OPPOSITE extremes,
#   which is the signature of an engaged respondent rather than a straight-liner.
#   A screen that cannot tell those apart is a worse screen. Note also that the
#   primary estimate is insensitive to the choice: OR = 1.72 (no exclusion),
#   1.68 (strict), 1.53 (relaxed), every interval containing 1.
#
# NON-FROZEN, POST-WINDOW
#   This is a NON-FROZEN ADDITIONAL PROGRAM written after the data-collection
#   window closed (2026-07-27), on data that had already been inspected. It is
#   NOT part of the registered frozen analysis code (commit 352792c:
#   prepare_cells.R / fit_glmm.R / test_alarms.R), which it does not modify and
#   only feeds. The exclusion rule it implements is registered; the threshold
#   operationalisation and this file are not. Anything downstream that cites the
#   screened results must carry that distinction.
#
# ANTI-DRIFT GUARD
#   The exclusion set is asserted against the two participant IDs adjudicated
#   under the criterion above. If a future rerun (changed export, changed
#   criterion, changed data) produces a different set, the script STOPS rather
#   than silently re-defining who was careless.
#
# Usage:
#   Rscript apply_screening.R [study_cells.csv] [questionnaire_items.csv] \
#                             [out_screened_cells.csv] [out_screening_log.csv]

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(here) == 0 || here == "") here <- "."
dd <- file.path(here, "data", "real-20260727")
cells_csv <- if (length(args) >= 1) args[[1]] else file.path(dd, "study_cells.csv")
quest_csv <- if (length(args) >= 2) args[[2]] else file.path(dd, "questionnaire_items.csv")
out_cells <- if (length(args) >= 3) args[[3]] else file.path(dd, "study_cells_screened.csv")
out_log   <- if (length(args) >= 4) args[[4]] else file.path(dd, "screening_log.csv")

# The adjudicated exclusion set (see ANTI-DRIFT GUARD above).
EXPECTED_EXCLUDED <- c("UCL0DHECO8C", "UCL5N72OPCK")

## ---- 1. load ---------------------------------------------------------------
# encoding / created_at handling mirrors prepare_cells.R: created_at must stay a
# raw ISO string so lexicographic min == chronological min for the de-dup below.
cells <- fread(cells_csv, na.strings = "", encoding = "UTF-8")
quest <- fread(quest_csv, na.strings = "", encoding = "UTF-8",
               colClasses = list(character = "created_at"))
stopifnot(nrow(cells) > 0, nrow(quest) > 0)
for (col in c("participant_id", "phase", "question", "answer", "created_at"))
  if (!col %in% names(quest)) stop(sprintf("questionnaire missing column: %s", col))
if (!"participant_id" %in% names(cells)) stop("cells missing column: participant_id")

post <- quest[phase == "posttask"]
if (nrow(post) == 0) stop("no posttask questionnaire rows found")

## ---- 2. de-dup to the earliest submission per participant ------------------
# Same registered rule prepare_cells.R applies to the pre-task battery: a
# participant with several submissions keeps the earliest. Without it a
# double-submit would pool two response vectors and could mask a flat one.
post[, created_at := as.character(created_at)]
stopifnot(all(grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}", post$created_at)))
if (uniqueN(sub("^.*(\\+\\d{2}:\\d{2}|Z)$", "\\1", post$created_at)) > 1)
  stop("mixed created_at UTC offsets: lexicographic de-dup would be wrong")
n_multi <- post[, uniqueN(created_at), by = participant_id][V1 > 1, .N]
first_sub <- post[, .(created_at = min(created_at)), by = participant_id]
post <- post[first_sub, on = .(participant_id, created_at)]

## ---- 3. numeric Likert responses only --------------------------------------
# The battery mixes Likert items with free-text items (technical-problem report,
# open comments, the H2 multi-select). Selecting purely numeric answers isolates
# the Likert responses. The selected question set is printed so that a wording or
# response-format change upstream is visible rather than silent.
is_numeric_answer <- function(x)
  !is.na(x) & grepl("^[[:space:]]*-?[0-9]+(\\.[0-9]+)?[[:space:]]*$", x)
lik <- post[is_numeric_answer(answer)]
if (nrow(lik) == 0) stop("no numeric posttask Likert responses found")
lik[, value := as.numeric(answer)]

## ---- 4. straight-lining screen ---------------------------------------------
scr <- lik[, .(n_items = .N, n_unique = uniqueN(value),
               values = paste(sort(unique(value)), collapse = "|")),
           by = participant_id]

cells_pids <- unique(cells$participant_id)
scr[, in_cells := participant_id %in% cells_pids]

flag_strict  <- scr[n_unique == 1L, participant_id]                      # applied
flag_relaxed <- scr[n_unique <= 2L, participant_id]                      # reference only

# The screen only removes rows that exist: a straight-liner with no behavioural
# cells is not part of the analysis sample and is already dropped by the
# cells-anchored join in prepare_cells.R. Restricting here keeps the two facts
# separate instead of conflating "flagged" with "excluded".
excluded <- sort(intersect(flag_strict, cells_pids))

if (!identical(excluded, sort(EXPECTED_EXCLUDED)))
  stop(sprintf(paste0("EXCLUSION SET DRIFT: strict criterion (n_unique==1, restricted to ",
                      "participants present in the cells export) yields {%s}; adjudicated set is {%s}. ",
                      "Stopping rather than silently re-defining the careless-responder set."),
               paste(excluded, collapse = ", "), paste(sort(EXPECTED_EXCLUDED), collapse = ", ")))

## ---- 5. emit screened cells + screening log --------------------------------
screened <- cells[!participant_id %in% excluded]
stopifnot(nrow(screened) < nrow(cells),
          uniqueN(screened$participant_id) == length(cells_pids) - length(excluded))
fwrite(screened, out_cells)

# One row per participant across the union of the cells export and the post-task
# battery, so a cells participant with no post-task submission (unscreenable) and
# a post-task participant with no cells (unexcludable) are both on the record.
log_dt <- merge(
  data.table(participant_id = union(cells_pids, scr$participant_id)),
  scr[, .(participant_id, n_likert_items = n_items, n_unique_values = n_unique, values)],
  by = "participant_id", all.x = TRUE)
log_dt[, in_cells := participant_id %in% cells_pids]
log_dt[, n_cells := cells[log_dt, on = "participant_id", .N, by = .EACHI]$N]
log_dt[, has_posttask := !is.na(n_unique_values)]
log_dt[, flagged_strict  := !is.na(n_unique_values) & n_unique_values == 1L]
log_dt[, flagged_relaxed := !is.na(n_unique_values) & n_unique_values <= 2L]
log_dt[, excluded := participant_id %in% excluded]
setorder(log_dt, -flagged_strict, n_unique_values, participant_id)
fwrite(log_dt, out_log)

## ---- 6. run log ------------------------------------------------------------
arm_n <- function(d) {
  a <- unique(d[, .(participant_id, condition)])[, .N, by = condition]
  paste(sprintf("%s=%d", a$condition, a$N), collapse = " ")
}
cat("[screening] AsPredicted #301918 Q6(a) careless-responder exclusion\n")
cat(sprintf("[screening] posttask participants=%d (%d with >1 submission, earliest kept)\n",
            uniqueN(post$participant_id), n_multi))
cat(sprintf("[screening] numeric Likert questions detected: %d\n", uniqueN(lik$question)))
for (qq in sort(unique(lik$question))) cat(sprintf("             - %s\n", substr(qq, 1, 78)))
cat("\n[screening] STRICT criterion (n_unique == 1) -- APPLIED\n")
if (length(flag_strict) == 0) cat("             (none flagged)\n")
for (p in sort(flag_strict)) {
  r <- scr[participant_id == p]
  cat(sprintf("             %s  items=%d  unique=%d  values={%s}  in_cells=%s\n",
              p, r$n_items, r$n_unique, r$values, r$in_cells))
}
cat(sprintf("[screening] -> excluded (flagged AND present in cells): %s\n",
            paste(excluded, collapse = ", ")))
only_flagged <- setdiff(flag_strict, cells_pids)
if (length(only_flagged))
  cat(sprintf("[screening] -> flagged but NOT in the cells export (no behavioural data; already outside the analysis sample): %s\n",
              paste(sort(only_flagged), collapse = ", ")))

cat("\n[screening] RELAXED criterion (n_unique <= 2) -- REFERENCE ONLY, NOT APPLIED\n")
for (p in sort(flag_relaxed)) {
  r <- scr[participant_id == p]
  cat(sprintf("             %s  items=%d  unique=%d  values={%s}  in_cells=%s\n",
              p, r$n_items, r$n_unique, r$values, r$in_cells))
}

no_post <- setdiff(cells_pids, scr$participant_id)
if (length(no_post))
  cat(sprintf("\n[screening] cells participant(s) with NO posttask battery (unscreenable, retained): %s\n",
              paste(sort(no_post), collapse = ", ")))

cat(sprintf("\n[screening] BEFORE: %d cells, %d participants  (%s)\n",
            nrow(cells), uniqueN(cells$participant_id), arm_n(cells)))
cat(sprintf("[screening] AFTER : %d cells, %d participants  (%s)\n",
            nrow(screened), uniqueN(screened$participant_id), arm_n(screened)))
cat(sprintf("[screening] removed %d cells / %d participants\n",
            nrow(cells) - nrow(screened), length(excluded)))
cat(sprintf("[screening] wrote %s\n[screening] wrote %s\n", out_cells, out_log))
