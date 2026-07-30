#!/usr/bin/env Rscript
# Post-window ADD-ON (NOT part of the frozen registered pipeline 352792c).
#
# Spec: POSTWINDOW-非冻结附加程序-断线格识别-20260726.md  (this directory)
#       — "识别规则（唯一判据，正向白名单）", written 2026-07-26, approved same day.
# Implements the registered Q6 coding rule the frozen exporter cannot express:
# a "Do it"-then-disconnect cell is indistinguishable from a genuine "Never mind"
# in study_cells.csv, because the DECISION VERB never reaches the export
# (export_study_cells.py codes both as decision_type='decline').
#
#   Rscript identify_disconnect_cells.R <raw_dir> <study_cells.csv> <out_csv>
#
# A (user, month, round) cell is a DISCONNECT cell iff, in raw UserAction:
#   (1) it holds a decision row (amount IS NULL) whose action is a Do-it verb:
#         "Follow"   -> forward   direction (Do it on a follow rec)
#         "Unfollow" -> step-back direction (Do it on an unfollow rec;
#                       "Never mind" is action="Decline" and never qualifies)
#   (2) the cell holds NO paired investment row — paired per chatbot/follow_ratio.py:
#         same portfolio (case-insensitive), direction-compliant verb
#         (follow -> Invest/Increase; unfollow -> Unfollow/Decrease),
#         amount IS NOT NULL
#       This is exactly export_study_cells._compliant_invested() returning None,
#       i.e. the cell landed in the exporter's `decline` branch. A cell with a
#       compliant amount==0 row is a zero_optout (the participant DID answer) and
#       is NOT a disconnect — see the sensitivity variants printed below.
#   (3) it passes the existing export filters: participant has a resolved
#       condition, username not TEST_USER*, decision row not adhoc-paired.
#
# Reads only. Writes only <out_csv>. Touches no frozen script.

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: Rscript identify_disconnect_cells.R <raw_dir> <study_cells.csv> <out_csv>")
raw_dir <- args[[1]]; cells_csv <- args[[2]]; out_csv <- args[[3]]

ADHOC_SOURCE <- "llm_adhoc"
DOIT_FORWARD   <- "Follow"     # Do it on a follow rec   (Never mind = "Decline")
DOIT_STEPBACK  <- "Unfollow"   # Do it on an unfollow rec (forward-Unfollow decision row)
FOLLOW_COMPLIANT   <- c("Invest", "Increase")
UNFOLLOW_COMPLIANT <- c("Unfollow", "Decrease")

ua   <- fread(file.path(raw_dir, "UserAction.csv"),     na.strings = "", encoding = "UTF-8",
              colClasses = list(character = c("recommendation", "client_nonce", "portfolio", "action")))
recs <- fread(file.path(raw_dir, "Recommendation.csv"), na.strings = "", encoding = "UTF-8")
part <- fread(file.path(raw_dir, "Participant.csv"),    na.strings = "", encoding = "UTF-8")
cells_export <- fread(cells_csv, na.strings = "", encoding = "UTF-8")
stopifnot(nrow(ua) > 0, nrow(recs) > 0, nrow(part) > 0, nrow(cells_export) > 0)
for (col in c("user", "month", "round", "action", "amount", "portfolio", "recommendation", "client_nonce"))
  if (!col %in% names(ua)) stop(sprintf("UserAction missing column: %s", col))

ua[, month := as.integer(month)]
ua[, round := as.integer(round)]
ua[, amount := as.numeric(amount)]
stopifnot(!anyNA(ua$month), !anyNA(ua$round), !anyNA(ua$user), !anyNA(ua$action))

## ---- 1. resolve the Recommendation FK by STABLE ID ------------------------
# The Django CSV export renders the FK through __str__ ("<user>: <direction> <profile>"),
# which is NOT unique (35/450 reprs collide here) — joining on it would misclassify.
# client_nonce carries the stable id as "rec-<id>[:suffix]" (execute_action writes it),
# so resolve on that and VERIFY the resolved row's __str__ reproduces the export column.
# A row that cannot be resolved this way stops the run rather than defaulting either way.
recs[, id := as.integer(id)]
recs[, str_repr := paste0(user, ": ", direction, " ", profile_name)]
adhoc_ids <- recs[source == ADHOC_SOURCE, id]

ua[, ua_key := .I]
ua[, rec_id := as.integer(sub("^rec-([0-9]+).*$", "\\1", client_nonce))]
ua[!grepl("^rec-[0-9]+", client_nonce), rec_id := NA_integer_]

# Free-path writes (portfolio UI rather than /execute_action/) set the FK but carry a
# UUID nonce, so the stable id is unavailable for them. Fall back to (__str__, month):
# a rec's month is in the export, and it breaks every collision present here. Resolve
# ONLY when exactly one candidate remains — 0 or 2+ stops the run, because guessing
# would silently move a row in or out of the adhoc exclusion.
need_fb <- ua[!is.na(recommendation) & is.na(rec_id),
              .(ua_key, month, str_repr = recommendation)]
if (nrow(need_fb) > 0) {
  cand <- merge(need_fb, recs[, .(str_repr, rmonth = as.integer(month), fb_id = id)],
                by = "str_repr", allow.cartesian = TRUE)[month == rmonth]
  hits <- cand[, .(n = .N, fb_id = fb_id[[1]]), by = ua_key]
  ambiguous <- hits[n != 1]
  missed    <- setdiff(need_fb$ua_key, hits$ua_key)
  if (nrow(ambiguous) > 0 || length(missed) > 0)
    stop(sprintf("FK fallback failed: %d row(s) ambiguous on (__str__, month), %d with no candidate",
                 nrow(ambiguous), length(missed)))
  ua[hits, rec_id := i.fb_id, on = "ua_key"]
  cat(sprintf("[fk] %d row(s) resolved by nonce, %d by (__str__, month) fallback, %d NULL-FK\n",
              ua[!is.na(recommendation) & !ua_key %in% need_fb$ua_key, .N],
              nrow(need_fb), ua[is.na(recommendation), .N]))
}
orphan <- setdiff(na.omit(ua$rec_id), recs$id)
if (length(orphan)) stop(sprintf("client_nonce references unknown Recommendation id(s): %s",
                                 paste(orphan, collapse = ", ")))
ua <- merge(ua, recs[, .(rec_id = id, rec_str = str_repr, rec_source = source,
                         rec_direction = direction)],
            by = "rec_id", all.x = TRUE, sort = FALSE)
bad_repr <- ua[!is.na(rec_id) & !is.na(recommendation) & rec_str != recommendation]
if (nrow(bad_repr) > 0)
  stop(sprintf("%d row(s) where the nonce-resolved rec does not match the exported FK string (e.g. %s vs %s)",
               nrow(bad_repr), bad_repr$rec_str[[1]], bad_repr$recommendation[[1]]))
# _NOT_ADHOC (follow_ratio.py): a NULL-FK row is KEPT; only an adhoc-paired row is dropped.
ua[, is_adhoc := !is.na(rec_id) & rec_id %in% adhoc_ids]

## ---- 2. existing population filter ----------------------------------------
# Mirrors export_study_cells._rows: Participant row required, condition resolved,
# username not TEST_USER*. (`user` is the username in this export.)
eligible <- part[!is.na(condition) & condition != "" & !startsWith(user, "TEST_USER"), unique(user)]
cat(sprintf("[filter] %d participants, %d eligible (resolved condition, non-TEST_USER)\n",
            uniqueN(part$user), length(eligible)))
dropped_users <- setdiff(unique(ua$user), eligible)
if (length(dropped_users))
  cat(sprintf("[filter] users in UserAction excluded from the study population: %s\n",
              paste(dropped_users, collapse = ", ")))

keep <- ua[user %in% eligible & !is_adhoc]

## ---- 3. cells, decision rows, and reconciliation against the export --------
decisions <- keep[is.na(amount)]
cells <- unique(decisions[, .(user, month, round)])
setorder(cells, user, month, round)
cat(sprintf("[cells] %d decision rows -> %d distinct cells (export CSV has %d)\n",
            nrow(decisions), nrow(cells), nrow(cells_export)))
# One decision row per cell is a write-layer invariant; more than one means the
# earliest governs (follow_ratio._decision_row). Report it rather than assume.
multi <- decisions[, .N, by = .(user, month, round)][N > 1]
if (nrow(multi) > 0) {
  cat(sprintf("[cells] WARNING: %d cell(s) hold >1 decision row; earliest-wins applied\n", nrow(multi)))
  print(multi)
  setorder(decisions, user, month, round, created_at)
  decisions <- decisions[, .SD[1], by = .(user, month, round)]
}

recon_export <- unique(cells_export[, .(user = participant_id,
                                        month = as.integer(month), round = as.integer(round))])
only_here   <- fsetdiff(cells, recon_export)
only_export <- fsetdiff(recon_export, cells)
recon_ok <- nrow(only_here) == 0 && nrow(only_export) == 0
if (recon_ok) {
  cat("[recon] OK: reconstructed cell set is identical to study_cells.csv\n")
} else {
  cat(sprintf("[recon] MISMATCH: %d cell(s) only in raw reconstruction, %d only in export\n",
              nrow(only_here), nrow(only_export)))
  if (nrow(only_here))   print(only_here)
  if (nrow(only_export)) print(only_export)
}

## ---- 4. paired-investment lookup ------------------------------------------
# Paired per follow_ratio.is_follow / export._compliant_invested: same cell, same
# portfolio (case-insensitive), direction-compliant verb, non-adhoc.
invest <- keep[!is.na(amount)]
invest[, pkey := tolower(portfolio)]
decisions[, pkey := tolower(portfolio)]

n_compliant <- function(drow, verbs, require_positive) {
  rows <- invest[user == drow$user & month == drow$month & round == drow$round &
                   pkey == drow$pkey & action %in% verbs]
  if (require_positive) rows <- rows[amount > 0]
  nrow(rows)
}
n_any_invest <- function(drow)
  invest[user == drow$user & month == drow$month & round == drow$round, .N]

doit <- decisions[action %in% c(DOIT_FORWARD, DOIT_STEPBACK)]
cat(sprintf("[decision verbs] %s\n",
            paste(sprintf("%s=%d", decisions[, .N, by = action]$action,
                          decisions[, .N, by = action]$N), collapse = "  ")))

res <- rbindlist(lapply(seq_len(nrow(doit)), function(i) {
  d <- doit[i]
  verbs <- if (d$action == DOIT_FORWARD) FOLLOW_COMPLIANT else UNFOLLOW_COMPLIANT
  data.table(
    participant_id = d$user, month = d$month, round = d$round,
    decision_action = d$action,
    disconnect_type = if (d$action == DOIT_FORWARD) "follow_doit" else "stepback_doit",
    portfolio = d$portfolio,
    rec_direction = if (is.na(d$rec_direction)) "" else d$rec_direction,
    decision_created_at = d$created_at,
    n_paired_any_amount = n_compliant(d, verbs, FALSE),   # primary rule uses this
    n_paired_positive   = n_compliant(d, verbs, TRUE),
    n_any_invest_row    = n_any_invest(d)
  )
}), fill = TRUE)

# Primary rule (spec ②): no paired investment row at all -> exporter coded `decline`
primary   <- res[n_paired_any_amount == 0]
# Sensitivity A (wider): no paired row with amount>0 — also sweeps in zero_optout cells
variant_a <- res[n_paired_positive == 0]
# Sensitivity B (narrowest literal): no investment row of ANY kind in the cell
variant_b <- res[n_any_invest_row == 0]

setorder(primary, participant_id, month, round)
fwrite(primary[, .(participant_id, month, round, decision_action, disconnect_type,
                   portfolio, rec_direction, decision_created_at,
                   n_paired_any_amount, n_paired_positive, n_any_invest_row)], out_csv)

## ---- 5. report -------------------------------------------------------------
cat("\n=== disconnect-cell identification ===\n")
cat(sprintf("Do-it decision rows examined : %d  (Follow=%d, forward-Unfollow=%d)\n",
            nrow(res), res[decision_action == DOIT_FORWARD, .N], res[decision_action == DOIT_STEPBACK, .N]))
cat(sprintf("PRIMARY  disconnect cells    : %d   (Follow-type=%d, step-back-type=%d)\n",
            nrow(primary),
            primary[disconnect_type == "follow_doit", .N],
            primary[disconnect_type == "stepback_doit", .N]))
cat(sprintf("  sensitivity A (no paired amount>0, includes zero_optout) : %d\n", nrow(variant_a)))
cat(sprintf("  sensitivity B (no investment row of any kind in the cell): %d\n", nrow(variant_b)))
if (nrow(primary) > 0) print(primary)
cat(sprintf("wrote -> %s  (%d rows)\n", out_csv, nrow(primary)))
if (!recon_ok) {
  cat("\n[recon] FAILED — the reconstructed population does not match the frozen export.\n")
  quit(status = 1)
}
