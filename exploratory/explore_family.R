#!/usr/bin/env Rscript
# =============================================================================
# explore_family.R -- AsPredicted #301918, Q5 exploratory family (canonical).
#
# Registration Q5 (verbatim):
#   "All other analyses (CAS descriptives, condition x portfolio-complexity
#    interaction, reasoning open-rate/dwell dose-response, mechanism battery
#    contrasts, MC1/MC2, single-mediator descriptive estimates) are
#    exploratory, carry no inferential weight, and are reported with 90% CIs;
#    no multiplicity correction is applied because no exploratory result will
#    be described with significance language."
#
# Task -> registered phrase map:
#   T1 CAS descriptives ................ Q5 "CAS descriptives"
#   T2 complexity interaction .......... Q5 "condition x portfolio-complexity interaction"
#   T3 open-rate / dwell ............... Q5 "reasoning open-rate/dwell dose-response"
#   T4 mechanism battery ............... Q5 "mechanism battery contrasts"
#   T5 MC1 / MC2 ....................... Q5 "MC1/MC2"
#   T6 single-mediator descriptives .... Q5 "single-mediator descriptive estimates"
#   T7 ITT vs per-protocol split ....... Q8 "Reasoning open-rate defines a mandatory
#                                          ITT-vs-per-protocol split within the reasoning arm"
#   T8 do-it-then-disconnect ........... Q6 "A 'Do it'-then-disconnect cell is treated as
#                                          missing; sensitivity analysis re-codes it as 0"
#   T9 absolute action size ............ Q3 "Secondary (descriptive only): Concordant
#                                          Allocation Share (...), and absolute action size"
#
# NO INFERENTIAL WEIGHT. Every interval below is a 90% interval reported as an
# estimate of magnitude and uncertainty. Nothing here is a test; no result is
# described with significance language; no multiplicity correction is applied.
#
# Analysis sets:
#   behavioural  N = 36 participants / 387 cells (merged_cells_screened.csv)
#   questionnaire N = 35 (UCLCEUS0U3T has no post-task submission)
#
# Interval machinery:
#   cell-level statistics (follow rates, mean CAS)  -> participant cluster
#       bootstrap, percentile, B = 2000, seed = 20260727, STRATIFIED BY ARM
#       (arm allocation is fixed by the between-subjects design, so resampling
#       within arm preserves it; the archived confirmatory bootstrap
#       posthoc-复算归档-20260729/cluster_bootstrap_canonical.R is unstratified
#       because it resamples a single pooled Firth fit -- see provenance.md).
#   participant-level single proportions -> Wilson 90% interval (participants
#       are independent units; no clustering to absorb).
#   medians / IQRs / dwell -> point descriptives, no interval.
#
# Usage:
#   Rscript explore_family.R [DATA_DIR] [FROZEN_DIR] [OUT_DIR]
# =============================================================================

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR   <- if (length(args) >= 1) args[[1]] else
  "/Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道/data/real-20260727"
FROZEN_DIR <- if (length(args) >= 2) args[[2]] else
  "/Users/wolpfor/Desktop/UCL毕业项目/SocialTradingChatbot/chatbot/frozen_data"
OUT_DIR    <- if (length(args) >= 3) args[[3]] else "/tmp/exploratory-family"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BOOT_B    <- 2000L
BOOT_SEED <- 20260727L
CONF      <- 0.90

RES <- list()   # collector: every reported number lands here
add <- function(section, label, ...) {
  RES[[length(RES) + 1L]] <<- c(list(section = section, label = label), list(...))
  invisible(NULL)
}
fmt <- function(x, d = 3) if (is.na(x)) "NA" else formatC(x, format = "f", digits = d)

# ---- interval helpers -------------------------------------------------------

# Stratified-by-arm participant cluster bootstrap over a participant-level
# sufficient-statistic table. Every cell-level statistic used below is a ratio
# of sums over participants, so resampling participant rows reproduces the
# cell-level cluster bootstrap exactly and cheaply.
boot_p <- function(psum, statfun, B = BOOT_B, seed = BOOT_SEED) {
  set.seed(seed)
  idx_by_arm <- split(seq_len(nrow(psum)), psum$condition)
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    take <- unlist(lapply(idx_by_arm, function(ix) sample(ix, length(ix), replace = TRUE)),
                   use.names = FALSE)
    out[b] <- tryCatch(statfun(psum[take]), error = function(e) NA_real_)
  }
  ok <- out[is.finite(out)]
  if (!length(ok)) return(c(lo = NA_real_, hi = NA_real_, n_eff = 0))
  c(lo = unname(quantile(ok, (1 - CONF) / 2)),
    hi = unname(quantile(ok, 1 - (1 - CONF) / 2)),
    n_eff = length(ok))
}

wilson <- function(x, n, conf = CONF) {
  if (n == 0) return(c(lo = NA_real_, hi = NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n; d <- 1 + z^2 / n
  ctr <- p + z^2 / (2 * n); hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lo = (ctr - hw) / d, hi = (ctr + hw) / d)
}

cronbach <- function(m) {                       # m: participants x items matrix
  m <- m[complete.cases(m), , drop = FALSE]
  k <- ncol(m); if (k < 3 || nrow(m) < 3) return(NA_real_)
  (k / (k - 1)) * (1 - sum(apply(m, 2, var)) / var(rowSums(m)))
}

# =============================================================================
# 0. LOAD
# =============================================================================
cells <- fread(file.path(DATA_DIR, "merged_cells_screened.csv"), na.strings = "",
               encoding = "UTF-8")
quest <- fread(file.path(DATA_DIR, "questionnaire_items.csv"), na.strings = "",
               encoding = "UTF-8", colClasses = list(character = "created_at"))
recs  <- fread(file.path(DATA_DIR, "raw", "Recommendation.csv"), na.strings = "",
               encoding = "UTF-8")
revt  <- fread(file.path(DATA_DIR, "raw", "ReasoningEvent.csv"), na.strings = "",
               encoding = "UTF-8")
fw    <- fread(file.path(DATA_DIR, "raw", "FactorWeight.csv"), na.strings = "",
               encoding = "UTF-8")
disc  <- fread(file.path(DATA_DIR, "disconnect_cells.csv"), na.strings = "",
               encoding = "UTF-8")

cells[, followed := as.integer(followed)]
cells[, month := as.integer(month)][, round := as.integer(round)]
cells[, cas := as.numeric(cas)][, invested := as.numeric(invested)]
PIDS <- sort(unique(cells$participant_id))
stopifnot(nrow(cells) == 387L, length(PIDS) == 36L)

arm_of <- unique(cells[, .(participant_id, condition)])
stopifnot(nrow(arm_of) == length(PIDS))

# =============================================================================
# T2-prep. PORTFOLIO COMPLEXITY TIER  (registered: "condition x portfolio-complexity")
#   Classification is NOT invented here. It is the `complexity_tier` column of
#   the frozen stimulus tables in the deployed repository:
#     chatbot/frozen_data/portfolio-compositions.csv      (traders 1-10, v1)
#     chatbot/frozen_data/portfolio-compositions-v2.csv   (traders 11-20, v2)
#   The deployed 10-trader set draws 4 from v1 and 6 from v2; see provenance.md.
#   Integrity gate: the live FactorWeight rows must reproduce the frozen
#   composition percentages exactly, else stop.
# =============================================================================
comp <- rbindlist(list(
  fread(file.path(FROZEN_DIR, "portfolio-compositions.csv"), na.strings = "", encoding = "UTF-8"),
  fread(file.path(FROZEN_DIR, "portfolio-compositions-v2.csv"), na.strings = "", encoding = "UTF-8")
))
FACTORS <- c("energy", "technology", "financials", "consumer_discretionary", "consumer_staples")
comp_long <- melt(comp, id.vars = c("portfolio", "complexity_tier"),
                  measure.vars = paste0(FACTORS, "_pct"),
                  variable.name = "factor", value.name = "pct")
comp_long[, factor := sub("_pct$", "", as.character(factor))]
comp_long <- comp_long[pct > 0]

deployed <- sort(unique(fw$profile))
tier_map <- unique(comp[portfolio %in% deployed, .(portfolio, complexity_tier)])
stopifnot(!anyDuplicated(tier_map$portfolio),
          setequal(tier_map$portfolio, deployed))

# integrity gate: live weights == frozen composition, trader by trader
chk <- merge(fw[, .(portfolio = profile, factor, weight = as.numeric(weight))],
             comp_long[portfolio %in% deployed, .(portfolio, factor, pct = as.numeric(pct))],
             by = c("portfolio", "factor"), all = TRUE)
if (nrow(chk[is.na(weight) | is.na(pct) | abs(weight - pct) > 1e-9]) > 0)
  stop("live FactorWeight does not match frozen portfolio-compositions: complexity tier unsafe")
cat(sprintf("[gate] complexity tier verified against live FactorWeight for %d deployed traders\n",
            length(deployed)))

# cell -> recommended trader -> tier  (scheduled cards only; ad-hoc excluded by design)
recs[, round_i := suppressWarnings(as.integer(sub("^.*R(\\d+)$", "\\1", round)))]
sched <- recs[source == "scheduled" & !is.na(round_i),
              .(participant_id = user, month = as.integer(month), round = round_i,
                profile_name, rec_direction = direction)]
stopifnot(!anyDuplicated(sched[, .(participant_id, month, round)]))
cells <- merge(cells, sched, by = c("participant_id", "month", "round"),
               all.x = TRUE, sort = FALSE)
if (anyNA(cells$profile_name)) stop("cells without a scheduled Recommendation row")
if (nrow(cells[direction != rec_direction]) > 0)
  stop("direction disagreement between cells and Recommendation")
cells <- merge(cells, tier_map, by.x = "profile_name", by.y = "portfolio",
               all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(cells$complexity_tier))
# binary collapse used by the manuscript's exploratory contrast (offsetting vs rest)
cells[, offsetting := as.integer(complexity_tier == "two-factor-offsetting")]

# =============================================================================
# T3-prep. REASONING PANEL TELEMETRY
#   ReasoningEvent event_type vocabulary (observed, whole table):
#     flash_shown        -- panel auto-flashed open on card delivery (1 per block)
#     flash_autocollapse -- the same auto-flash closing itself
#     open               -- participant-initiated open
#     collapse           -- close of a participant-initiated open; carries dwell_ms
#   block_id = one reasoning panel = one recommendation card. surface is always
#   'main'; visible always True; condition_snapshot always llm-reasoning.
# =============================================================================
revt[, month := as.integer(month)][, dwell_ms := as.numeric(dwell_ms)]
revt <- revt[user %in% PIDS]                       # analysis set only
REASONERS <- arm_of[condition == "llm-reasoning", participant_id]
stopifnot(all(revt$user %in% REASONERS))

blocks <- revt[event_type == "flash_shown", .(n_flash = .N), by = .(user, block_id, month)]
opens  <- revt[event_type == "open",         .(n_open  = .N), by = .(user, block_id, month)]
blocks <- merge(blocks, opens, by = c("user", "block_id", "month"), all = TRUE)
blocks[is.na(n_open), n_open := 0L][is.na(n_flash), n_flash := 0L]
blocks[, opened := as.integer(n_open > 0)]
dwell <- revt[event_type == "collapse" & is.finite(dwell_ms), .(user, month, block_id, dwell_ms)]

open_p <- data.table(participant_id = REASONERS)
open_p <- merge(open_p, blocks[, .(n_blocks = .N, n_blocks_opened = sum(opened)),
                               by = .(participant_id = user)],
                by = "participant_id", all.x = TRUE)
open_p[is.na(n_blocks), n_blocks := 0L][is.na(n_blocks_opened), n_blocks_opened := 0L]
open_p[, open_rate := ifelse(n_blocks > 0, n_blocks_opened / n_blocks, NA_real_)]
open_p[, ever_opened := as.integer(n_blocks_opened > 0)]

# (participant, month) open indicator, for the month-level dose reading
pm_open <- blocks[, .(month_opened = as.integer(sum(opened) > 0)), by = .(participant_id = user, month)]

# =============================================================================
# T4-prep. QUESTIONNAIRE SCALES  (main.tex:1886 four-block on-mechanism battery)
#   Block SP  items 1-5  self-authored social-presence (descriptive, outside mechanism)
#   Block T   B1-B3      Hoffman-2023 Explanation Satisfaction  [reasoning arm ONLY]
#   Block U   C1-C3      Madsen-Gregor understandability        [both arms]
#   Block R   C4-C5      Madsen-Gregor reliability proxy        [both arms]
#   Block TR  C6-C7      Korber-2019 trust/competence proxy     [both arms]
#   MC1       H1         perceived-reasoning awareness          [both arms]
#   MC2       H2         named-reason recall                    [reasoning arm ONLY]
# =============================================================================
post <- quest[phase == "posttask" & participant_id %in% PIDS]
# De-dup: same registered rule prepare_cells.R applies to the pre-task instrument --
# a participant with more than one submission keeps the EARLIEST. One participant
# (UCLGVXX8SZQ) double-submitted the post-task form 0.68 s apart. created_at is
# ISO-8601 on one fixed UTC offset, so lexicographic min == chronological min;
# both properties are asserted before being relied on.
post[, created_at := as.character(created_at)]
stopifnot(all(grepl("^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}", post$created_at)))
if (uniqueN(sub("^.*(\\+\\d{2}:\\d{2}|Z)$", "\\1", post$created_at)) > 1)
  stop("mixed created_at UTC offsets: lexicographic de-dup would be wrong")
n_pre_dedup <- nrow(post)
post <- post[post[, .(created_at = min(created_at)), by = participant_id],
             on = .(participant_id, created_at)]
if (n_pre_dedup != nrow(post))
  cat(sprintf("[dedup] post-task: dropped %d duplicate-submission rows (%d -> %d)\n",
              n_pre_dedup - nrow(post), n_pre_dedup, nrow(post)))
stopifnot(post[, uniqueN(created_at), by = participant_id][, all(V1 == 1L)])
item_of <- function(prefix) {
  x <- post[startsWith(question, prefix), .(participant_id, answer)]
  if (anyDuplicated(x$participant_id)) stop(sprintf("duplicate item rows for %s", prefix))
  setnames(x, "answer", prefix)[]
}
SP_ITEMS <- c("1. The chatbot", "2. The chatbot", "3. The chatbot",
              "4. The chatbot", "5. The chatbot")
LIK <- list(
  SP = SP_ITEMS,
  T  = c("B1.", "B2.", "B3."),
  U  = c("C1.", "C2.", "C3."),
  R  = c("C4.", "C5."),
  TR = c("C6.", "C7."),
  MC1 = c("H1.")
)
qwide <- data.table(participant_id = sort(unique(post$participant_id)))
for (p in unique(unlist(LIK))) qwide <- merge(qwide, item_of(p), by = "participant_id", all.x = TRUE)
for (p in unique(unlist(LIK))) set(qwide, j = p, value = suppressWarnings(as.numeric(qwide[[p]])))
for (s in names(LIK)) {
  cols <- LIK[[s]]
  qwide[, (s) := rowMeans(as.matrix(.SD), na.rm = FALSE), .SDcols = cols]
}
mc2 <- item_of("H2. ")
setnames(mc2, "H2. ", "mc2_raw")
MC2_CORRECT <- "a central bank interest-rate rise"   # the only option naming a real
                                                     # MarketEvent (M1); see provenance.md
mc2[, mc2_correct := as.integer(mc2_raw == MC2_CORRECT)]
qwide <- merge(qwide, mc2, by = "participant_id", all.x = TRUE)
qwide <- merge(qwide, arm_of, by = "participant_id", all.x = TRUE)
N_Q <- nrow(qwide)

# =============================================================================
# PARTICIPANT-LEVEL SUFFICIENT STATISTICS (the bootstrap unit)
# =============================================================================
psum <- cells[, .(n_cells = .N, n_follow = sum(followed),
                  cas_sum = sum(cas[decision_type == "follow"], na.rm = TRUE),
                  cas_n   = sum(decision_type == "follow" & !is.na(cas)),
                  # T9 absolute action size: `invested` is the GBP size of the capital
                  # move on the cell. Two denominators, see the T9 block below.
                  inv_sum_move = sum(invested[!is.na(invested) & invested > 0]),
                  inv_n_move   = sum(!is.na(invested) & invested > 0),
                  inv_sum_fol  = sum(invested[decision_type == "follow"]),
                  inv_n_fol    = sum(decision_type == "follow"),
                  n_off   = sum(offsetting), n_off_follow = sum(followed * offsetting),
                  n_non   = sum(1L - offsetting), n_non_follow = sum(followed * (1L - offsetting)),
                  n_sf    = sum(complexity_tier == "single-factor"),
                  n_sf_follow = sum(followed * (complexity_tier == "single-factor")),
                  n_2fo   = sum(complexity_tier == "two-factor-offsetting"),
                  n_2fo_follow = sum(followed * (complexity_tier == "two-factor-offsetting")),
                  n_3fm   = sum(complexity_tier == "three-factor-mixed"),
                  n_3fm_follow = sum(followed * (complexity_tier == "three-factor-mixed"))),
              by = .(participant_id, condition)]
psum <- merge(psum, open_p[, .(participant_id, n_blocks, n_blocks_opened, open_rate, ever_opened)],
              by = "participant_id", all.x = TRUE)
psum <- merge(psum, qwide[, .(participant_id, SP, T, U, R, TR, MC1, mc2_correct)],
              by = "participant_id", all.x = TRUE)
psum[, follow_ratio := n_follow / n_cells]
setorder(psum, condition, participant_id)
stopifnot(nrow(psum) == 36L)

rate  <- function(d, num, den) sum(d[[num]]) / sum(d[[den]])
arm   <- function(d, a) d[condition == a]
d_rate <- function(d, num, den)
  rate(arm(d, "llm-reasoning"), num, den) - rate(arm(d, "llm-control"), num, den)
d_mean <- function(d, col) mean(arm(d, "llm-reasoning")[[col]], na.rm = TRUE) -
                           mean(arm(d, "llm-control")[[col]],  na.rm = TRUE)

N_R <- psum[condition == "llm-reasoning", .N]; N_C <- psum[condition == "llm-control", .N]
cat(sprintf("[set] behavioural N=%d (reasoning %d / control %d), cells=%d; questionnaire N=%d\n",
            nrow(psum), N_R, N_C, nrow(cells), N_Q))

# =============================================================================
# T1. CAS DESCRIPTIVES              [Q5 "CAS descriptives"]   N=36 behavioural
#   CAS is defined ONLY on cells with decision_type == 'follow' (a
#   conditional-on-follow share of the available balance), matching
#   fit_glmm.R's口径 and registration Q3.
# =============================================================================
cas_cells <- cells[decision_type == "follow" & !is.na(cas)]
for (a in c("llm-reasoning", "llm-control")) {
  v <- cas_cells[condition == a, cas]
  add("T1 CAS", sprintf("CAS %s", a), n_cells = length(v),
      n_participants = uniqueN(cas_cells[condition == a, participant_id]),
      mean = mean(v), median = median(v),
      q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)))
}
ci <- boot_p(psum, function(d) {
  r <- arm(d, "llm-reasoning"); c0 <- arm(d, "llm-control")
  sum(r$cas_sum) / sum(r$cas_n) - sum(c0$cas_sum) / sum(c0$cas_n)
})
add("T1 CAS", "CAS mean difference (reasoning - control)",
    est = sum(arm(psum, "llm-reasoning")$cas_sum) / sum(arm(psum, "llm-reasoning")$cas_n) -
          sum(arm(psum, "llm-control")$cas_sum)  / sum(arm(psum, "llm-control")$cas_n),
    lo = ci[["lo"]], hi = ci[["hi"]], method = "participant cluster bootstrap, stratified by arm")

# =============================================================================
# T9. ABSOLUTE ACTION SIZE   [Q3 "Secondary (descriptive only): Concordant
#     Allocation Share (...), and absolute action size"]        N=36 behavioural
#
#   Carrier: `invested` in merged_cells_screened.csv -- the GBP size of the
#   capital move executed on that cell (see provenance.md §5). Unlike CAS it is
#   NOT restricted to follow-direction rounds, so two denominators are reported:
#     PRIMARY   every cell carrying a positive capital move (follow adoptions +
#               step-back/divest adoptions). This is what "absolute" adds over
#               CAS: an unsigned magnitude that CAS structurally cannot hold,
#               because a share of the available balance is undefined for a
#               withdrawal (apply_screening.R, structural-NA rule).
#     SECONDARY follow cells only -- the same subset CAS lives on, so the two
#               secondaries can be read side by side on one denominator.
#   The 4 zero_optout cells (accepted, allocated GBP 0) sit in neither: they
#   carry no capital move. Adding them as zeros would move the primary mean by
#   at most ~1.4% of cells; they are non-adherence under the DV coding.
#   Declines carry invested = NA (structural, no move made).
# =============================================================================
mv_cells  <- cells[!is.na(invested) & invested > 0]
fol_cells <- cells[decision_type == "follow"]
for (g in list(list("PRIMARY: all cells with a positive capital move", mv_cells),
               list("SECONDARY: follow cells only (CAS denominator)", fol_cells))) {
  for (a in c("llm-reasoning", "llm-control")) {
    v <- g[[2]][condition == a, invested]
    add("T9 action size", sprintf("action size GBP | %s | %s", g[[1]], a),
        n_cells = length(v), n_participants = uniqueN(g[[2]][condition == a, participant_id]),
        mean = mean(v), median = median(v),
        q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)))
  }
}
for (g in list(c("PRIMARY: all cells with a positive capital move", "inv_sum_move", "inv_n_move"),
               c("SECONDARY: follow cells only (CAS denominator)", "inv_sum_fol", "inv_n_fol"))) {
  st <- function(d) rate(arm(d, "llm-reasoning"), g[2], g[3]) - rate(arm(d, "llm-control"), g[2], g[3])
  ci <- boot_p(psum, st)
  add("T9 action size", sprintf("action size GBP | mean difference (reasoning - control) | %s", g[1]),
      est = st(psum), lo = ci[["lo"]], hi = ci[["hi"]],
      method = "participant cluster bootstrap, stratified by arm, same seed/B as T1 CAS")
}
add("T9 action size", "cells carrying no capital move (excluded from both denominators)",
    n_decline_NA_invested = cells[decision_type == "decline", .N],
    n_zero_optout = cells[decision_type == "zero_optout", .N],
    method = "decline rows carry invested = NA (structural, no move); zero_optout rows carry invested = 0")

# =============================================================================
# T2. CONDITION x PORTFOLIO COMPLEXITY   [Q5 "condition x portfolio-complexity interaction"]
# =============================================================================
tier_counts <- cells[, .(n_cells = .N, n_traders = uniqueN(profile_name)), by = complexity_tier]
for (i in seq_len(nrow(tier_counts)))
  add("T2 complexity", sprintf("stimulus coverage: %s", tier_counts$complexity_tier[i]),
      n_cells = tier_counts$n_cells[i], n_traders = tier_counts$n_traders[i])

for (tg in list(c("single-factor", "n_sf", "n_sf_follow"),
                c("two-factor-offsetting", "n_2fo", "n_2fo_follow"),
                c("three-factor-mixed", "n_3fm", "n_3fm_follow"))) {
  for (a in c("llm-reasoning", "llm-control")) {
    d <- arm(psum, a)
    add("T2 complexity", sprintf("follow rate | %s | %s", tg[1], a),
        n_cells = sum(d[[tg[2]]]), n_follow = sum(d[[tg[3]]]),
        est = sum(d[[tg[3]]]) / sum(d[[tg[2]]]))
  }
  ci <- boot_p(psum, function(d) d_rate(d, tg[3], tg[2]))
  add("T2 complexity", sprintf("arm difference | %s", tg[1]),
      est = d_rate(psum, tg[3], tg[2]), lo = ci[["lo"]], hi = ci[["hi"]],
      method = "cluster bootstrap")
}
for (bg in list(c("offsetting", "n_off", "n_off_follow"),
                c("non-offsetting (single-factor + three-factor-mixed)", "n_non", "n_non_follow"))) {
  for (a in c("llm-reasoning", "llm-control")) {
    d <- arm(psum, a)
    add("T2 complexity", sprintf("follow rate | %s | %s", bg[1], a),
        n_cells = sum(d[[bg[2]]]), n_follow = sum(d[[bg[3]]]),
        est = sum(d[[bg[3]]]) / sum(d[[bg[2]]]))
  }
  ci <- boot_p(psum, function(d) d_rate(d, bg[3], bg[2]))
  add("T2 complexity", sprintf("arm difference | %s", bg[1]),
      est = d_rate(psum, bg[3], bg[2]), lo = ci[["lo"]], hi = ci[["hi"]],
      method = "cluster bootstrap")
}
did <- function(d) d_rate(d, "n_off_follow", "n_off") - d_rate(d, "n_non_follow", "n_non")
ci <- boot_p(psum, did)
add("T2 complexity", "INTERACTION contrast: (arm diff on offsetting) - (arm diff on non-offsetting)",
    est = did(psum), lo = ci[["lo"]], hi = ci[["hi"]],
    method = "difference-in-differences of follow rates, cluster bootstrap")

# GLMM with the interaction term, attempted as a secondary read. The confirmatory
# analysis already fell back off the GLMM at this N, so convergence status is
# reported rather than assumed.
glmm_note <- tryCatch({
  suppressMessages(library(lme4))
  cells[, cond_f := relevel(factor(condition), ref = "llm-control")]
  cells[, off_f  := factor(offsetting, levels = c(0, 1), labels = c("non-offsetting", "offsetting"))]
  m <- glmer(followed ~ cond_f * off_f + literacy + is_HCI + (1 | participant_id),
             data = cells, family = binomial,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  s <- summary(m)
  msgs <- m@optinfo$conv$lme4$messages
  cf <- s$coefficients
  term <- grep(":", rownames(cf), value = TRUE)[1]
  b <- cf[term, "Estimate"]; se <- cf[term, "Std. Error"]
  z <- qnorm(1 - (1 - CONF) / 2)
  add("T2 complexity", "GLMM interaction OR (secondary read)",
      term = term, est_OR = exp(b), lo = exp(b - z * se), hi = exp(b + z * se),
      method = "binomial GLMM, Wald 90% interval",
      convergence = if (is.null(msgs)) "no lme4 convergence message" else paste(msgs, collapse = "; "))
  if (is.null(msgs)) "converged (no message)" else paste(msgs, collapse = "; ")
}, error = function(e) {
  add("T2 complexity", "GLMM interaction OR (secondary read)", est = NA_real_,
      method = "binomial GLMM", convergence = paste("ERROR:", conditionMessage(e)))
  paste("ERROR:", conditionMessage(e))
})
cat(sprintf("[T2] GLMM interaction status: %s\n", glmm_note))

# =============================================================================
# T3. OPEN-RATE / DWELL DOSE-RESPONSE  [Q5 "reasoning open-rate/dwell dose-response"]
#   Reasoning arm only. Card-level unit = one reasoning block_id.
# =============================================================================
n_blk <- nrow(blocks); n_blk_open <- sum(blocks$opened)
add("T3 open/dwell", "card-level open rate (panels opened / panels rendered)",
    n_panels = n_blk, n_opened = n_blk_open, est = n_blk_open / n_blk,
    method = paste("DENOMINATOR = reasoning panel RENDERS, not unique scheduled cards.",
                   "chat.js:496 emits one flash_shown per reasoning-bearing bot message,",
                   "so 273 panels sit above the 226 scheduled recommendations these 18",
                   "participants received (follow-ups and ad-hoc cards also carry a panel).",
                   "No interval: panels are clustered in participants; the participant-level",
                   "figures below are the units that carry uncertainty."))
add("T3 open/dwell", "participants with zero opens (reasoning arm)",
    n = open_p[ever_opened == 0, .N], of = nrow(open_p),
    est = open_p[ever_opened == 0, .N] / nrow(open_p),
    lo = wilson(open_p[ever_opened == 0, .N], nrow(open_p))[["lo"]],
    hi = wilson(open_p[ever_opened == 0, .N], nrow(open_p))[["hi"]],
    method = "Wilson 90% interval")
add("T3 open/dwell", "per-participant open rate (reasoning arm)",
    n = nrow(open_p), mean = mean(open_p$open_rate, na.rm = TRUE),
    median = median(open_p$open_rate, na.rm = TRUE),
    q1 = unname(quantile(open_p$open_rate, .25, na.rm = TRUE)),
    q3 = unname(quantile(open_p$open_rate, .75, na.rm = TRUE)))
add("T3 open/dwell", "blocks flashed per participant",
    n = nrow(open_p), median = median(open_p$n_blocks),
    q1 = unname(quantile(open_p$n_blocks, .25)), q3 = unname(quantile(open_p$n_blocks, .75)))
add("T3 open/dwell", "dwell_ms per open-collapse pair",
    n_pairs = nrow(dwell), median = median(dwell$dwell_ms),
    q1 = unname(quantile(dwell$dwell_ms, .25)), q3 = unname(quantile(dwell$dwell_ms, .75)),
    min = min(dwell$dwell_ms), max = max(dwell$dwell_ms),
    method = "interval telemetry cannot separate a glance from a read (main.tex:1884); read conservatively")

# month-level dose: cells in a (participant, month) with >=1 open vs 0 opens
cells_r <- merge(cells[condition == "llm-reasoning"], pm_open,
                 by = c("participant_id", "month"), all.x = TRUE)
cells_r[is.na(month_opened), month_opened := 0L]
pm_sum <- cells_r[, .(n_op = sum(month_opened), n_op_follow = sum(followed * month_opened),
                      n_no = sum(1L - month_opened), n_no_follow = sum(followed * (1L - month_opened))),
                  by = .(participant_id, condition)]
for (g in list(c("month with >=1 panel open", "n_op", "n_op_follow"),
               c("month with 0 panel opens", "n_no", "n_no_follow")))
  add("T3 open/dwell", sprintf("follow rate | %s (reasoning arm)", g[1]),
      n_cells = sum(pm_sum[[g[2]]]), n_follow = sum(pm_sum[[g[3]]]),
      est = sum(pm_sum[[g[3]]]) / sum(pm_sum[[g[2]]]))
ci <- boot_p(cbind(pm_sum, dummy = 0), function(d)
  sum(d$n_op_follow) / sum(d$n_op) - sum(d$n_no_follow) / sum(d$n_no))
add("T3 open/dwell", "within-arm month-level dose contrast (open month - no-open month)",
    est = sum(pm_sum$n_op_follow) / sum(pm_sum$n_op) - sum(pm_sum$n_no_follow) / sum(pm_sum$n_no),
    lo = ci[["lo"]], hi = ci[["hi"]],
    method = "cluster bootstrap within reasoning arm; within-participant self-selected exposure, not causal")

# participant-level open-rate -> follow-ratio association (reasoning arm)
pr <- psum[condition == "llm-reasoning" & !is.na(open_rate)]
r_est <- cor(pr$open_rate, pr$follow_ratio)
set.seed(BOOT_SEED)
rb <- replicate(BOOT_B, { i <- sample.int(nrow(pr), nrow(pr), TRUE)
  suppressWarnings(cor(pr$open_rate[i], pr$follow_ratio[i])) })
rb <- rb[is.finite(rb)]
add("T3 open/dwell", "participant-level correlation: open rate vs follow ratio (reasoning arm)",
    n = nrow(pr), est_r = r_est,
    lo = unname(quantile(rb, .05)), hi = unname(quantile(rb, .95)),
    method = "Pearson r, participant bootstrap 90% percentile interval")

# =============================================================================
# T4. MECHANISM BATTERY CONTRASTS   [Q5 "mechanism battery contrasts"]  N=35
#   Only blocks rendered to BOTH arms can be contrasted: SP, U, R, TR.
#   Block T (B1-B3) is reasoning-arm-only -> single-arm description.
# =============================================================================
SCALES <- c(U = "U", R = "R", TR = "TR", SP = "SP")
SCALE_LABEL <- c(U  = "Block U understandability (C1-C3, Madsen-Gregor)",
                 R  = "Block R reliability proxy (C4-C5, Madsen-Gregor)",
                 TR = "Block TR trust/competence proxy (C6-C7, Korber)",
                 SP = "Social presence (items 1-5, self-authored, unvalidated)")
for (s in names(SCALES)) {
  for (a in c("llm-reasoning", "llm-control")) {
    v <- qwide[condition == a, get(s)]; v <- v[!is.na(v)]
    add("T4 mechanism", sprintf("%s | %s", SCALE_LABEL[[s]], a),
        n = length(v), mean = mean(v), sd = sd(v), median = median(v),
        q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)))
  }
  ci <- boot_p(psum[!is.na(get(s))], function(d) d_mean(d, s))
  add("T4 mechanism", sprintf("%s | arm difference (reasoning - control)", SCALE_LABEL[[s]]),
      est = d_mean(psum, s), lo = ci[["lo"]], hi = ci[["hi"]],
      method = "participant bootstrap, stratified by arm")
  m <- as.matrix(qwide[, LIK[[s]], with = FALSE])
  a_all <- cronbach(m)
  if (!is.na(a_all))
    add("T4 mechanism", sprintf("%s | Cronbach alpha (both arms pooled)", SCALE_LABEL[[s]]),
        n = sum(complete.cases(m)), est = a_all,
        method = "k>=3 items only; 2-item blocks report inter-item correlation instead")
  else {
    mm <- m[complete.cases(m), , drop = FALSE]
    add("T4 mechanism", sprintf("%s | inter-item correlation (2-item block)", SCALE_LABEL[[s]]),
        n = nrow(mm), est = cor(mm[, 1], mm[, 2]))
  }
}
vT <- qwide[condition == "llm-reasoning" & !is.na(T), T]
add("T4 mechanism", "Block T explanation satisfaction (B1-B3, Hoffman) | reasoning arm ONLY",
    n = length(vT), mean = mean(vT), sd = sd(vT), median = median(vT),
    q1 = unname(quantile(vT, .25)), q3 = unname(quantile(vT, .75)),
    method = "single-arm description; control arm saw no explanation to rate, no contrast exists")
mT <- as.matrix(qwide[condition == "llm-reasoning", LIK$T, with = FALSE])
add("T4 mechanism", "Block T | Cronbach alpha (reasoning arm)",
    n = sum(complete.cases(mT)), est = cronbach(mT))

# =============================================================================
# T5. MC1 / MC2                     [Q5 "MC1/MC2"]
# =============================================================================
for (a in c("llm-reasoning", "llm-control")) {
  v <- qwide[condition == a & !is.na(MC1), MC1]
  add("T5 MC", sprintf("MC1 (H1 perceived reasoning awareness, 1-7) | %s", a),
      n = length(v), mean = mean(v), sd = sd(v), median = median(v),
      q1 = unname(quantile(v, .25)), q3 = unname(quantile(v, .75)))
}
ci <- boot_p(psum[!is.na(MC1)], function(d) d_mean(d, "MC1"))
add("T5 MC", "MC1 | arm difference (reasoning - control)",
    est = d_mean(psum, "MC1"), lo = ci[["lo"]], hi = ci[["hi"]],
    method = "participant bootstrap, stratified by arm")

mc1_str <- qwide[condition == "llm-reasoning"]
mc1_str <- merge(mc1_str, open_p[, .(participant_id, ever_opened)], by = "participant_id", all.x = TRUE)
for (g in c(1L, 0L)) {
  v <- mc1_str[ever_opened == g & !is.na(MC1), MC1]
  add("T5 MC", sprintf("MC1 | reasoning arm, %s", ifelse(g == 1L, "opened >=1 panel", "never opened")),
      n = length(v), mean = mean(v), sd = sd(v), median = median(v))
}
vc <- qwide[condition == "llm-control" & !is.na(MC1), MC1]
add("T5 MC", "MC1 | control arm reference", n = length(vc), mean = mean(vc), sd = sd(vc))

mc2v <- qwide[condition == "llm-reasoning" & !is.na(mc2_correct)]
w <- wilson(sum(mc2v$mc2_correct), nrow(mc2v))
add("T5 MC", "MC2 (H2 named-reason recall) | reasoning arm correct rate",
    n = nrow(mc2v), n_correct = sum(mc2v$mc2_correct),
    est = mean(mc2v$mc2_correct), lo = w[["lo"]], hi = w[["hi"]],
    method = "Wilson 90% interval; chance = 1/4 across the four fixed options (1/3 excluding 'I do not remember')")
mc2_tab <- qwide[condition == "llm-reasoning", .N, by = mc2_raw]
for (i in seq_len(nrow(mc2_tab)))
  add("T5 MC", sprintf("MC2 option chosen: %s", mc2_tab$mc2_raw[i]), n = mc2_tab$N[i])

# =============================================================================
# T6. SINGLE-MEDIATOR DESCRIPTIVE ESTIMATES [Q5 "single-mediator descriptive estimates"]
#   OPERATIONALISATION (ours, stated as such): for each between-arms mechanism
#   scale, the participant-level association between the scale score and the
#   participant's Follow Ratio, pooled across arms, as Pearson r with a
#   participant bootstrap 90% interval. Rationale: the registration commits only
#   to "single-mediator DESCRIPTIVE estimates"; a formal mediation model is
#   explicitly out of scope (main.tex:1890 -- neither manipulation check enters
#   the mediation analysis) and would be unidentified at N=35 with a null-ish
#   a-path. This reports the b-path-shaped association only, and it is a
#   cross-sectional correlation, not a causal path.
# =============================================================================
med_dt <- psum[!is.na(U)]
for (s in c("U", "R", "TR", "SP", "MC1")) {
  d <- med_dt[!is.na(get(s))]
  r_est <- cor(d[[s]], d$follow_ratio)
  set.seed(BOOT_SEED)
  rb <- replicate(BOOT_B, { i <- sample.int(nrow(d), nrow(d), TRUE)
    suppressWarnings(cor(d[[s]][i], d$follow_ratio[i])) })
  rb <- rb[is.finite(rb)]
  add("T6 mediator-shaped descriptives",
      sprintf("%s vs Follow Ratio (participant level, arms pooled)",
              ifelse(s == "MC1", "MC1 (H1)", SCALE_LABEL[[s]])),
      n = nrow(d), est_r = r_est,
      lo = unname(quantile(rb, .05)), hi = unname(quantile(rb, .95)),
      method = "Pearson r, participant bootstrap 90% percentile interval; DESCRIPTIVE ASSOCIATION, not a mediation path")
  # a-path shaped: arm difference on the scale is already in T4; here the c-path
  # analogue within arm is reported so the pair can be read together.
}

# =============================================================================
# T7. ITT vs PER-PROTOCOL SPLIT     [Q8 mandatory]
# =============================================================================
add("T7 ITT/PP", "ITT: reasoning arm (all allocated)",
    n_participants = N_R, n_cells = sum(arm(psum, "llm-reasoning")$n_cells),
    n_follow = sum(arm(psum, "llm-reasoning")$n_follow),
    est = rate(arm(psum, "llm-reasoning"), "n_follow", "n_cells"))
add("T7 ITT/PP", "control arm reference",
    n_participants = N_C, n_cells = sum(arm(psum, "llm-control")$n_cells),
    n_follow = sum(arm(psum, "llm-control")$n_follow),
    est = rate(arm(psum, "llm-control"), "n_follow", "n_cells"))
ci <- boot_p(psum, function(d) d_rate(d, "n_follow", "n_cells"))
add("T7 ITT/PP", "ITT arm difference (reasoning - control)",
    est = d_rate(psum, "n_follow", "n_cells"), lo = ci[["lo"]], hi = ci[["hi"]],
    method = "cluster bootstrap; this is the exploratory descriptive twin of the confirmatory GLMM/Firth estimate, NOT a re-run of it")
for (g in c(1L, 0L)) {
  d <- psum[condition == "llm-reasoning" & ever_opened == g]
  add("T7 ITT/PP", sprintf("PP stratum: reasoning arm, %s",
                           ifelse(g == 1L, "opened >=1 panel", "never opened")),
      n_participants = nrow(d), n_cells = sum(d$n_cells), n_follow = sum(d$n_follow),
      est = sum(d$n_follow) / sum(d$n_cells))
  sub <- psum[(condition == "llm-reasoning" & ever_opened == g) | condition == "llm-control"]
  ci <- boot_p(sub, function(dd) d_rate(dd, "n_follow", "n_cells"))
  add("T7 ITT/PP", sprintf("PP contrast vs control: %s - control",
                           ifelse(g == 1L, "opened", "never-opened")),
      est = d_rate(sub, "n_follow", "n_cells"), lo = ci[["lo"]], hi = ci[["hi"]],
      method = "cluster bootstrap; PP strata are POST-ALLOCATION SELF-SELECTION, randomisation does not hold within them -- not a causal read")
}

# =============================================================================
# T8. 'DO IT'-THEN-DISCONNECT       [Q6]
# =============================================================================
add("T8 disconnect", "Do-it-then-disconnect cells identified",
    n = nrow(disc),
    method = paste("identify_disconnect_cells.R output (disconnect_cells.csv) is empty:",
                   "no cell was accepted-then-abandoned, so the registered 'treat as missing'",
                   "primary coding and the 'recode as 0' sensitivity coding are the same dataset;",
                   "the sensitivity analysis is vacuous, not skipped"))

# =============================================================================
# T10. DATA-QUALITY / FLOW MANIFEST                      [Results 5.1 source table]
#   Not a registered analysis: this is the provenance table behind the flow and
#   quality sentences. Counts only, no contrasts, no intervals.
# =============================================================================
qr_raw  <- fread(file.path(DATA_DIR, "raw", "QuestionnaireResponse.csv"), na.strings = "", encoding = "UTF-8")
part_raw <- fread(file.path(DATA_DIR, "raw", "Participant.csv"), na.strings = "", encoding = "UTF-8")
ua_raw  <- fread(file.path(DATA_DIR, "raw", "UserAction.csv"), na.strings = "", encoding = "UTF-8",
                 colClasses = list(character = c("recommendation", "client_nonce", "portfolio", "action")))
cells_unscreened <- fread(file.path(DATA_DIR, "merged_cells.csv"), na.strings = "", encoding = "UTF-8")

# ---- 10.1 funnel ----
add("T10 manifest", "funnel: QuestionnaireResponse unique users (any submission)",
    n = uniqueN(qr_raw$user))
add("T10 manifest", "funnel: users with a pre-task submission",
    n = uniqueN(quest[phase == "pretask", participant_id]))
add("T10 manifest", "funnel: Participant rows (allocated, condition resolved)",
    n = uniqueN(part_raw$user),
    n_reasoning = part_raw[condition == "llm-reasoning", uniqueN(user)],
    n_control   = part_raw[condition == "llm-control",  uniqueN(user)])
add("T10 manifest", "funnel: cells export before screening",
    n_participants = uniqueN(cells_unscreened$participant_id), n_cells = nrow(cells_unscreened))
add("T10 manifest", "funnel: analysis set after straight-lining screening",
    n_participants = length(PIDS), n_cells = nrow(cells),
    n_reasoning = N_R, n_control = N_C)
add("T10 manifest", "funnel: participants with a post-task submission (analysis set)",
    n = N_Q, of = length(PIDS),
    method = "UCLCEUS0U3T submitted no post-task form; questionnaire analyses run at N=35")

# ---- 10.2 cells per participant ----
npc <- psum$n_cells
add("T10 manifest", "cells per participant (analysis set)",
    n_participants = length(npc), median = median(npc),
    q1 = unname(quantile(npc, .25)), q3 = unname(quantile(npc, .75)),
    min = min(npc), max = max(npc),
    method = "unequal by design: the schedule advances only on decision-bearing rounds, so a participant who exhausted fewer rounds contributes fewer cells")

# ---- 10.3 Q6(b) shortest decision latency ----
#   Decision row = UserAction with amount IS NULL, eligible participant, not
#   ad-hoc paired -- the口径 of identify_disconnect_cells.R lines 119-132
#   (earliest row wins when a cell holds more than one). Latency = that row's
#   created_at minus its Recommendation's created_at. The FK is resolved by the
#   STABLE id in client_nonce ("rec-<id>"), with the (__str__, month) fallback,
#   because the exported FK string collides (identify_disconnect_cells.R:58-96).
#   口径 note: this runs on the UNSCREENED full set (411 cells / 38 participants),
#   matching the registered Q6(b) sub-second exclusion rule, which is applied
#   BEFORE the straight-lining screen.
recs_raw <- fread(file.path(DATA_DIR, "raw", "Recommendation.csv"), na.strings = "", encoding = "UTF-8")
recs_raw[, id := as.integer(id)][, str_repr := paste0(user, ": ", direction, " ", profile_name)]
adhoc_ids <- recs_raw[source == "llm_adhoc", id]
ua_raw[, `:=`(month = as.integer(month), round = as.integer(round), amount = as.numeric(amount), ua_key = .I)]
# suppressWarnings: a free-path write carries a UUID nonce, so the sub() leaves a
# non-numeric string and as.integer warns. The next line sets those to NA anyway,
# and the (__str__, month) fallback below resolves them. Silencing the warning here
# keeps a real coercion problem elsewhere in the script visible.
ua_raw[, rec_id := suppressWarnings(as.integer(sub("^rec-([0-9]+).*$", "\\1", client_nonce)))]
ua_raw[!grepl("^rec-[0-9]+", client_nonce), rec_id := NA_integer_]
need_fb <- ua_raw[!is.na(recommendation) & is.na(rec_id), .(ua_key, month, str_repr = recommendation)]
if (nrow(need_fb) > 0) {
  cand <- merge(need_fb, recs_raw[, .(str_repr, rmonth = as.integer(month), fb_id = id)],
                by = "str_repr", allow.cartesian = TRUE)[month == rmonth]
  hits <- cand[, .(n = .N, fb_id = fb_id[[1]]), by = ua_key]
  if (nrow(hits[n != 1]) > 0 || length(setdiff(need_fb$ua_key, hits$ua_key)) > 0)
    stop("FK fallback ambiguous or unmatched: latency pairing unsafe")
  ua_raw[hits, rec_id := i.fb_id, on = "ua_key"]
}
if (length(setdiff(na.omit(ua_raw$rec_id), recs_raw$id))) stop("client_nonce references unknown Recommendation id")
eligible <- part_raw[!is.na(condition) & condition != "" & !startsWith(user, "TEST_USER"), unique(user)]
dec <- ua_raw[user %in% eligible & is.na(amount) & !(!is.na(rec_id) & rec_id %in% adhoc_ids)]
setorder(dec, user, month, round, created_at)
dec <- dec[, .SD[1], by = .(user, month, round)]
dec <- merge(dec, recs_raw[, .(rec_id = id, rec_created_at = created_at)], by = "rec_id", all.x = TRUE)
dec[, latency_s := as.numeric(difftime(as.POSIXct(created_at, tz = "UTC"),
                                       as.POSIXct(rec_created_at, tz = "UTC"), units = "secs"))]
paired <- dec[!is.na(latency_s)]
add("T10 manifest", "Q6(b) decision-latency pairing completeness (UNSCREENED full set)",
    n_paired = nrow(paired), n_cells = nrow(dec),
    method = "decision row (amount IS NULL, non-adhoc, eligible participant) paired to its Recommendation by stable id; 口径 = unscreened 411-cell set, matching the pre-screen position of the Q6(b) rule")
add("T10 manifest", "Q6(b) shortest decision latency observed",
    min_seconds = min(paired$latency_s), median_seconds = median(paired$latency_s),
    q1 = unname(quantile(paired$latency_s, .25)), q3 = unname(quantile(paired$latency_s, .75)))
add("T10 manifest", "Q6(b) sub-second decisions excluded",
    n = paired[latency_s < 1, .N],
    method = "the registered sub-second settlement-window exclusion found nothing to exclude: the fastest decision is well clear of the threshold, so the rule is vacuous rather than skipped")

# ---- 10.4 post-task duplicate submissions ----
#   ATTRIBUTION: the earliest-submission de-dup rule is NOT a registration clause.
#   It comes from the frozen export/analysis code -- prepare_cells.R applies
#   exactly this rule to the PRE-task instrument ("Pre-registered口径: keep the
#   EARLIEST submission per participant"). T4/T5/T6 extend the same rule to the
#   post-task instrument, which the frozen pipeline never had to handle because
#   it consumes pre-task items only. Cite it as an analysis-stage convention
#   inherited from the frozen code, never as a registered rule.
dup_pids <- quest[phase == "posttask" & participant_id %in% PIDS,
                  .(n_sub = uniqueN(created_at)), by = participant_id][n_sub > 1]
add("T10 manifest", "post-task duplicate submissions detected",
    n_participants = nrow(dup_pids), participants = paste(dup_pids$participant_id, collapse = ","),
    n_rows_dropped = n_pre_dedup - nrow(post),
    method = "earliest submission kept; rule inherited from prepare_cells.R's pre-task de-dup, NOT a registration clause")

# ---- 10.5 technical-fault ledger ----
#   NOT computed to a headline number. The free-text item admits no mechanical
#   fault/no-fault rule (answers range from "no" through "Not really. Sometimes
#   the chatbot wouldn't let me ask about a particular person" to pure UX
#   complaints such as "command palette feels very restricted"). What is
#   mechanically defensible is the three-way tally below. Any headline count of
#   participants "reporting a technical problem" is a HUMAN ADJUDICATION and must
#   carry its own citation to the adjudication record -- it is not derived here.
tech <- post[startsWith(question, "<br>Did you encounter any technical problems")]
tech[, a := tolower(trimws(gsub("\\.$", "", answer)))]
tech[is.na(a), a := ""]
NEGATION <- c("", "no", "n/a", "na", "none", "nope")
add("T10 manifest", "technical-problem item: response tally (analysis set, N=35)",
    n_answered = nrow(tech), n_blank = tech[a == "", .N],
    n_literal_negation = tech[a != "" & a %in% NEGATION, .N],
    n_substantive_text = tech[!(a %in% NEGATION), .N],
    method = paste("literal string rule only (case/whitespace/trailing-period normalised,",
                   "whitelist: blank/no/n\\a/na/none/nope). 'substantive text' is NOT a fault count:",
                   "it mixes genuine faults with UX complaints. A headline fault count is a human",
                   "adjudication and must cite the adjudication record, not this row."))

# =============================================================================
# WRITE
# =============================================================================
lines <- c("# results-exploratory.md -- machine-generated by explore_family.R", "",
           sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
           sprintf("Bootstrap: B = %d, seed = %d, percentile, 90%% intervals.", BOOT_B, BOOT_SEED),
           "", "All results below are EXPLORATORY: no inferential weight, no multiplicity",
           "correction, no significance language. Intervals are 90% estimation intervals.", "")
cur <- ""
for (r in RES) {
  if (r$section != cur) { cur <- r$section; lines <- c(lines, "", sprintf("## %s", cur), "") }
  kv <- r[!(names(r) %in% c("section", "label"))]
  parts <- vapply(names(kv), function(k) {
    v <- kv[[k]]
    sprintf("%s=%s", k, if (is.numeric(v)) fmt(v, if (abs(v) >= 1000) 0 else 3) else as.character(v))
  }, character(1))
  lines <- c(lines, sprintf("- **%s** — %s", r$label, paste(parts, collapse = "; ")))
}
writeLines(lines, file.path(OUT_DIR, "results-exploratory.md"), useBytes = TRUE)
fwrite(psum, file.path(OUT_DIR, "participant_summary.csv"))
fwrite(cells[, .(participant_id, condition, month, round, followed, decision_type, direction,
                 cas, profile_name, complexity_tier, offsetting)],
       file.path(OUT_DIR, "cells_with_complexity.csv"))
fwrite(blocks, file.path(OUT_DIR, "reasoning_blocks.csv"))
fwrite(dec[, .(participant_id = user, month, round, action, rec_id,
               decision_created_at = created_at, rec_created_at, latency_s)],
       file.path(OUT_DIR, "decision_latency.csv"))
cat(sprintf("[done] %d result rows -> %s/results-exploratory.md\n", length(RES), OUT_DIR))
