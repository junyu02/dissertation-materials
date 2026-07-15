#!/usr/bin/env Rscript
# Follow-Ratio confirmatory analysis (B2).
#
# Primary (confirmatory): binomial GLMM
#   followed ~ condition + literacy + is_HCI + factor(month) + factor(round) + (1 | participant_id)
#   logit link, Laplace (nAGQ = 1);
#   report OR + 90% CI (Wald primary, profile-likelihood as sensitivity).
#
# `followed` is the direction-aware DV (Plan X): a cell is followed==1 iff the
# participant took the direction-compliant adoption, i.e. decision_type is
# 'follow' (buy/hold advice adopted) OR 'unfollow' (divest advice adopted).
# `direction` labels the advice, not the choice, so a declined unfollow keeps
# direction=='unfollow' with decision_type=='decline'. cas is the
# conditional-on-follow share and is defined ONLY when decision_type=='follow'.
#
# Pre-registered covariates (outline B7): literacy (Lusardi Big-Three, 0-3) and
# is_HCI (0/1) are participant-level; month enters as FIXED dummies via
# factor(month) rather than a second random effect because it has < 8 levels
# (5 here) -- the pre-registered fallback rule for low-cardinality grouping
# factors (outline SS 3.9). On the real path prepare_cells.R derives literacy /
# is_HCI from the pre-task questionnaire before this script runs.
#
# Pre-registered fallback: if the GLMM fails to converge under CODED criteria
#   (lme4 convergence message, max scaled gradient > tol, or non-zero optimizer
#   return code) OR errors outright, collapse to participant level (per-participant
#   follow counts) and fit a Firth-penalised logistic regression (logistf). On
#   the fallback path the primary CI is logistf's penalised-profile CI (the
#   recommended interval under Firth), with Wald kept as a parity macro. The
#   chosen path is recorded in \FRfitPath so downstream text never has to guess.
#
# CAS (secondary, descriptive only): conditional-on-follow share, per-arm
#   mean / median / IQR. No model is fit to CAS.
#
# Set FORCE_FALLBACK=1 in the environment to bypass the GLMM and exercise the
# Firth branch (used to prove both paths emit an identical macro set).
#
# Usage: Rscript fit_glmm.R [input_csv] [output_tex]

suppressMessages({
  library(data.table)
  library(lme4)
  library(logistf)
})

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(here) == 0 || here == "") here <- "."
in_csv  <- if (length(args) >= 1) args[[1]] else file.path(here, "data", "synth_cells.csv")
out_tex <- if (length(args) >= 2) args[[2]] else file.path(here, "results-values.tex")

GRAD_TOL <- 0.002   # threshold borrowed from lme4 checkConv's default tol for the relative
                    # gradient (solve(H, g)); see README for the joint convergence criterion
Z90 <- qnorm(0.95)  # 1.6449: 90% two-sided Wald multiplier

## ---- 1. load ---------------------------------------------------------------
dt <- fread(in_csv, na.strings = "")
stopifnot(nrow(dt) > 0)
req <- c("participant_id", "condition", "month", "round", "followed", "decision_type",
         "direction", "available_snapshot", "invested", "cas", "literacy", "is_HCI")
missing_cols <- setdiff(req, names(dt))
if (length(missing_cols)) stop(sprintf("input missing columns: %s", paste(missing_cols, collapse = ", ")))

# Domain checks BEFORE integer coercion: as.integer(1.7) would silently truncate
# to 1 and then pass a post-coercion %in% c(0,1) check (R2 cross-check finding).
stopifnot(
  all(as.numeric(dt$followed) %in% c(0, 1)),
  all(as.numeric(dt$literacy) %in% 0:3),
  all(as.numeric(dt$is_HCI) %in% c(0, 1))
)
dt[, followed := as.integer(followed)]
dt[, cas := as.numeric(cas)]
dt[, invested := as.numeric(invested)]
dt[, available_snapshot := as.numeric(available_snapshot)]
dt[, month := as.integer(month)]
dt[, round := as.integer(round)]   # enters the formula as factor(round): NA would be na.omit'ed silently
dt[, literacy := as.integer(literacy)]
dt[, is_HCI := as.integer(is_HCI)]

## ---- 1b. hard input guards ------------------------------------------------
# data.table `dt[cond]` treats NA in `cond` as FALSE, so a NA in ANY audited
# column would slip past the structural audit below (an NA decision_type row
# appears in NO per-type count; an NA direction skips the consistency checks)
# and then be silently dropped by glmer's na.omit. Fail loud here, before any
# i-expression, and whitelist every enum so an unknown label cannot hide either.
DECISION_TYPES <- c("follow", "unfollow", "decline", "zero_optout")
DIRECTIONS     <- c("follow", "unfollow")
stopifnot(
  !anyNA(dt$followed),
  all(dt$followed %in% c(0L, 1L)),
  !anyNA(dt$condition),
  setequal(unique(dt$condition), c("llm-reasoning", "llm-control")),
  !anyNA(dt$participant_id),
  !anyNA(dt$month),
  !anyNA(dt$round),
  !anyNA(dt$decision_type), all(dt$decision_type %in% DECISION_TYPES),
  !anyNA(dt$direction),     all(dt$direction %in% DIRECTIONS),
  # available == 0 is a LEGAL state (a fully-invested participant declining a
  # recommendation — real data has such rows); only negative/NA/Inf are corrupt.
  # A pathological follow at available==0 is caught downstream: cas = invested/0
  # = Inf trips the (0,1] range alarm.
  !anyNA(dt$available_snapshot), all(is.finite(dt$available_snapshot)),
  all(dt$available_snapshot >= 0),
  !anyNA(dt$literacy), all(dt$literacy %in% 0:3),
  !anyNA(dt$is_HCI),   all(dt$is_HCI %in% c(0L, 1L))
)
# One row per (participant, month, round): a duplicated cell (double export,
# upstream double-write) would silently double-count inside the GLMM.
if (anyDuplicated(dt, by = c("participant_id", "month", "round")) > 0)
  stop("duplicate (participant_id, month, round) cells in input")
if (dt[, uniqueN(condition), by = participant_id][V1 > 1, .N] > 0)
  stop("participant(s) span more than one condition")
# literacy / is_HCI are participant-level covariates: exactly one value per participant.
if (dt[, uniqueN(literacy), by = participant_id][V1 > 1, .N] > 0)
  stop("participant(s) span more than one literacy value")
if (dt[, uniqueN(is_HCI), by = participant_id][V1 > 1, .N] > 0)
  stop("participant(s) span more than one is_HCI value")

## ---- 2. structural-NA vs true-missing audit -------------------------------
# Bidirectional contract (Plan X):
#   * followed == 1  iff  decision_type in {follow, unfollow}   (both are adoptions)
#   * cas is defined iff  decision_type == 'follow'             (a conditional-on-follow share)
# A NULL cas on any non-follow row (decline / zero_optout / unfollow) is a
# STRUCTURAL NA -- an active choice, never missing data. `unfollow` records a
# POSITIVE divestment `invested` but no cas (a share of available is undefined
# for a withdrawal). A NULL cas on a follow row is TRUE MISSING and must not
# occur in behavioural data.
n_follow   <- dt[decision_type == "follow", .N]
n_unfollow <- dt[decision_type == "unfollow", .N]
n_decline  <- dt[decision_type == "decline", .N]
n_optout   <- dt[decision_type == "zero_optout", .N]
# structural NA = any non-follow (decline+zero_optout+unfollow) with a NULL cas
n_structural_na <- dt[decision_type != "follow" & is.na(cas), .N]
n_true_missing  <- dt[decision_type == "follow" & is.na(cas), .N]
n_contradiction <- dt[decision_type != "follow" & !is.na(cas), .N]   # share recorded off a follow

# decision_type must agree with followed, and each type must carry the right
# `invested` encoding: decline -> NULL, zero_optout -> 0, unfollow -> positive.
n_type_mismatch <- dt[(followed == 1) != (decision_type %in% c("follow", "unfollow")), .N]
n_decline_bad   <- dt[decision_type == "decline" & !is.na(invested), .N]
n_optout_bad    <- dt[decision_type == "zero_optout" & (is.na(invested) | invested != 0), .N]
n_unfollow_bad  <- dt[decision_type == "unfollow" & (is.na(invested) | invested <= 0), .N]
# a follow adoption must itself carry a positive invested (checked directly, not
# only via the cas identity, which an NA invested would NA-skip)
n_follow_bad    <- dt[decision_type == "follow" & (is.na(invested) | invested <= 0), .N]
# contract-undefined combination: the deployed platform has never produced a
# zero_optout on an unfollow-direction cell; if real data ever does, stop and
# adjudicate rather than silently coding it (R2 cross-check finding).
n_undef_combo   <- dt[decision_type == "zero_optout" & direction == "unfollow", .N]

# direction must agree with decision_type on the two adoption types (an adopted
# recommendation cannot point the opposite way). decline / zero_optout can carry
# either direction (they decline whatever was advised), so they are not checked.
n_dir_unfollow_bad <- dt[decision_type == "unfollow" & direction != "unfollow", .N]
n_dir_follow_bad   <- dt[decision_type == "follow"   & direction != "follow",   .N]

# a present cas (follow rows only) must be a valid share and equal invested/available
n_cas_range_bad <- dt[decision_type == "follow" & !is.na(cas) & (cas <= 0 | cas > 1), .N]
n_cas_calc_bad  <- dt[decision_type == "follow" & !is.na(cas) &
                        abs(cas - invested / available_snapshot) >= 1e-9, .N]

alarms <- c()
if (n_true_missing     > 0) alarms <- c(alarms, sprintf("%d follow rows with NULL cas (true missing)", n_true_missing))
if (n_contradiction    > 0) alarms <- c(alarms, sprintf("%d non-follow rows with non-NULL cas", n_contradiction))
if (n_type_mismatch    > 0) alarms <- c(alarms, sprintf("%d rows where followed disagrees with decision_type", n_type_mismatch))
if (n_decline_bad      > 0) alarms <- c(alarms, sprintf("%d decline rows with non-NULL invested", n_decline_bad))
if (n_optout_bad       > 0) alarms <- c(alarms, sprintf("%d zero_optout rows without invested==0", n_optout_bad))
if (n_unfollow_bad     > 0) alarms <- c(alarms, sprintf("%d unfollow rows without a positive invested", n_unfollow_bad))
if (n_follow_bad       > 0) alarms <- c(alarms, sprintf("%d follow rows without a positive invested", n_follow_bad))
if (n_undef_combo      > 0) alarms <- c(alarms, sprintf("%d zero_optout rows on an unfollow-direction cell (contract-undefined combination)", n_undef_combo))
if (n_dir_unfollow_bad > 0) alarms <- c(alarms, sprintf("%d unfollow rows with direction != unfollow", n_dir_unfollow_bad))
if (n_dir_follow_bad   > 0) alarms <- c(alarms, sprintf("%d follow rows with direction != follow", n_dir_follow_bad))
if (n_cas_range_bad    > 0) alarms <- c(alarms, sprintf("%d rows with cas outside (0,1]", n_cas_range_bad))
if (n_cas_calc_bad     > 0) alarms <- c(alarms, sprintf("%d rows where cas != invested/available_snapshot", n_cas_calc_bad))
if (length(alarms)) stop("DATA INTEGRITY ALARM:\n  - ", paste(alarms, collapse = "\n  - "))
cat(sprintf("[audit] follow=%d  unfollow=%d  decline=%d  zero_optout=%d  structural-NA cas=%d  true-missing cas=%d  (clean)\n",
            n_follow, n_unfollow, n_decline, n_optout, n_structural_na, n_true_missing))

## ---- 3. model matrix -------------------------------------------------------
dt[, condition := relevel(factor(condition), ref = "llm-control")]
dt[, participant_id := factor(participant_id)]
TERM <- "conditionllm-reasoning"   # fixed-effect term: reasoning vs control
n_total <- nlevels(dt$participant_id)
n_cells <- nrow(dt)

## ---- 4. primary GLMM with coded convergence gate --------------------------
force_fallback <- identical(Sys.getenv("FORCE_FALLBACK"), "1")

fit_glmm <- function() {
  # month enters as fixed dummies via factor(month): < 8 levels, so it is a
  # fixed low-cardinality control, not a second random effect (pre-registered,
  # outline SS 3.9). round enters as fixed dummies too — SS 3.9's non-stationarity
  # /learning clause ("round modelled as a factor"), adopted 2026-07-15 after the
  # R2 cross-check flagged the outline-internal tension. literacy / is_HCI are
  # participant-level fixed covariates.
  m <- glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) + (1 | participant_id),
             data = dt, family = binomial, nAGQ = 1)
  conv_msgs <- m@optinfo$conv$lme4$messages
  opt_code  <- m@optinfo$conv$opt   # optimizer return code (0 = success); not carried in $lme4$messages
  dd <- m@optinfo$derivs
  scaled_grad <- tryCatch(max(abs(solve(dd$Hessian, dd$gradient))),
                          error = function(e) NA_real_)
  # NOTE (pre-registered, coded criterion): a singular/boundary fit (RE variance
  # = 0) does NOT trigger the fallback — it degenerates to a pooled logistic,
  # which is statistically defensible; singularity is computed and printed for
  # the record but is deliberately not part of the convergence gate.
  converged <- length(conv_msgs) == 0 &&
               is.finite(scaled_grad) && scaled_grad <= GRAD_TOL &&
               isTRUE(opt_code == 0)
  list(model = m, converged = converged, scaled_grad = scaled_grad,
       conv_msgs = conv_msgs, opt_code = opt_code, singular = isSingular(m))
}

glmm <- NULL
if (!force_fallback) {
  glmm <- tryCatch(fit_glmm(), error = function(e) {
    cat(sprintf("[glmm] hard error: %s\n", conditionMessage(e)))
    NULL
  })
}

use_glmm <- !force_fallback && !is.null(glmm) && isTRUE(glmm$converged)

results <- list()
if (use_glmm) {
  m <- glmm$model
  beta <- fixef(m)[[TERM]]
  se   <- sqrt(diag(vcov(m)))[[TERM]]
  prof <- tryCatch(
    withCallingHandlers(
      confint(m, parm = TERM, method = "profile", level = 0.90),
      warning = function(w) {   # keep the estimate, but LOG the warning (a non-monotone
        cat(sprintf("[profile-CI warning] %s\n", conditionMessage(w)))   # profile is worth seeing) rather than silently muffling it
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) c(NA_real_, NA_real_)
  )
  results$fit_path     <- "GLMM (Laplace)"
  results$or_point     <- exp(beta)
  results$or_low       <- exp(beta - Z90 * se)   # primary CI = Wald 90%
  results$or_high      <- exp(beta + Z90 * se)
  results$or_low_wald  <- results$or_low          # explicit Wald pair (== primary on this path)
  results$or_high_wald <- results$or_high
  results$or_low_pf    <- exp(prof[[1]])          # likelihood-profile 90% (sensitivity)
  results$or_high_pf   <- exp(prof[[2]])
  results$grad_max     <- glmm$scaled_grad
  cat(sprintf("[glmm] converged  scaled_grad=%.2e  singular=%s\n",
              glmm$scaled_grad, glmm$singular))
} else {
  # Firth fallback: collapse to per-participant follow counts, then Firth logistic.
  # condition is between-subjects (constant within participant), so the aggregated
  # binomial recovers the arm contrast; weights carry the per-participant trial counts.
  reason <- if (force_fallback) "FORCE_FALLBACK=1"
            else if (is.null(glmm)) "GLMM errored"
            else sprintf("non-convergence (scaled_grad=%.2e, msgs=%d, opt=%s)",
                         glmm$scaled_grad, length(glmm$conv_msgs), format(glmm$opt_code))
  cat(sprintf("[fallback] Firth path taken: %s\n", reason))

  # condition / literacy / is_HCI are all constant within participant, so grouping
  # by them keeps one row per participant. month and round dummies CANNOT survive the
  # collapse (they vary within participant) -- their loss is a known degradation
  # of this degraded path; the fallback trades month adjustment for guaranteed
  # convergence.
  agg <- dt[, .(n_follow = sum(followed), n_total = .N),
            by = .(participant_id, condition, literacy, is_HCI)]
  firth_long <- rbindlist(lapply(seq_len(nrow(agg)), function(i) {
    data.table(y = c(1L, 0L),
               w = c(agg$n_follow[i], agg$n_total[i] - agg$n_follow[i]),
               condition = agg$condition[i],
               literacy = agg$literacy[i],
               is_HCI = agg$is_HCI[i])
  }))
  firth_long <- firth_long[w > 0]
  firth_long[, condition := relevel(factor(condition), ref = "llm-control")]
  fit <- logistf(y ~ condition + literacy + is_HCI, data = firth_long, weights = w, alpha = 0.10)
  idx  <- which(names(coef(fit)) == TERM)
  beta <- coef(fit)[[idx]]
  se   <- sqrt(diag(vcov(fit)))[[idx]]
  results$fit_path     <- "Firth (participant-collapsed)"
  results$or_point     <- exp(beta)
  results$or_low       <- exp(fit$ci.lower[[idx]])     # primary CI = logistf penalised-profile 90%
  results$or_high      <- exp(fit$ci.upper[[idx]])
  results$or_low_pf    <- results$or_low               # profile pair (== primary on this path)
  results$or_high_pf   <- results$or_high
  results$or_low_wald  <- exp(beta - Z90 * se)         # Wald 90% (parity / sensitivity)
  results$or_high_wald <- exp(beta + Z90 * se)
  results$grad_max     <- NA_real_

  # pre-registered clustering sensitivity (see README analysis plan): compare
  # per-participant follow proportions with a Wilcoxon rank-sum. Reported only
  # when the fallback fires; the pooled Firth CI above is not clustering-adjusted.
  agg[, prop := n_follow / n_total]
  wcx <- tryCatch(wilcox.test(prop ~ condition, data = agg, exact = FALSE),
                  error = function(e) NULL)
  if (!is.null(wcx))
    cat(sprintf("[fallback] participant-level Wilcoxon (clustering sensitivity): W=%.1f p=%.4f\n",
                wcx$statistic, wcx$p.value))
}

## ---- 5. descriptive: follow rates + CAS ------------------------------------
rate <- function(arm) dt[condition == arm, mean(followed)]
rate_reasoning <- rate("llm-reasoning")
rate_control   <- rate("llm-control")

# cas is conditional-on-follow only: followed==1 now also covers unfollow
# adoptions (NULL cas), so filter on decision_type, not followed.
cas_dt <- dt[decision_type == "follow" & !is.na(cas)]
cas_summary <- function(arm) {
  v <- cas_dt[condition == arm, cas]
  if (length(v) == 0)   # an arm with zero follow adoptions: emit NA macros, don't NaN/error
    return(list(mean = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_, n = 0L))
  list(mean = mean(v), median = median(v),
       q1 = as.numeric(quantile(v, 0.25)), q3 = as.numeric(quantile(v, 0.75)), n = length(v))
}
cas_r <- cas_summary("llm-reasoning")
cas_c <- cas_summary("llm-control")

## ---- 6. emit LaTeX macros --------------------------------------------------
fnum <- function(x, d = 2) if (is.na(x)) "NA" else formatC(x, format = "f", digits = d)
macro <- function(name, value) sprintf("\\newcommand{\\%s}{%s}", name, value)

lines <- c(
  "% results-values.tex -- auto-generated by fit_glmm.R. Do not edit by hand.",
  sprintf("%% source: %s", basename(in_csv)),
  sprintf("%% generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  macro("FRfitPath",  results$fit_path),
  macro("FRnTotal",   n_total),
  macro("FRnCells",   n_cells),
  "",
  macro("FRORpoint",  fnum(results$or_point, 2)),
  macro("FRORlow",    fnum(results$or_low, 2)),
  macro("FRORhigh",   fnum(results$or_high, 2)),
  macro("FRORlowProfile",  fnum(results$or_low_pf, 2)),
  macro("FRORhighProfile", fnum(results$or_high_pf, 2)),
  macro("FRORlowWald",     fnum(results$or_low_wald, 2)),
  macro("FRORhighWald",    fnum(results$or_high_wald, 2)),
  macro("FRgradMax",  if (is.na(results$grad_max)) "NA" else formatC(results$grad_max, format = "e", digits = 2)),
  "",
  macro("FRrateReasoning", fnum(rate_reasoning, 3)),
  macro("FRrateControl",   fnum(rate_control, 3)),
  "",
  macro("CASmeanReasoning",   fnum(cas_r$mean, 3)),
  macro("CASmeanControl",     fnum(cas_c$mean, 3)),
  macro("CASmedianReasoning", fnum(cas_r$median, 3)),
  macro("CASmedianControl",   fnum(cas_c$median, 3)),
  macro("CASqOneReasoning",   fnum(cas_r$q1, 3)),
  macro("CASqThreeReasoning", fnum(cas_r$q3, 3)),
  macro("CASqOneControl",     fnum(cas_c$q1, 3)),
  macro("CASqThreeControl",   fnum(cas_c$q3, 3)),
  "",
  macro("FRnUnfollow",     n_unfollow),
  macro("FRnDecline",      n_decline),
  macro("FRnOptout",       n_optout),
  macro("FRnStructuralNA", n_structural_na),
  macro("FRnTrueMissing",  n_true_missing)
)
writeLines(lines, out_tex)

cat("\n=== results ===\n")
cat(sprintf("fit path : %s\n", results$fit_path))
cat(sprintf("OR       : %s  (90%% primary CI %s, %s)  [Wald on GLMM, penalised-profile on Firth]\n",
            fnum(results$or_point), fnum(results$or_low), fnum(results$or_high)))
cat(sprintf("OR 90%% profile CI: %s, %s   Wald CI: %s, %s\n",
            fnum(results$or_low_pf), fnum(results$or_high_pf),
            fnum(results$or_low_wald), fnum(results$or_high_wald)))
cat(sprintf("follow rate : reasoning=%s control=%s\n",
            fnum(rate_reasoning, 3), fnum(rate_control, 3)))
cat(sprintf("CAS mean    : reasoning=%s control=%s\n",
            fnum(cas_r$mean, 3), fnum(cas_c$mean, 3)))
cat(sprintf("wrote macros -> %s\n", out_tex))
