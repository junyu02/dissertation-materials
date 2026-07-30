#!/usr/bin/env Rscript
# cluster_bootstrap_canonical.R -- participant-level cluster bootstrap for the
# Firth fallback estimate (POST-HOC sensitivity; NOT part of the frozen
# registered pipeline commit 352792c).
#
# PROVENANCE
#   Algorithm identical to the two independent implementations that produced
#   the reported intervals (executor run 2026-07-27; claude-worker independent
#   rewrite, archived here as boot_seed.R, which reproduced them bit-for-bit
#   under the same seed). This canonical version differs ONLY in that input/
#   output paths are arguments instead of hard-coded /tmp paths (Codex round-2
#   M3 requirement: archived code must carry seed, B, and no /tmp hard-coding).
#
# REPORTED NUMBERS THIS SCRIPT REPRODUCES (seed 20260727, B = 2000):
#   N=38  merged_cells.csv           -> 90% percentile CI [0.859, 2.999]
#   N=36  merged_cells_screened.csv  -> 90% percentile CI [0.821, 2.967]
#   (Monte-Carlo jitter at B=2000 is roughly +/-0.02-0.04: do not quote more
#    than two decimals in the manuscript.)
#
# Usage:
#   Rscript cluster_bootstrap_canonical.R <merged_cells csv> [B] [seed]
#   e.g. Rscript cluster_bootstrap_canonical.R \
#          ../data/real-20260727/merged_cells_screened.csv 2000 20260727

suppressMessages({library(data.table); library(logistf)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: cluster_bootstrap_canonical.R <merged_cells.csv> [B] [seed]")
infile <- args[[1]]
B      <- if (length(args) >= 2) as.integer(args[[2]]) else 2000L
seed   <- if (length(args) >= 3) as.integer(args[[3]]) else 20260727L

TERM <- "conditionllm-reasoning"

d <- fread(infile, na.strings = "")
d[, followed := as.integer(followed)]
d[, condition := relevel(factor(condition), ref = "llm-control")]
d[, participant_id := factor(participant_id)]

# Firth fit on participant-collapsed counts (mirrors the registered fallback:
# collapse to per-participant follow counts, keep participant-level covariates,
# expand to weighted success/failure rows; month/round cannot survive collapse).
firth_of <- function(dd) {
  agg <- dd[, .(n_follow = sum(followed), n_total = .N),
            by = .(participant_id, condition, literacy, is_HCI)]
  long <- rbindlist(lapply(seq_len(nrow(agg)), function(i)
    data.table(y = c(1L, 0L),
               w = c(agg$n_follow[i], agg$n_total[i] - agg$n_follow[i]),
               condition = agg$condition[i],
               literacy = agg$literacy[i], is_HCI = agg$is_HCI[i])))
  long <- long[w > 0]
  long[, condition := relevel(factor(condition), ref = "llm-control")]
  f <- logistf(y ~ condition + literacy + is_HCI, data = long,
               weights = w, alpha = 0.10)
  coef(f)[[which(names(coef(f)) == TERM)]]
}

set.seed(seed)
ids <- unique(d$participant_id)
out <- rep(NA_real_, B)
for (b in seq_len(B)) {
  take <- sample(ids, length(ids), replace = TRUE)
  # re-label resampled clusters so a participant drawn twice contributes twice
  boot <- rbindlist(lapply(seq_along(take), function(i) {
    x <- d[participant_id == take[i]]
    x[, participant_id := paste0(take[i], "_", i)]
    x
  }))
  out[b] <- tryCatch(firth_of(boot), error = function(e) NA_real_)
}

ok <- !is.na(out)
qs <- quantile(exp(out[ok]), c(.05, .50, .95))
cat(sprintf("input=%s  B=%d  seed=%d  effective=%d\n", infile, B, seed, sum(ok)))
cat(sprintf("OR percentile 05/50/95 = %.6f / %.6f / %.6f\n", qs[1], qs[2], qs[3]))
