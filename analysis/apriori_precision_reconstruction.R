#!/usr/bin/env Rscript
# ============================================================================
# A-priori precision calculation -- REPRODUCIBLE RECONSTRUCTION
#
# WHAT THIS IS. The a-priori precision calculation for the primary Follow-Ratio
# contrast was performed and recorded BEFORE data collection: the expected 90%
# CI half-width on the log-odds scale ("~ +/- 0.6") is written into the frozen
# study outline v3.1, dated 2026-06-27, sections 3.5 and 3.9. Collection began
# 2026-07-14/15. The calculation itself was done outside version control and no
# script was archived at the time. This file re-derives that number from the
# same literature-prior inputs so the reported figure is reproducible.
#
# WHAT THIS IS NOT. This script is NOT part of the registered, frozen analysis
# code (repo commit 352792c, AsPredicted filing 2026-07-15). It fits no model,
# reads no file, and touches nothing the confirmatory pipeline produces. Every
# input below is a published prior or a fixed design constant; NO collected data
# of any kind enters this file. It is documentation of a design-stage
# calculation, added after the fact for transparency, and it changes no
# pre-specified analytic degree of freedom.
#
# Usage: Rscript apriori_precision_reconstruction.R
# Requires: base R only (no packages).
# ============================================================================

## ---- inputs: literature priors --------------------------------------------
# Milana (2023), same-platform Follow-Ratio base rates, verbatim from the PDF:
P_BASE  <- 0.68   # mean Follow Ratio (proportion of recommendations adopted)
SD_OBS  <- 0.22   # SD of the per-participant Follow Ratio (observed, between-participant)

## ---- inputs: design constants (this study) --------------------------------
N_ARM   <- c(15, 20)  # target N per arm, between-subjects 1x2 (outline v3.1 SS 3.9)
# K_CELLS is the SCHEDULED round count, used here as an upper bound on the
# decision cells a participant can contribute. The deployed Follow Ratio counts
# only decision-bearing rounds, so a participant's real denominator is <= 20.
# This affects the Monte Carlo cross-check below, which is therefore an
# idealised sanity check and NOT isomorphic to the deployed outcome; the
# analytical delta-method half-widths reported in the manuscript do not depend
# on K_CELLS. (Caveat added 2026-07-26 after external review.)
K_CELLS <- 20         # scheduled rounds per participant = 5 months x 4 rounds
                      # (as-deployed; cf. analysis pipeline README, "40 participants
                      #  x 20 cells = 800 rows")
CONF    <- 0.90       # 90% CI, matching the confirmatory report (fit_glmm.R Z90)
Z       <- qnorm(1 - (1 - CONF) / 2)   # 1.6449

## ---- estimand --------------------------------------------------------------
# The confirmatory contrast is the reasoning-vs-control difference in Follow
# Ratio on the LOGIT scale (fit_glmm.R reports exp(beta) with a 90% Wald CI, so
# beta itself is the log-odds contrast and the half-width below is exactly
# Z * SE(beta)).
#
# The sampling unit is the PARTICIPANT, not the cell. SD_OBS is the spread of
# participant-level Follow Ratios, so SD_OBS/sqrt(n) is already the
# cluster-honest SE of an arm mean: it absorbs both between-participant
# heterogeneity and the within-participant binomial noise across the K_CELLS
# rounds. This is the same stance the frozen pipeline codes as
# (1 | participant_id) and states as "the condition effect is estimated on ~N
# clusters, not 600-800 rows".

## ---- method A: delta method (analytic) -------------------------------------
# d logit(p) / dp = 1 / (p(1-p))  =>  SE(logit p_bar) ~ SE(p_bar) / (p(1-p)).
# Two independent arms, equal n and (under the planning assumption of a small
# effect near the base rate) equal SD  =>  SE(contrast) = sqrt(2) * SE(logit p_bar).
delta_hw <- function(n, p = P_BASE, sd = SD_OBS) {
  se_prop  <- sd / sqrt(n)                 # SE of one arm's mean Follow Ratio
  se_logit <- se_prop / (p * (1 - p))      # delta-method transfer to log-odds
  Z * sqrt(2) * se_logit                   # 90% half-width of the two-arm contrast
}

# Same quantity left on the raw proportion scale, purely to identify which scale
# a quoted number lives on (a risk-difference half-width is ~4.6x smaller here).
delta_hw_prop <- function(n, sd = SD_OBS) Z * sqrt(2) * sd / sqrt(n)

## ---- method B: Monte Carlo (same estimator, no normal approximation) --------
# Beta-binomial participants: p_i ~ Beta(mean = P_BASE), y_i ~ Binom(K_CELLS, p_i).
# The latent Beta SD is set so the SIMULATED observed participant-ratio SD equals
# SD_OBS -- otherwise the binomial layer would stack extra noise on top of a
# figure that already includes it:
#   var_obs = var_latent + E[p(1-p)]/K = var_latent*(1 - 1/K) + P_BASE(1-P_BASE)/K
var_obs    <- SD_OBS^2
var_latent <- (var_obs - P_BASE * (1 - P_BASE) / K_CELLS) / (1 - 1 / K_CELLS)
stopifnot(var_latent > 0, var_latent < P_BASE * (1 - P_BASE))
conc  <- P_BASE * (1 - P_BASE) / var_latent - 1   # Beta concentration a+b
BETA_A <- P_BASE * conc
BETA_B <- (1 - P_BASE) * conc

sim_hw <- function(n, reps = 20000L, seed = 20260627L) {
  set.seed(seed)   # seeded with the outline-v3.1 date; fixed for reproducibility
  draw_arm <- function() {
    p <- matrix(rbeta(n * reps, BETA_A, BETA_B), nrow = reps)
    matrix(rbinom(n * reps, K_CELLS, p), nrow = reps) / K_CELLS   # per-participant Follow Ratio
  }
  a <- draw_arm(); b <- draw_arm()
  ma <- rowMeans(a); mb <- rowMeans(b)
  # per-arm SE of the mean ratio, then delta-transform each arm at its own mean
  se <- function(x) apply(x, 1, sd) / sqrt(n)
  se_logit_a <- se(a) / (ma * (1 - ma))
  se_logit_b <- se(b) / (mb * (1 - mb))
  hw <- Z * sqrt(se_logit_a^2 + se_logit_b^2)
  hw[is.finite(hw)]
}

## ---- report ----------------------------------------------------------------
cat(sprintf("a-priori precision reconstruction  (%.0f%% CI, log-odds scale)\n", 100 * CONF))
cat(sprintf("inputs: Milana-2023 base rate p=%.2f, participant-level SD=%.2f; %d cells/participant; z=%.4f\n",
            P_BASE, SD_OBS, K_CELLS, Z))
cat(sprintf("simulation: Beta(%.4f, %.4f) x Binom(%d), latent SD=%.4f -> observed ratio SD=%.2f\n\n",
            BETA_A, BETA_B, K_CELLS, sqrt(var_latent), SD_OBS))

cat(sprintf("%-8s %-14s %-24s %-16s %s\n",
            "N/arm", "delta (logit)", "simulated median (90% band)", "as OR factor", "proportion scale"))
for (n in N_ARM) {
  d <- delta_hw(n)
  s <- sim_hw(n)
  cat(sprintf("%-8d %-14.3f %-24s %-16s %.3f\n",
              n, d,
              sprintf("%.3f (%.3f-%.3f)", median(s), quantile(s, 0.05), quantile(s, 0.95)),
              sprintf("x/ %.2f", exp(d)),
              delta_hw_prop(n)))
}

## ---- self-check ------------------------------------------------------------
# The two methods estimate the same quantity by different routes; a >15%
# disagreement means one of them is wrong, not that the design changed.
for (n in N_ARM) {
  d <- delta_hw(n); s <- median(sim_hw(n))
  stopifnot(abs(s - d) / d < 0.15)
}
# Guard the reported claim itself: the reconstruction must land near the
# +/- 0.6 recorded in outline v3.1 SS 3.5 at the lower end of the N range.
stopifnot(abs(delta_hw(15) - 0.6) < 0.2)
cat("\n[self-check] delta vs simulation agree within 15%; N=15 half-width within 0.2 of the outline's +/-0.6\n")
