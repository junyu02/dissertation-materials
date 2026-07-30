# F3 -- follow rate by portfolio-complexity tier.
# Claim: the point-estimate ordering follows the design hinge (offsetting largest),
#        and every arm-difference interval contains 0.
source("/tmp/figures-build/_common.R")

cwc  <- load_cwc(); psum <- load_psum()
stopifnot(nrow(cwc) == 387L)

TIER <- c("single-factor", "two-factor-offsetting", "three-factor-mixed")
TLAB <- c("single-factor" = "Single\nfactor",
          "two-factor-offsetting" = "Two factor,\noffsetting",
          "three-factor-mixed" = "Three factor,\nmixed")

rate <- cwc[, .(n = .N, k = sum(followed), r = sum(followed) / .N),
            by = .(complexity_tier, condition)]
rate[, tier := factor(complexity_tier, levels = TIER, labels = TLAB[TIER])]

chk <- function(t, a) rate[complexity_tier == t & condition == a, r]
stopifnot(near(chk("single-factor",         "llm-reasoning"), 0.826, 1e-3),
          near(chk("single-factor",         "llm-control"),   0.789, 1e-3),
          near(chk("two-factor-offsetting", "llm-reasoning"), 0.724, 1e-3),
          near(chk("two-factor-offsetting", "llm-control"),   0.582, 1e-3),
          near(chk("three-factor-mixed",    "llm-reasoning"), 0.750, 1e-3),
          near(chk("three-factor-mixed",    "llm-control"),   0.679, 1e-3))

# Arm differences: same estimator, same seed/B as explore_family.R T2, so these
# reproduce results-exploratory.md exactly rather than being retyped.
arm  <- function(d, a) d[condition == a]
drate <- function(d, num, den)
  sum(arm(d, "llm-reasoning")[[num]]) / sum(arm(d, "llm-reasoning")[[den]]) -
  sum(arm(d, "llm-control")[[num]])   / sum(arm(d, "llm-control")[[den]])

COLS <- list("single-factor"         = c("n_sf_follow",  "n_sf"),
             "two-factor-offsetting" = c("n_2fo_follow", "n_2fo"),
             "three-factor-mixed"    = c("n_3fm_follow", "n_3fm"))
dif <- rbindlist(lapply(TIER, function(t) {
  cc <- COLS[[t]]
  ci <- boot_p(psum, function(d) drate(d, cc[1], cc[2]))
  data.table(complexity_tier = t, est = drate(psum, cc[1], cc[2]),
             lo = ci[["lo"]], hi = ci[["hi"]])
}))
dif[, tier := factor(complexity_tier, levels = TIER, labels = TLAB[TIER])]

# assertions against exploratory-注册族-20260729/results-exploratory.md, T2
stopifnot(near(dif$est, c(0.037, 0.142, 0.071), 1e-3),
          near(dif$lo,  c(-0.109, -0.035, -0.093), 1e-3),
          near(dif$hi,  c(0.172, 0.331, 0.242), 1e-3))

pL <- ggplot(rate, aes(tier, r, colour = condition, group = condition)) +
  geom_line(linewidth = .5, alpha = .8) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * r),
                vjust = fifelse(condition == "llm-reasoning", -1.2, 2.0)),
            size = 2.5, show.legend = FALSE) +
  scale_discrete_identity(aesthetics = "vjust") +
  scale_colour_manual(values = ARM_COL, labels = ARM_LAB) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.5, 0.90),
                     breaks = seq(0.5, 0.9, .1)) +
  labs(x = NULL, y = "Follow rate (cell level)", tag = "A") +
  theme_fig() + theme(plot.margin = margin(4, 4, 4, 4))

pR <- ggplot(dif, aes(tier, est)) +
  geom_hline(yintercept = 0, linewidth = .4, colour = "grey35") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .10, linewidth = .55, colour = C_REASON) +
  geom_point(size = 2.4, colour = C_REASON) +
  scale_y_continuous(labels = function(x) sprintf("%+.0f", 100 * x),
                     limits = c(-0.15, 0.36), breaks = seq(-0.1, 0.3, .1)) +
  labs(x = NULL, y = "Arm difference, reasoning - control (pp, 90% CI)", tag = "B") +
  theme_fig() +
  theme(axis.title.y = element_text(size = 8), plot.margin = margin(4, 6, 4, 4),
        legend.position = "top", legend.text = element_text(colour = "white"),
        legend.key = element_blank())

# pR carries no legend of its own; the blank strip keeps its panel top aligned
# with pL, which does.
pR <- pR + guides(colour = "none") +
  theme(legend.position = "none", plot.margin = margin(16, 6, 4, 4))

save_fig2(pL, pR, "fig-complexity", w = 5.5, h = 3.0, widths = c(1.05, 1))
cat("F3 ok | tier rates and arm-difference intervals match results-exploratory.md T2\n")
print(dif[, .(complexity_tier, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3))])
