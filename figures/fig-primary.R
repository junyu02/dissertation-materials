# F1 -- participant-level follow ratio, both arms.
# Claim: the raw arm gap (+8.1pp) is small next to between-participant spread;
#        the two arms overlap heavily.
source("/tmp/figures-build/_common.R")

cells <- load_cells(); psum <- load_psum()

cell_rate <- cells[, .(r = sum(followed) / .N), by = condition]
stopifnot(near(cell_rate[condition == "llm-reasoning", r], 0.7740),
          near(cell_rate[condition == "llm-control",   r], 0.6927))

arm <- psum[, .(m = mean(follow_ratio), n = .N), by = condition]
stopifnot(near(arm[condition == "llm-reasoning", m], 0.7417),
          near(arm[condition == "llm-control",   m], 0.6648),
          arm[, all(n == 18L)])

ci <- rbindlist(lapply(c("llm-reasoning", "llm-control"), function(a) {
  b <- boot_g(psum[condition == a], function(d) mean(d$follow_ratio))
  data.table(condition = a, lo = b[["lo"]], hi = b[["hi"]])
}))
arm <- merge(arm, ci, by = "condition")

psum[, xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]
arm[,  xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]

set.seed(BOOT_SEED)
p <- ggplot(psum, aes(xarm, follow_ratio, colour = condition)) +
  geom_point(aes(size = n_cells), alpha = .45, stroke = 0,
             position = position_jitter(width = .11, height = 0, seed = BOOT_SEED)) +
  geom_errorbar(data = arm, aes(y = m, ymin = lo, ymax = hi), width = .07, linewidth = .6,
                position = position_nudge(x = .28)) +
  geom_point(data = arm, aes(y = m), shape = 18, size = 3.4,
             position = position_nudge(x = .28)) +
  scale_colour_manual(values = ARM_COL, labels = ARM_LAB, guide = "none") +
  scale_size_continuous(range = c(1.1, 3.4), breaks = c(6, 10, 14),
                        name = "cells per participant") +
  scale_x_discrete(labels = ARM_LAB) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.02),
                     breaks = seq(0, 1, .25)) +
  labs(x = NULL, y = "Follow ratio (participant level)",
       caption = sprintf(paste("Diamond: mean of participant follow ratios (%.3f control,",
                               "%.3f reasoning), bars: 90%% CI\n(participant bootstrap,",
                               "B = 2000). Cell-level rates: %.3f and %.3f, a raw gap of",
                               "%+.1f pp."),
                         arm[condition == "llm-control", m], arm[condition == "llm-reasoning", m],
                         cell_rate[condition == "llm-control", r],
                         cell_rate[condition == "llm-reasoning", r],
                         100 * (cell_rate[condition == "llm-reasoning", r] -
                                cell_rate[condition == "llm-control", r]))) +
  theme_fig() +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0,
                                    margin = margin(t = 3)))

save_fig(p, "fig-primary", w = 5.5, h = 3.2)
cat(sprintf("F1 ok | reasoning %.4f [%.3f, %.3f] | control %.4f [%.3f, %.3f]\n",
            arm[condition == "llm-reasoning", m], arm[condition == "llm-reasoning", lo],
            arm[condition == "llm-reasoning", hi], arm[condition == "llm-control", m],
            arm[condition == "llm-control", lo], arm[condition == "llm-control", hi]))
