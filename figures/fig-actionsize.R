# F7 -- absolute action size, primary denominator (every cell carrying a positive
# capital move), by arm.
# Claim: descriptively, positive capital moves are larger in the reasoning arm.
source("/tmp/figures-build/_common.R")

cells <- load_cells()
mv <- cells[!is.na(invested) & invested > 0]

s <- mv[, .(n = .N, np = uniqueN(participant_id), mean = mean(invested),
            median = median(invested)), by = condition]
stopifnot(s[condition == "llm-reasoning", n] == 161L,
          s[condition == "llm-control",   n] == 124L,
          near(s[condition == "llm-reasoning", mean],   460.660, 1e-2),
          near(s[condition == "llm-control",   mean],   333.077, 1e-2),
          near(s[condition == "llm-reasoning", median], 300.000, 1e-2),
          near(s[condition == "llm-control",   median], 277.330, 1e-2))

pmean <- mv[, .(m = mean(invested)), by = .(participant_id, condition)]

mv[,    xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]
pmean[, xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]
s[,     xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]

# participant means sit just to the right of each arm's box; the offset is
# precomputed because ggplot2 cannot compose jitter with nudge in one position.
set.seed(BOOT_SEED)
pmean[, xnum := as.numeric(xarm) + 0.25 + runif(.N, -0.055, 0.055)]

p <- ggplot(mv, aes(xarm, invested, colour = condition, fill = condition)) +
  geom_violin(width = .82, alpha = .16, linewidth = .35, colour = NA, trim = TRUE) +
  geom_boxplot(width = .17, outlier.shape = NA, alpha = 0, linewidth = .45, fatten = 1.6) +
  geom_point(data = pmean, aes(x = xnum, y = m), shape = 21, size = 1.8, stroke = .7,
             fill = "white") +
  geom_point(data = s, aes(y = mean), shape = 18, size = 3.2) +
  scale_colour_manual(values = ARM_COL, guide = "none") +
  scale_fill_manual(values = ARM_COL, guide = "none") +
  scale_x_discrete(labels = c("llm-control" = sprintf("Control\n%d cells, %d participants",
                                                      s[condition == "llm-control", n],
                                                      s[condition == "llm-control", np]),
                              "llm-reasoning" = sprintf("Reasoning\n%d cells, %d participants",
                                                      s[condition == "llm-reasoning", n],
                                                      s[condition == "llm-reasoning", np])),
                   expand = expansion(add = 0.62)) +
  scale_y_log10(breaks = c(10, 50, 100, 250, 500, 1000, 3000),
                labels = function(x) paste0("\u00a3", comma(x)),
                limits = c(6, 5200), expand = c(0, 0)) +
  labs(x = NULL, y = "Size of the capital move (log scale)",
       caption = paste("Diamond: arm mean. Box: median and IQR. Open circles:",
                       "participant means.\nDenominator: every executed capital move",
                       "with amount > \u00a30, including adopted step-back moves.")) +
  theme_fig() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(linewidth = .2, colour = "grey93"),
        plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0,
                                    margin = margin(t = 3)))

save_fig(p, "fig-actionsize", w = 5.0, h = 3.3)
cat("F7 ok\n"); print(s)
