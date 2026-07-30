# F6 -- capital allocation share (CAS) on follow cells, by arm.
# Claim: descriptively, the share of available capital committed on a follow cell
#        is higher in the reasoning arm.
source("/tmp/figures-build/_common.R")

cells <- load_cells()
cas <- cells[decision_type == "follow" & !is.na(cas)]

s <- cas[, .(n = .N, np = uniqueN(participant_id), mean = mean(cas), median = median(cas)),
         by = condition]
stopifnot(s[condition == "llm-reasoning", n] == 108L,
          s[condition == "llm-control",   n] ==  91L,
          near(s[condition == "llm-reasoning", mean],   0.505, 1e-3),
          near(s[condition == "llm-control",   mean],   0.315, 1e-3),
          near(s[condition == "llm-reasoning", median], 0.441, 1e-3),
          near(s[condition == "llm-control",   median], 0.229, 1e-3))

pmean <- cas[, .(m = mean(cas)), by = .(participant_id, condition)]

cas[, xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]
pmean[, xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]
s[, xarm := factor(condition, levels = c("llm-control", "llm-reasoning"))]

# participant means sit just to the right of each arm's box; the offset is
# precomputed because ggplot2 cannot compose jitter with nudge in one position.
set.seed(BOOT_SEED)
pmean[, xnum := as.numeric(xarm) + 0.25 + runif(.N, -0.055, 0.055)]

p <- ggplot(cas, aes(xarm, cas, colour = condition, fill = condition)) +
  geom_violin(width = .82, alpha = .16, linewidth = .35, colour = NA, trim = TRUE) +
  geom_boxplot(width = .17, outlier.shape = NA, alpha = 0, linewidth = .45,
               fatten = 1.6, position = position_nudge(x = 0)) +
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
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.02),
                     breaks = seq(0, 1, .25), expand = c(0, 0)) +
  labs(x = NULL, y = "Capital allocation share on the cell",
       caption = paste("Diamond: arm mean. Box: median and IQR. Open circles:",
                       "participant means.\nCAS is defined only on cells adopted in the follow direction.")) +
  theme_fig() +
  theme(legend.position = "none",
        plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0,
                                    margin = margin(t = 3)))

save_fig(p, "fig-cas", w = 5.0, h = 3.3)
cat("F6 ok\n"); print(s)
