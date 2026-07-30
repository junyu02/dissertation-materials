# F2 -- primary estimate and its sensitivity siblings, odds-ratio scale.
# Claim: every interval that adjusts for participant clustering contains OR = 1.
#
# These five rows are LOCKED VALUES supplied with the figure specification, not
# recomputed here (the GLMM / Firth / CR1 fits live in the confirmatory pipeline).
# Source of record: 00-论文交付/分析管道 confirmatory outputs as quoted in the figure
# spec. The cluster-adjusted twin of row 1 in descriptive form is
# exploratory-注册族-20260729/results-exploratory.md, T7 ITT (+0.081 [-0.043, 0.202]).
source("/tmp/figures-build/_common.R")

# y runs top (largest) to bottom; header rows carry no estimate.
rows <- data.table(
  y      = 9:1,
  kind   = c("h", "e", "h", "e", "h", "e", "e", "h", "e"),
  label  = c("Registered primary",
             "GLMM, Wald interval",
             "Registered sensitivity",
             "GLMM, profile interval",
             "Post-hoc sensitivity",
             "Cluster-bootstrap Firth",
             "CR1 cluster-robust, pooled",
             "Registered fallback (strict-screen N = 36)",
             "Firth, collapsed"),
  or = c(NA, 1.68, NA, 1.68, NA, 1.51, 1.53, NA, 1.51),
  lo = c(NA, 0.83, NA, 0.82, NA, 0.82, 0.82, NA, 1.02),
  hi = c(NA, 3.42, NA, 3.60, NA, 2.97, 2.83, NA, 2.24))
rows[, clustered := kind == "e" & y != 1]

est <- rows[kind == "e"]
stopifnot(nrow(est) == 5L,
          est[clustered == TRUE, all(lo < 1 & hi > 1)],
          est[clustered == FALSE, lo > 1 & near(or, 1.51, 1e-9)])

YL <- c(-0.15, 9.6)

# --- left panel: row labels + numeric column ---------------------------------
pL <- ggplot(rows, aes(y = y)) +
  geom_text(data = rows[kind == "h"], aes(x = 0, label = label), hjust = 0,
            size = 2.55, fontface = "bold", colour = "grey30") +
  geom_text(data = est, aes(x = 0.05, label = label), hjust = 0,
            size = 2.75, colour = "grey10") +
  geom_text(data = est, aes(x = 1, label = sprintf("%.2f [%.2f, %.2f]", or, lo, hi)),
            hjust = 1, size = 2.65, colour = "grey10", family = "mono") +
  annotate("text", x = 0.05, y = 0.35, hjust = 0, size = 2.35, fontface = "italic",
           colour = "grey35",
           label = "Open point, dashed bar: different estimand,\nno cluster adjustment.") +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  labs(x = "Odds ratio, reasoning vs control (90% CI)", y = NULL) +
  theme_fig() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), panel.grid = element_blank(),
        axis.title.x = element_text(colour = "white", size = 8),  # reserves the same
        plot.margin = margin(4, 2, 4, 4))               # baseline height as pR

# --- right panel: the intervals ----------------------------------------------
pR <- ggplot(est, aes(or, y)) +
  geom_vline(xintercept = 1, linewidth = .4, colour = "grey35") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, linetype = clustered),
                 height = .18, linewidth = .55, colour = C_REASON) +
  geom_point(aes(shape = clustered), size = 2.3, colour = C_REASON,
             fill = "white", stroke = .8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21), guide = "none") +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "22"), guide = "none") +
  scale_x_log10(breaks = c(0.8, 1, 1.5, 2, 3, 4),
                labels = c("0.8", "1", "1.5", "2", "3", "4"),
                limits = c(0.72, 4.0), expand = c(0, 0)) +
  scale_y_continuous(limits = YL, expand = c(0, 0)) +
  labs(x = "Odds ratio, reasoning vs control (90% CI)", y = NULL) +
  theme_fig() +
  theme(axis.title.x = element_text(size = 8),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = .25, colour = "grey92"),
        plot.margin = margin(4, 6, 4, 2))

save_fig2(pL, pR, "fig-forest", w = 5.5, h = 2.9, widths = c(1.22, 1))
cat("F2 ok | 5 rows; 4 cluster-adjusted intervals contain 1; fallback excludes 1\n")
