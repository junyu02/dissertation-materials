# F4 -- reasoning-panel engagement: how often panels were voluntarily reopened,
# and how long they stayed open after that event.
# Claim: panels were rarely reopened; dwell time is reported without assigning a
# reading-depth label.
source("/tmp/figures-build/_common.R")

cells <- load_cells(); blocks <- load_blocks()
dwell <- load_dwell(unique(cells$participant_id))

stopifnot(nrow(blocks) == 273L, sum(blocks$opened) == 51L,
          near(sum(blocks$opened) / nrow(blocks), 0.187, 1e-3))

mo <- blocks[, .(n = .N, k = sum(opened), r = sum(opened) / .N), by = month][order(month)]
stopifnot(nrow(mo) == 5L, sum(mo$n) == 273L, sum(mo$k) == 51L)

stopifnot(nrow(dwell) == 52L,
          median(dwell$dwell_ms) == 11896,
          round(quantile(dwell$dwell_ms, .25)) == 5263,
          quantile(dwell$dwell_ms, .75) == 46749)

dwell[, sec := dwell_ms / 1000]
q <- quantile(dwell$sec, c(.25, .5, .75))

pL <- ggplot(mo, aes(factor(month), r)) +
  geom_col(fill = C_REASON, width = .62, alpha = .85) +
  geom_text(aes(label = sprintf("%d/%d", k, n)), vjust = -0.6, size = 2.5, colour = "grey20") +
  geom_hline(yintercept = sum(mo$k) / sum(mo$n), linetype = "22", linewidth = .45,
             colour = "grey30") +
  annotate("text", x = 3.0, y = sum(mo$k) / sum(mo$n) + .018, hjust = .5, size = 2.4,
           colour = "grey30", label = "overall 51 / 273 = 18.7%") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.365),
                     breaks = seq(0, .3, .1), expand = c(0, 0)) +
  labs(x = "Simulated month", y = "Panels voluntarily reopened / rendered", tag = "A") +
  theme_fig() +
  theme(legend.position = "none", axis.title = element_text(size = 8),
        plot.margin = margin(6, 4, 4, 4))

pR <- ggplot(dwell, aes(sec)) +
  annotate("rect", xmin = q[[1]], xmax = q[[3]], ymin = 0, ymax = Inf,
           fill = C_REASON, alpha = .13) +
  geom_histogram(bins = 14, fill = C_REASON, colour = "white", linewidth = .25, alpha = .85) +
  geom_vline(xintercept = q[[2]], linewidth = .5, colour = "grey20") +
  geom_rug(sides = "b", linewidth = .25, colour = "grey35", alpha = .7) +
  annotate("text", x = 700, y = 11.6, hjust = 1, size = 2.5, colour = "grey20",
           label = sprintf("median %.1f s", q[[2]])) +
  annotate("text", x = 700, y = 10.3, hjust = 1, size = 2.4, colour = "grey35",
           label = sprintf("IQR %.1f-%.1f s", q[[1]], q[[3]])) +
  scale_x_log10(breaks = c(1, 10, 100),
                labels = c("1 s", "10 s", "100 s"),
                expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 12, 2),
                     expand = expansion(mult = c(0, .04))) +
  coord_cartesian(xlim = c(0.6, 800)) +
  labs(x = "Dwell per reopen-collapse pair (log scale)",
       y = sprintf("Reopen-collapse pairs (n = %d)", nrow(dwell)), tag = "B") +
  theme_fig() +
  theme(legend.position = "none", axis.title = element_text(size = 8),
        panel.grid.major.x = element_line(linewidth = .2, colour = "grey93"),
        plot.margin = margin(6, 6, 4, 4))

save_fig2(pL, pR, "fig-engagement", w = 5.5, h = 2.9, widths = c(1, 1.05))
cat("F4 ok | monthly re-opened-panel rates:\n"); print(mo)
cat(sprintf("dwell n=%d median=%.0f ms IQR [%.0f, %.0f]\n", nrow(dwell),
            q[[2]] * 1000, q[[1]] * 1000, q[[3]] * 1000))
