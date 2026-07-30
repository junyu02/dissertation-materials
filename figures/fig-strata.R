# F5 -- behaviour against perception, across the engagement strata.
# Claim: as engagement rises the perceived "it explained itself" reading climbs,
#        while the follow rate does not.
#
# The strata are post-allocation self-selection: randomisation does not hold
# within them, so nothing here is a causal read.
source("/tmp/figures-build/_common.R")

psum <- load_psum()
stopifnot(nrow(psum) == 36L)

psum[, grp := fifelse(condition == "llm-control", "control",
              fifelse(ever_opened == 1L, "opened", "never"))]
LEV  <- c("control", "never", "opened")
GLAB <- c(control = "Control", never = "Reasoning\n0 re-opens", opened = "Reasoning\n\u22651 re-open")
GCOL <- c(control = C_CONTROL, never = "#7FB8DB", opened = C_REASON)
GFIL <- c(control = C_CONTROL, never = "white",   opened = C_REASON)
psum[, grp := factor(grp, levels = LEV)]

# --- panel A: follow rate, cell level, participant cluster bootstrap ---------
fr <- rbindlist(lapply(LEV, function(g) {
  d  <- psum[grp == g]
  ci <- boot_g(d, function(x) sum(x$n_follow) / sum(x$n_cells))
  data.table(grp = g, np = nrow(d), nc = sum(d$n_cells),
             est = sum(d$n_follow) / sum(d$n_cells), lo = ci[["lo"]], hi = ci[["hi"]])
}))
stopifnot(near(fr$est, c(0.693, 0.795, 0.758), 1e-3),
          fr$np == c(18L, 9L, 9L))

# --- panel B: MC1 (H1 perceived-reasoning awareness, 1-7) -------------------
mc <- rbindlist(lapply(LEV, function(g) {
  d  <- psum[grp == g & !is.na(MC1)]
  ci <- boot_g(d, function(x) mean(x$MC1))
  data.table(grp = g, n = nrow(d), est = mean(d$MC1), lo = ci[["lo"]], hi = ci[["hi"]])
}))
stopifnot(near(mc$est, c(2.444, 3.875, 6.222), 1e-3),
          mc$n == c(18L, 8L, 9L))

fr[, grp := factor(grp, levels = LEV)]; mc[, grp := factor(grp, levels = LEV)]

base <- function(p) p +
  scale_colour_manual(values = GCOL, guide = "none") +
  scale_fill_manual(values = GFIL, guide = "none") +
  scale_x_discrete(labels = GLAB[LEV]) +
  theme_fig() +
  theme(legend.position = "none", axis.text.x = element_text(size = 7.5),
        axis.title.y = element_text(size = 8.5))

pL <- base(ggplot(fr, aes(grp, est, colour = grp, fill = grp)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .11, linewidth = .55) +
  geom_point(size = 3, shape = 21, stroke = .9) +
  geom_text(aes(y = hi, label = sprintf("%.3f", est)), vjust = -0.9, size = 2.5,
            colour = "grey20") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.55, 0.95),
                     breaks = seq(.6, .9, .1)) +
  labs(x = NULL, y = "Follow rate (cell level, 90% CI)", tag = "A"))

pR <- base(ggplot(mc, aes(grp, est, colour = grp, fill = grp)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .11, linewidth = .55) +
  geom_point(size = 3, shape = 21, stroke = .9) +
  geom_text(aes(y = hi, label = sprintf("%.2f", est)), vjust = -0.9, size = 2.5,
            colour = "grey20") +
  scale_y_continuous(limits = c(1, 7.6), breaks = 1:7) +
  labs(x = NULL, y = "MC1 perceived reasoning, 1-7 (90% CI)", tag = "B"))

# group n's on the axis of each panel, since the two panels have different n
pL <- pL + annotate("text", x = 1:3, y = 0.565, size = 2.3, colour = "grey40",
                    label = sprintf("n = %d", fr$np))
pR <- pR + annotate("text", x = 1:3, y = 1.15, size = 2.3, colour = "grey40",
                    label = sprintf("n = %d", mc$n))

save_fig2(pL, pR, "fig-strata", w = 5.5, h = 3.1, widths = c(1, 1))
cat("F5 ok\n"); print(fr); print(mc)
