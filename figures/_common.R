# Shared style + data loaders for the 7 thesis result figures.
# Sourced by every fig-*.R. Read-only over the analysis pipeline.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(grid)
})

D    <- "/Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道"
EXP  <- file.path(D, "exploratory-注册族-20260729")
RAW  <- file.path(D, "data/real-20260727/raw")
OUT  <- "/tmp/figures-build"

BOOT_B <- 2000L; BOOT_SEED <- 20260727L; CONF <- 0.90

# Okabe-Ito, colour-blind safe
C_REASON <- "#0072B2"; C_CONTROL <- "#E69C00"; C_THIRD <- "#009E73"
ARM_COL  <- c("llm-reasoning" = C_REASON, "llm-control" = C_CONTROL)
ARM_LAB  <- c("llm-reasoning" = "Reasoning", "llm-control" = "Control")

theme_fig <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.line = element_line(linewidth = .3, colour = "grey30"),
          axis.ticks = element_line(linewidth = .3, colour = "grey30"),
          legend.position = "top", legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
          plot.tag = element_text(face = "bold", size = base + 1),
          plot.margin = margin(4, 6, 4, 4))
}

# ---- data -------------------------------------------------------------------
load_cells <- function()
  fread(file.path(D, "data/real-20260727/merged_cells_screened.csv"),
        na.strings = "", encoding = "UTF-8")
load_psum <- function()
  fread(file.path(EXP, "participant_summary.csv"), na.strings = "", encoding = "UTF-8")
load_cwc <- function()
  fread(file.path(EXP, "cells_with_complexity.csv"), na.strings = "", encoding = "UTF-8")
load_blocks <- function()
  fread(file.path(EXP, "reasoning_blocks.csv"), na.strings = "", encoding = "UTF-8")

# dwell: collapse events carrying dwell_ms, restricted to the 36-participant set
load_dwell <- function(pids) {
  revt <- fread(file.path(RAW, "ReasoningEvent.csv"), na.strings = "", encoding = "UTF-8")
  revt <- revt[user %in% pids]
  revt[, dwell_ms := as.numeric(dwell_ms)]
  revt[event_type == "collapse" & is.finite(dwell_ms), .(user, month, block_id, dwell_ms)]
}

# ---- bootstrap --------------------------------------------------------------
# Verbatim port of explore_family.R::boot_p -- participant cluster bootstrap,
# resampled within arm, so the exploratory-family intervals reproduce exactly.
boot_p <- function(psum, statfun, B = BOOT_B, seed = BOOT_SEED) {
  set.seed(seed)
  idx_by_arm <- split(seq_len(nrow(psum)), psum$condition)
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    take <- unlist(lapply(idx_by_arm, function(ix) sample(ix, length(ix), replace = TRUE)),
                   use.names = FALSE)
    out[b] <- tryCatch(statfun(psum[take]), error = function(e) NA_real_)
  }
  ok <- out[is.finite(out)]
  c(lo = unname(quantile(ok, (1 - CONF) / 2)), hi = unname(quantile(ok, 1 - (1 - CONF) / 2)))
}

# Single-group participant bootstrap (each group resampled on its own).
boot_g <- function(d, statfun, B = BOOT_B, seed = BOOT_SEED) {
  set.seed(seed)
  out <- replicate(B, tryCatch(statfun(d[sample.int(nrow(d), nrow(d), TRUE)]),
                               error = function(e) NA_real_))
  ok <- out[is.finite(out)]
  c(lo = unname(quantile(ok, (1 - CONF) / 2)), hi = unname(quantile(ok, 1 - (1 - CONF) / 2)))
}

near <- function(x, y, tol = 5e-4) all(abs(x - y) < tol)

# ---- output -----------------------------------------------------------------
save_fig <- function(p, name, w = 5.5, h = 3.4) {
  # cairo_pdf: keeps the vector output while handling non-ASCII glyphs (GBP sign,
  # en dash) that the base pdf() encoding drops.
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, type = "cairo")
}

# Two panels side by side. patchwork/cowplot are not installed in this R; grid
# viewports do the same job with no new dependency.
save_fig2 <- function(pL, pR, name, w = 5.5, h = 3.2, widths = c(1, 1)) {
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(widths, "null"))))
    print(pL, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pR, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    popViewport()
  }
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h)
  draw(); invisible(dev.off())
  png(file.path(OUT, paste0(name, ".png")), width = w * 300, height = h * 300, res = 300,
      type = "cairo")
  draw(); invisible(dev.off())
}
