#!/usr/bin/env Rscript
# Supervisor verification archive for the manuscript sentence
# "Welch and participant-bootstrap checks left both conclusions unchanged"
# (S2 = all cells; S3 = excluding Month 1 cells). Fixed seed; no new
# manuscript numbers are introduced -- this archives the crossing-zero status only.

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
data_dir <- if (length(args) >= 1) args[[1]] else file.path(here, "data", "real-20260727")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(here, "sensitivity-outputs-20260805")
input <- file.path(data_dir, "merged_cells_screened.csv")
stopifnot(!grepl("synth", input, ignore.case = TRUE), file.exists(input))

d <- fread(input, na.strings = "")
stopifnot(nrow(d) == 387L, uniqueN(d$participant_id) == 36L)
d[, followed := as.integer(followed)]
by_cond <- d[, .(n_participants = uniqueN(participant_id), n_follow = sum(followed), n_cells = .N), by = condition]
stopifnot(all(by_cond$n_participants == 18L),
          by_cond[condition == "llm-reasoning", n_follow] == 161L,
          by_cond[condition == "llm-control", n_follow] == 124L)

check_variant <- function(dd, label) {
  rates <- dd[, .(ratio = mean(followed)), by = .(participant_id, condition)]
  r <- rates[condition == "llm-reasoning", ratio]
  c <- rates[condition == "llm-control", ratio]
  stopifnot(length(r) == 18L, length(c) == 18L)
  welch <- t.test(r, c, var.equal = FALSE, conf.level = 0.90)
  set.seed(20260805)
  B <- 2000L
  boot <- replicate(B, mean(sample(r, replace = TRUE)) - mean(sample(c, replace = TRUE)))
  boot_ci <- unname(quantile(boot, c(0.05, 0.95)))
  data.table(
    variant = label,
    welch_ci_low = welch$conf.int[[1L]], welch_ci_high = welch$conf.int[[2L]],
    welch_crosses_zero = welch$conf.int[[1L]] <= 0 && welch$conf.int[[2L]] >= 0,
    boot_ci_low = boot_ci[[1L]], boot_ci_high = boot_ci[[2L]],
    boot_crosses_zero = boot_ci[[1L]] <= 0 && boot_ci[[2L]] >= 0
  )
}

res <- rbindlist(list(
  check_variant(d, "S2_all_cells"),
  check_variant(d[month >= 2], "S3_excluding_month1")
))
res[, conclusion_unchanged := welch_crosses_zero & boot_crosses_zero]
fwrite(res, file.path(out_dir, "wb_check_results.csv"))
status <- if (all(res$conclusion_unchanged)) "PASS" else "FUSED"
writeLines(status, file.path(out_dir, "wb_check_status.txt"))
print(res)
cat("[status] ", status, "\n", sep = "")
if (status != "PASS") stop("Welch/bootstrap conclusion-unchanged claim not reproduced")
