#!/usr/bin/env Rscript
# Negative tests for fit_glmm.R's integrity gate (R2 cross-check finding: the
# alarm branches were never exercised — a gate no test has ever tripped is a
# gate you only THINK you have).
#
# Each case takes the clean synthetic csv, corrupts exactly one thing, runs
# fit_glmm.R on it and asserts the pipeline STOPS with the expected message.
# A final control case asserts the clean file still passes.
#
#   Rscript test_alarms.R          # uses data/synth_cells.csv (run synth_data.py first)

suppressMessages(library(data.table))

here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 error = function(e) ".")
if (length(here) == 0 || here == "") here <- "."
base_csv <- file.path(here, "data", "synth_cells.csv")
fit_r    <- file.path(here, "fit_glmm.R")
stopifnot(file.exists(base_csv), file.exists(fit_r))

base <- fread(base_csv, na.strings = "")
tmpdir <- tempdir()

run_case <- function(name, mutate, expect_pattern, expect_fail = TRUE) {
  dt <- copy(base)
  dt <- mutate(dt)
  in_csv  <- file.path(tmpdir, paste0("case_", gsub("[^a-z0-9]+", "_", tolower(name)), ".csv"))
  out_tex <- file.path(tmpdir, "out.tex")
  fwrite(dt, in_csv)
  res <- suppressWarnings(system2("Rscript", c(fit_r, in_csv, out_tex),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status"); if (is.null(status)) status <- 0L
  txt <- paste(res, collapse = "\n")
  ok_fail <- if (expect_fail) status != 0L else status == 0L
  ok_msg  <- grepl(expect_pattern, txt, fixed = FALSE)
  if (ok_fail && ok_msg) {
    cat(sprintf("  PASS  %-38s (exit=%d)\n", name, status))
    TRUE
  } else {
    cat(sprintf("  FAIL  %-38s exit=%d msg-match=%s\n", name, status, ok_msg))
    cat(sprintf("        expected pattern: %s\n", expect_pattern))
    FALSE
  }
}

# pick reference rows of each kind (guaranteed by synth_data's self-check)
i_follow   <- which(base$decision_type == "follow")[1]
i_unfollow <- which(base$decision_type == "unfollow")[1]
i_decline  <- which(base$decision_type == "decline")[1]

cat("=== fit_glmm.R integrity-gate negative tests ===\n")
results <- c(
  run_case("NA decision_type slips nothing",
           function(d) { d$decision_type[i_decline] <- NA; d },
           "decision_type"),
  run_case("unknown decision_type enum",
           function(d) { d$decision_type[i_decline] <- "timeout"; d },
           "decision_type"),
  run_case("NA direction",
           function(d) { d$direction[i_unfollow] <- NA; d },
           "direction"),
  run_case("follow row with NULL cas (true missing)",
           function(d) { d$cas[i_follow] <- NA; d },
           "true missing"),
  run_case("unfollow row with a cas value",
           function(d) { d$cas[i_unfollow] <- 0.5; d },
           "non-follow rows with non-NULL cas"),
  run_case("followed=1 on a decline",
           function(d) { d$followed[i_decline] <- 1; d },
           "disagrees with decision_type"),
  run_case("unfollow with invested 0",
           function(d) { d$invested[i_unfollow] <- 0; d },
           "unfollow rows without a positive invested"),
  run_case("follow with NA invested",
           function(d) { d$invested[i_follow] <- NA; d },
           "follow rows without a positive invested|cas !="),
  run_case("duplicated (participant,month,round)",
           function(d) rbind(d, d[i_follow]),
           "duplicate \\(participant_id, month, round\\)"),
  run_case("fractional followed 1.7 not truncated",
           function(d) { d$followed[i_follow] <- 1.7; d },
           "followed"),
  run_case("negative available_snapshot",
           function(d) { d$available_snapshot[i_decline] <- -5; d },
           "available_snapshot"),
  run_case("zero_optout on an unfollow-direction cell",
           function(d) { d$decision_type[i_unfollow] <- "zero_optout"
                         d$followed[i_unfollow] <- 0
                         d$invested[i_unfollow] <- 0
                         d$cas[i_unfollow] <- NA; d },
           "contract-undefined combination"),
  run_case("clean file still passes (control)",
           function(d) d,
           "=== results ===", expect_fail = FALSE)
)

cat(sprintf("=== %d/%d cases passed ===\n", sum(results), length(results)))
if (!all(results)) quit(status = 1)
