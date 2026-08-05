#!/usr/bin/env Rscript

suppressMessages({
  library(data.table)
  library(lme4)
})

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
data_dir <- if (length(args) >= 1) args[[1]] else file.path(here, "data", "real-20260727")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(here, "sensitivity-outputs-20260805")
input <- file.path(data_dir, "merged_cells_screened.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stdout_path <- file.path(out_dir, "lopo_stdout.log")
sink(stdout_path, split = TRUE)
on.exit(sink(), add = TRUE)

status_path <- file.path(out_dir, "lopo_run_status.txt")
writeLines("RUNNING", status_path)

stopifnot(!grepl("synth", input, ignore.case = TRUE), file.exists(input))
sha_line <- system2("/usr/bin/shasum", c("-a", "256", shQuote(input)), stdout = TRUE)
input_sha256 <- strsplit(sha_line, "[[:space:]]+")[[1]][1]

dt <- fread(input, na.strings = "")
required <- c("participant_id", "condition", "month", "round", "followed", "literacy", "is_HCI")
stopifnot(length(setdiff(required, names(dt))) == 0L,
          nrow(dt) == 387L,
          uniqueN(dt$participant_id) == 36L,
          !anyNA(dt[, ..required]))

dt[, `:=`(
  followed = as.integer(followed),
  literacy = as.integer(literacy),
  is_HCI = as.integer(is_HCI),
  month = as.integer(month),
  round = as.integer(round)
)]
stopifnot(all(dt$followed %in% 0:1), all(dt$literacy %in% 0:3), all(dt$is_HCI %in% 0:1))

arm_gate <- unique(dt[, .(participant_id, condition)])[, .N, by = condition][order(condition)]
cell_gate <- dt[, .(cells = .N, followed = sum(followed)), by = condition][order(condition)]
participant_means <- dt[, .(ratio = mean(followed)), by = .(participant_id, condition)][,
  .(mean_ratio = mean(ratio)), by = condition][order(condition)]

stopifnot(identical(arm_gate$condition, c("llm-control", "llm-reasoning")),
          identical(arm_gate$N, c(18L, 18L)),
          identical(cell_gate$cells, c(179L, 208L)),
          identical(cell_gate$followed, c(124L, 161L)),
          identical(round(participant_means$mean_ratio, 4), c(0.6648, 0.7417)),
          identical(input_sha256, "6397721d0ae7af5f93044b905f71373b5c06bd2d2389e53ae0ff3a6f6ddb00c8"))

fwrite(data.table(
  check = c("input_not_synthetic", "participants", "cells", "control_participants",
            "reasoning_participants", "control_followed_cells", "reasoning_followed_cells",
            "control_mean_ratio_4dp", "reasoning_mean_ratio_4dp"),
  observed = c("TRUE", "36", "387", "18", "18", "124/179", "161/208", "0.6648", "0.7417"),
  target = c("TRUE", "36", "387", "18", "18", "124/179", "161/208", "0.6648", "0.7417"),
  pass = TRUE
), file.path(out_dir, "lopo_data_gate.csv"))
fwrite(data.table(file = basename(input), sha256 = input_sha256),
       file.path(out_dir, "lopo_input_sha256.csv"))

dt[, condition := relevel(factor(condition), ref = "llm-control")]
ids <- sort(unique(dt$participant_id))
z90 <- qnorm(0.95)

fit_once <- function(d, control = NULL) {
  warnings <- character()
  model <- withCallingHandlers(
    if (is.null(control)) {
      glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) +
              (1 | participant_id), data = d, family = binomial, nAGQ = 1)
    } else {
      glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) +
              (1 | participant_id), data = d, family = binomial, nAGQ = 1,
            control = control)
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  term <- "conditionllm-reasoning"
  beta <- fixef(model)[[term]]
  se <- sqrt(diag(vcov(model)))[[term]]
  list(beta = beta, se = se, warnings = warnings)
}

bobyqa_control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
rows <- rbindlist(lapply(seq_along(ids), function(i) {
  d <- droplevels(dt[participant_id != ids[[i]]])
  default <- fit_once(d)
  bobyqa <- fit_once(d, bobyqa_control)
  data.table(
    omitted = sprintf("leaveout_%02d", i),
    n_participants = uniqueN(d$participant_id),
    n_cells = nrow(d),
    or_default = exp(default$beta),
    ci90_low_default = exp(default$beta - z90 * default$se),
    ci90_high_default = exp(default$beta + z90 * default$se),
    default_warning = length(default$warnings) > 0L,
    default_warning_text = paste(default$warnings, collapse = " | "),
    or_bobyqa = exp(bobyqa$beta),
    bobyqa_warning = length(bobyqa$warnings) > 0L,
    bobyqa_warning_text = paste(bobyqa$warnings, collapse = " | "),
    abs_log_or_difference = abs(default$beta - bobyqa$beta)
  )
}))
fwrite(rows, file.path(out_dir, "lopo_results.csv"))

summary <- data.table(
  min_or = min(rows$or_default),
  max_or = max(rows$or_default),
  intervals_excluding_1 = sum(rows$ci90_low_default > 1 | rows$ci90_high_default < 1),
  default_warning_refits = sum(rows$default_warning),
  bobyqa_warning_refits = sum(rows$bobyqa_warning),
  max_abs_log_or_difference = max(rows$abs_log_or_difference)
)
summary[, `:=`(
  claim_or_range = round(min_or, 2) == 1.41 && round(max_or, 2) == 1.99,
  claim_one_interval = intervals_excluding_1 == 1L,
  claim_optimizer = default_warning_refits == 25L && bobyqa_warning_refits == 0L &&
    max_abs_log_or_difference < 0.01
)]
summary[, all_claims_reproduced := claim_or_range && claim_one_interval && claim_optimizer]
fwrite(summary, file.path(out_dir, "lopo_summary.csv"))

cat(sprintf("[gate] participants=36 cells=387 conditions=18:18 sha256=%s\n", input_sha256))
cat(sprintf("[LOPO] OR range %.10f--%.10f; intervals excluding 1: %d/36\n",
            summary$min_or, summary$max_or, summary$intervals_excluding_1))
cat(sprintf("[LOPO] default warning refits=%d; bobyqa warning refits=%d; max |delta log(OR)|=%.10g\n",
            summary$default_warning_refits, summary$bobyqa_warning_refits,
            summary$max_abs_log_or_difference))

capture.output(sessionInfo(), file = file.path(out_dir, "lopo_session_info.txt"))
if (!summary$all_claims_reproduced) {
  writeLines("FUSED", status_path)
  stop("LOPO claim gate failed; see lopo_summary.csv")
}
writeLines("PASS", status_path)
cat("[status] PASS: all three manuscript claims reproduced\n")
