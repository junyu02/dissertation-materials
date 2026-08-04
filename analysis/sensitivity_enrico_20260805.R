#!/usr/bin/env Rscript
# Enrico sensitivity analyses, frozen before the single production run.
#
# S1a/S1b: population-averaged binary GEE, logit link, participant clusters,
# exchangeable working correlation, Mancl-DeRouen small-sample robust covariance,
# two-sided 90% Wald CI. S1a has condition only; S1b adds financial literacy,
# is_HCI, and month/round fixed dummies.
# S1c: condition-only binomial GLMM with participant random intercept, Laplace
# approximation, registered convergence gate, OR + Wald 90% CI + random-intercept SD.
# S2: participant Follow Ratio, equal-variance two-sample t test, df=34, 90% CI.
# S3: S2 after excluding month 1 cells; every participant must retain >=1 cell.
# S3b: registered adjusted GLMM restricted to months 2-5; report the registered
# convergence gate as observed and do not change optimizer/specification.
#
# Fuse rules: target rounding mismatch, a confidence interval excluding its null,
# or any data/S3 assertion failure stops the run after complete outputs are written.

suppressMessages({
  library(data.table)
  library(lme4)
})

args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!length(here) || !nzchar(here)) here <- "."
data_dir <- if (length(args) >= 1) args[[1]] else file.path(here, "data", "real-20260727")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(here, "sensitivity-outputs-20260805")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Z90 <- qnorm(0.95)
GRAD_TOL <- 0.002
cells_path <- file.path(data_dir, "merged_cells_screened.csv")
public_cells_path <- file.path(data_dir, "study_cells_screened.csv")
input_paths <- c(cells_path, public_cells_path)
stopifnot(!any(grepl("synth", input_paths, ignore.case = TRUE)), all(file.exists(input_paths)))

sha256 <- function(path) {
  line <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (length(line) != 1L) stop("SHA-256 failed for ", path)
  strsplit(line, "[[:space:]]+")[[1L]][1L]
}

hashes <- data.table(
  file = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  sha256 = vapply(input_paths, sha256, character(1))
)
fwrite(hashes, file.path(out_dir, "sensitivity_input_sha256.csv"))

gates <- data.table(assertion = character(), observed = character(), expected = character(), passed = logical())
gate <- function(name, observed, expected, passed) {
  gates <<- rbind(gates, data.table(assertion = name, observed = as.character(observed),
                                    expected = as.character(expected), passed = isTRUE(passed)))
  if (!isTRUE(passed)) stop(sprintf("DATA GATE FAILED: %s; observed=%s expected=%s", name, observed, expected))
}

d <- fread(cells_path, na.strings = "")
required <- c("participant_id", "condition", "month", "round", "followed", "literacy", "is_HCI")
gate("required columns", paste(setdiff(required, names(d)), collapse = ","), "none missing",
     length(setdiff(required, names(d))) == 0L)
d[, `:=`(month = as.integer(month), round = as.integer(round), followed = as.integer(followed),
         literacy = as.integer(literacy), is_HCI = as.integer(is_HCI))]
gate("participant N", uniqueN(d$participant_id), 36L, uniqueN(d$participant_id) == 36L)
gate("decision cells", nrow(d), 387L, nrow(d) == 387L)
by_condition <- d[, .(n_participants = uniqueN(participant_id), n_cells = .N,
                      n_follow = sum(followed)), by = condition]
gate("condition levels", paste(sort(by_condition$condition), collapse = ","),
     "llm-control,llm-reasoning",
     setequal(by_condition$condition, c("llm-control", "llm-reasoning")))
gate("participants per condition", paste(sort(by_condition$n_participants), collapse = ","),
     "18,18", all(by_condition$n_participants == 18L))
gate("reasoning cells/follows",
     paste(by_condition[condition == "llm-reasoning", n_follow],
           by_condition[condition == "llm-reasoning", n_cells], sep = "/"),
     "161/208",
     identical(by_condition[condition == "llm-reasoning", c(n_follow, n_cells)], c(161L, 208L)))
gate("control cells/follows",
     paste(by_condition[condition == "llm-control", n_follow],
           by_condition[condition == "llm-control", n_cells], sep = "/"),
     "124/179",
     identical(by_condition[condition == "llm-control", c(n_follow, n_cells)], c(124L, 179L)))
participant_rates <- d[, .(follow_ratio = mean(followed), n_cells = .N),
                       by = .(participant_id, condition)]
mean_reasoning <- participant_rates[condition == "llm-reasoning", mean(follow_ratio)]
mean_control <- participant_rates[condition == "llm-control", mean(follow_ratio)]
gate("participant mean reasoning", sprintf("%.10f", mean_reasoning), "round4=0.7417",
     round(mean_reasoning, 4) == 0.7417)
gate("participant mean control", sprintf("%.10f", mean_control), "round4=0.6648",
     round(mean_control, 4) == 0.6648)
gate("public cells SHA-256", hashes[file == normalizePath(public_cells_path), sha256],
     "934ad28634a1f379a5cf269d2a333efcd5089c50b3c9195fec8e8250539a7c8f",
     hashes[file == normalizePath(public_cells_path), sha256] ==
       "934ad28634a1f379a5cf269d2a333efcd5089c50b3c9195fec8e8250539a7c8f")

d[, condition := relevel(factor(condition), ref = "llm-control")]
d[, participant_id := factor(participant_id)]
TERM <- "conditionllm-reasoning"

fit_exchangeable_gee <- function(formula, data, tol = 1e-10, max_iter = 100L) {
  X <- model.matrix(formula, data)
  y <- data$followed
  clusters <- split(seq_along(y), data$participant_id)
  beta <- coef(glm.fit(X, y, family = binomial()))
  p <- ncol(X)

  build_state <- function(beta) {
    mu <- plogis(drop(X %*% beta))
    pearson <- (y - mu) / sqrt(mu * (1 - mu))
    numerator <- sum(vapply(clusters, function(ii) {
      products <- outer(pearson[ii], pearson[ii])
      sum(products[upper.tri(products)])
    }, numeric(1)))
    denominator <- sum(vapply(clusters, function(ii) choose(length(ii), 2), numeric(1))) - p
    alpha <- numerator / denominator
    lower <- -1 / (max(lengths(clusters)) - 1) + 1e-6
    alpha <- min(0.99, max(lower, alpha))
    bread <- matrix(0, p, p)
    score <- numeric(p)
    pieces <- vector("list", length(clusters))
    for (k in seq_along(clusters)) {
      ii <- clusters[[k]]
      mui <- mu[ii]
      ni <- length(ii)
      A <- diag(mui * (1 - mui), nrow = ni)
      Ahalf <- diag(sqrt(mui * (1 - mui)), nrow = ni)
      R <- matrix(alpha, ni, ni)
      diag(R) <- 1
      V <- Ahalf %*% R %*% Ahalf
      W <- solve(V)
      D <- A %*% X[ii, , drop = FALSE]
      bread <- bread + t(D) %*% W %*% D
      score <- score + drop(t(D) %*% W %*% (y[ii] - mui))
      pieces[[k]] <- list(ii = ii, D = D, W = W)
    }
    list(mu = mu, alpha = alpha, bread = bread, score = score, pieces = pieces)
  }

  converged <- FALSE
  for (iteration in seq_len(max_iter)) {
    state <- build_state(beta)
    beta_new <- beta + solve(state$bread, state$score)
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }
  if (!converged) stop("GEE did not converge in ", max_iter, " iterations")

  state <- build_state(beta)
  bread_inv <- solve(state$bread)
  meat_md <- matrix(0, p, p)
  for (piece in state$pieces) {
    ii <- piece$ii
    H <- piece$D %*% bread_inv %*% t(piece$D) %*% piece$W
    residual_md <- solve(diag(length(ii)) - H, y[ii] - state$mu[ii])
    cluster_score <- drop(t(piece$D) %*% piece$W %*% residual_md)
    meat_md <- meat_md + tcrossprod(cluster_score)
  }
  vcov_md <- bread_inv %*% meat_md %*% bread_inv
  list(coefficients = beta, vcov = vcov_md, alpha = state$alpha,
       iterations = iteration, n_cluster = length(clusters), n_parameters = p)
}

gee_row <- function(id, formula) {
  fit <- fit_exchangeable_gee(formula, d)
  beta <- unname(fit$coefficients[[TERM]])
  se <- sqrt(diag(fit$vcov))[[TERM]]
  data.table(
    analysis = id, scale = "odds ratio", model = paste(deparse(formula, width.cutoff = 500L), collapse = " "),
    estimate = exp(beta), ci_low = exp(beta - Z90 * se), ci_high = exp(beta + Z90 * se),
    p_value = 2 * pnorm(-abs(beta / se)), standard_error_log = se,
    tau = NA_real_, working_alpha = fit$alpha, n_participants = fit$n_cluster,
    n_cells = nrow(d), df = NA_real_, converged = TRUE,
    convergence_detail = sprintf("GEE iterations=%d; MD robust covariance", fit$iterations)
  )
}

glmm_row <- function(id, data, formula) {
  captured_warnings <- character()
  m <- withCallingHandlers(
    glmer(formula, data = data, family = binomial, nAGQ = 1),
    warning = function(w) captured_warnings <<- c(captured_warnings, conditionMessage(w))
  )
  conv_msgs <- m@optinfo$conv$lme4$messages
  if (is.null(conv_msgs)) conv_msgs <- character()
  opt_code <- m@optinfo$conv$opt
  derivs <- m@optinfo$derivs
  scaled_grad <- tryCatch(max(abs(solve(derivs$Hessian, derivs$gradient))),
                          error = function(e) NA_real_)
  converged <- length(conv_msgs) == 0L && is.finite(scaled_grad) &&
    scaled_grad <= GRAD_TOL && isTRUE(opt_code == 0)
  beta <- unname(fixef(m)[[TERM]])
  se <- sqrt(diag(vcov(m)))[[TERM]]
  tau <- unname(attr(VarCorr(m)$participant_id, "stddev")[[1L]])
  detail <- paste(unique(c(sprintf("scaled_grad=%.12g", scaled_grad),
                           sprintf("optimizer_code=%s", paste(opt_code, collapse = ",")),
                           sprintf("singular=%s", isSingular(m)), conv_msgs, captured_warnings)),
                  collapse = " | ")
  data.table(
    analysis = id, scale = "odds ratio", model = paste(deparse(formula, width.cutoff = 500L), collapse = " "),
    estimate = exp(beta), ci_low = exp(beta - Z90 * se), ci_high = exp(beta + Z90 * se),
    p_value = 2 * pnorm(-abs(beta / se)), standard_error_log = se,
    tau = tau, working_alpha = NA_real_, n_participants = uniqueN(data$participant_id),
    n_cells = nrow(data), df = NA_real_, converged = converged,
    convergence_detail = detail
  )
}

t_row <- function(id, data) {
  rates <- data[, .(follow_ratio = mean(followed), n_cells = .N),
                by = .(participant_id, condition)]
  reasoning <- rates[condition == "llm-reasoning", follow_ratio]
  control <- rates[condition == "llm-control", follow_ratio]
  fit <- t.test(reasoning, control, var.equal = TRUE, conf.level = 0.90)
  data.table(
    analysis = id, scale = "proportion difference", model = "equal-variance two-sample t",
    estimate = mean(reasoning) - mean(control), ci_low = unname(fit$conf.int[[1L]]),
    ci_high = unname(fit$conf.int[[2L]]), p_value = fit$p.value,
    standard_error_log = NA_real_, tau = NA_real_, working_alpha = NA_real_,
    n_participants = length(reasoning) + length(control), n_cells = nrow(data),
    df = unname(fit$parameter), converged = TRUE, convergence_detail = "not applicable"
  )
}

s1a <- gee_row("S1a", ~ condition)
s1b <- gee_row("S1b", ~ condition + literacy + is_HCI + factor(month) + factor(round))
s1c <- glmm_row("S1c", d, followed ~ condition + (1 | participant_id))
s2 <- t_row("S2", d)

d_m2_5 <- d[month >= 2]
retained <- d_m2_5[, .N, by = .(participant_id, condition)]
gate("S3 participants retained", nrow(retained), 36L, nrow(retained) == 36L)
gate("S3 positive denominator", min(retained$N), ">0", all(retained$N > 0L))
gate("S3 conditions retained", paste(retained[, .N, by = condition]$N, collapse = ","),
     "18,18", all(retained[, .N, by = condition]$N == 18L))
s3 <- t_row("S3", d_m2_5)
s3b <- glmm_row("S3b", d_m2_5,
                followed ~ condition + literacy + is_HCI + factor(month) + factor(round) +
                  (1 | participant_id))

results <- rbindlist(list(s1a, s1b, s1c, s2, s3, s3b), fill = TRUE)
results[, `:=`(
  estimate_rounded = ifelse(scale == "odds ratio", sprintf("%.2f", estimate),
                            sprintf("%.1f pp", 100 * estimate)),
  ci_rounded = ifelse(scale == "odds ratio",
                      sprintf("[%.2f, %.2f]", ci_low, ci_high),
                      sprintf("[%.1f, %.1f] pp", 100 * ci_low, 100 * ci_high)),
  p_rounded = sprintf("%.3f", p_value),
  interval_includes_null = ifelse(scale == "odds ratio",
                                  ci_low <= 1 & ci_high >= 1,
                                  ci_low <= 0 & ci_high >= 0)
)]

descriptives <- rbindlist(list(
  participant_rates[, .(n_participants = .N, mean_follow_ratio = mean(follow_ratio),
                        sd_follow_ratio = sd(follow_ratio), min_cells = min(n_cells), max_cells = max(n_cells)),
                    by = condition][, analysis := "S2"],
  d_m2_5[, .(follow_ratio = mean(followed), n_cells = .N), by = .(participant_id, condition)][,
    .(n_participants = .N, mean_follow_ratio = mean(follow_ratio), sd_follow_ratio = sd(follow_ratio),
      min_cells = min(n_cells), max_cells = max(n_cells)), by = condition][, analysis := "S3"]
))
setcolorder(descriptives, c("analysis", "condition", setdiff(names(descriptives), c("analysis", "condition"))))

diagnostics <- rbindlist(lapply(c("S2", "S3"), function(id) {
  x <- if (id == "S2") participant_rates else
    d_m2_5[, .(follow_ratio = mean(followed)), by = .(participant_id, condition)]
  r <- x[condition == "llm-reasoning", follow_ratio]
  c <- x[condition == "llm-control", follow_ratio]
  data.table(analysis = id,
             shapiro_p_reasoning = shapiro.test(r)$p.value,
             shapiro_p_control = shapiro.test(c)$p.value,
             variance_ratio_reasoning_over_control = var(r) / var(c),
             variance_test_p = var.test(r, c)$p.value)
}))

target_checks <- data.table(
  item = c("S1a OR", "S1b OR", "S2 reasoning mean", "S2 control mean", "S2 difference",
           "S2 CI low", "S2 CI high", "S2 p", "S3 reasoning mean", "S3 control mean",
           "S3 difference", "S3 CI low", "S3 CI high", "S3 p"),
  observed = c(round(s1a$estimate, 2), round(s1b$estimate, 2),
               round(100 * descriptives[analysis == "S2" & condition == "llm-reasoning", mean_follow_ratio], 1),
               round(100 * descriptives[analysis == "S2" & condition == "llm-control", mean_follow_ratio], 1),
               round(100 * s2$estimate, 1), round(100 * s2$ci_low, 1), round(100 * s2$ci_high, 1),
               round(s2$p_value, 3),
               round(100 * descriptives[analysis == "S3" & condition == "llm-reasoning", mean_follow_ratio], 1),
               round(100 * descriptives[analysis == "S3" & condition == "llm-control", mean_follow_ratio], 1),
               round(100 * s3$estimate, 1), round(100 * s3$ci_low, 1), round(100 * s3$ci_high, 1),
               round(s3$p_value, 3)),
  target = c(1.48, 1.50, 74.2, 66.5, 7.7, -7.0, 22.4, 0.382,
             75.7, 67.6, 8.1, -8.0, 24.1, 0.401)
)
target_checks[, matched := observed == target]

fwrite(gates, file.path(out_dir, "data_gate_assertions.csv"))
fwrite(results, file.path(out_dir, "sensitivity_results.csv"))
fwrite(descriptives, file.path(out_dir, "participant_descriptives.csv"))
fwrite(diagnostics, file.path(out_dir, "assumption_diagnostics.csv"))
fwrite(target_checks, file.path(out_dir, "target_comparison.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

fuse_reasons <- character()
if (any(!target_checks$matched))
  fuse_reasons <- c(fuse_reasons, paste("target mismatch:", paste(target_checks[matched == FALSE, item], collapse = ", ")))
if (any(!results$interval_includes_null))
  fuse_reasons <- c(fuse_reasons, paste("90% interval excludes null:",
                                        paste(results[interval_includes_null == FALSE, analysis], collapse = ", ")))
status <- if (length(fuse_reasons)) c("FUSED", fuse_reasons) else "PASS"
writeLines(status, file.path(out_dir, "sensitivity_run_status.txt"))

cat("[data gates]\n")
print(gates)
cat("\n[results: unrounded and rounded]\n")
print(results)
cat("\n[target comparison]\n")
print(target_checks)
cat("\n[status] ", status[[1L]], "\n", sep = "")
if (length(fuse_reasons)) stop(paste(fuse_reasons, collapse = "; "))
