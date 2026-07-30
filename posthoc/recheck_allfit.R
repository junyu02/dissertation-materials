#!/usr/bin/env Rscript
# Read-only recheck: allFit on N=38 and N=36; subsample stability check.
suppressMessages({ library(data.table); library(lme4) })

DIR <- "/Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道/data/real-20260727"
GRAD_TOL <- 0.002
TERM <- "conditionllm-reasoning"

prep <- function(f) {
  dt <- fread(file.path(DIR, f), na.strings = "")
  dt[, round := as.integer(round)]
  dt[, condition := relevel(factor(condition), ref = "llm-control")]
  dt[, participant_id := factor(participant_id)]
  dt[]
}

fit_one <- function(dt) {
  glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) +
          (1 | participant_id), data = dt, family = binomial, nAGQ = 1)
}

sg <- function(m) {
  dd <- m@optinfo$derivs
  tryCatch(max(abs(solve(dd$Hessian, dd$gradient))), error = function(e) NA_real_)
}

report <- function(f, label) {
  cat("\n\n################ ", label, " (", f, ") ################\n", sep = "")
  dt <- prep(f)
  cat(sprintf("rows=%d  participants=%d\n", nrow(dt), nlevels(dt$participant_id)))
  m <- fit_one(dt)
  cat(sprintf("base fit: OR=%.9f  scaled_grad=%.6e  msgs=%d  opt_code=%s  singular=%s\n",
              exp(fixef(m)[[TERM]]), sg(m), length(m@optinfo$conv$lme4$messages),
              format(m@optinfo$conv$opt), isSingular(m)))
  for (msg in m@optinfo$conv$lme4$messages) cat("   base msg: ", msg, "\n", sep = "")

  cat("\n--- allFit ---\n")
  set.seed(20260728)
  af <- suppressWarnings(suppressMessages(allFit(m, verbose = FALSE)))
  nm <- names(af)
  cat(sprintf("attempted (length of allFit object) = %d\n", length(af)))
  cat("optimizer names: ", paste(nm, collapse = ", "), "\n")

  ok <- sapply(af, is, "merMod")
  cat(sprintf("returned a fit object (is merMod)   = %d\n", sum(ok)))

  rows <- list()
  for (i in seq_along(af)) {
    o <- af[[i]]
    if (!ok[i]) {
      cat(sprintf("\n[%s] DID NOT RETURN A FIT. class=%s\n", nm[i], paste(class(o), collapse = "/")))
      cat("   error text: ", paste(as.character(o), collapse = " | "), "\n")
      rows[[length(rows)+1]] <- data.table(optimizer = nm[i], returned_fit = FALSE,
        n_msgs = NA_integer_, msgs = NA_character_, OR = NA_real_,
        logLik = NA_real_, scaled_grad = NA_real_, opt_code = NA_character_)
      next
    }
    ms <- o@optinfo$conv$lme4$messages
    if (is.null(ms)) ms <- character(0)
    orv <- exp(fixef(o)[[TERM]])
    ll  <- as.numeric(logLik(o))
    cat(sprintf("\n[%s] returned=TRUE  n_msgs=%d  OR=%.9f  logLik=%.9f  scaled_grad=%.6e  opt_code=%s  singular=%s\n",
                nm[i], length(ms), orv, ll, sg(o), format(o@optinfo$conv$opt), isSingular(o)))
    if (length(ms)) for (k in seq_along(ms)) cat(sprintf("   msg[%d]: %s\n", k, ms[k]))
    rows[[length(rows)+1]] <- data.table(optimizer = nm[i], returned_fit = TRUE,
      n_msgs = length(ms), msgs = paste(ms, collapse = " || "), OR = orv,
      logLik = ll, scaled_grad = sg(o), opt_code = format(o@optinfo$conv$opt))
  }
  tab <- rbindlist(rows)

  cat("\n--- summary(allFit) internals ---\n")
  s <- summary(af)
  cat("which.OK: ", paste(names(s$which.OK), s$which.OK, sep = "=", collapse = "  "), "\n")
  cat("sum(which.OK) = ", sum(s$which.OK), "\n")
  cat("msgs element (summary(af)$msgs):\n"); print(s$msgs)

  ret <- tab[returned_fit == TRUE]
  clean <- ret[n_msgs == 0]
  cat("\n--- HEADLINE NUMBERS ---\n")
  cat(sprintf("attempted                      : %d\n", length(af)))
  cat(sprintf("returned a fit                 : %d\n", nrow(ret)))
  cat(sprintf("returned AND zero conv messages: %d\n", nrow(clean)))
  cat(sprintf("OR range over ALL returned fits    : %.9f  (min=%.9f max=%.9f)\n",
              max(ret$OR)-min(ret$OR), min(ret$OR), max(ret$OR)))
  cat(sprintf("logLik range over ALL returned fits: %.9f\n", max(ret$logLik)-min(ret$logLik)))
  if (nrow(clean) > 1) {
    cat(sprintf("OR range over MESSAGE-FREE fits    : %.9f  (min=%.9f max=%.9f)\n",
                max(clean$OR)-min(clean$OR), min(clean$OR), max(clean$OR)))
    cat(sprintf("logLik range over MESSAGE-FREE fits: %.9f\n", max(clean$logLik)-min(clean$logLik)))
  } else cat("MESSAGE-FREE subset has <2 fits; no range.\n")

  cat("\n--- per-optimizer table ---\n")
  print(tab[, .(optimizer, returned_fit, n_msgs, OR = sprintf("%.9f", OR),
                logLik = sprintf("%.9f", logLik), scaled_grad = sprintf("%.3e", scaled_grad), opt_code)])
  cat("\n--- full message text per optimizer ---\n")
  for (i in seq_len(nrow(tab))) cat(sprintf("%-32s | %s\n", tab$optimizer[i],
       if (is.na(tab$msgs[i])) "<no fit>" else if (nchar(tab$msgs[i])==0) "<none>" else tab$msgs[i]))
  invisible(tab)
}

t38 <- report("merged_cells.csv", "N=38 (unscreened)")
t36 <- report("merged_cells_screened.csv", "N=36 (screened)")

## ---- Task 3: subsample stability ------------------------------------------
cat("\n\n################ TASK 3: subsample 34/36, 40 draws, seed 20260728 ################\n")
dt36 <- prep("merged_cells_screened.csv")
ids <- levels(dt36$participant_id)
cat(sprintf("N participants in screened set = %d\n", length(ids)))
set.seed(20260728)
res <- rbindlist(lapply(1:40, function(k) {
  keep <- sample(ids, 34)
  d <- droplevels(dt36[participant_id %in% keep])
  out <- tryCatch({
    m <- fit_one(d)
    ms <- m@optinfo$conv$lme4$messages; if (is.null(ms)) ms <- character(0)
    g <- sg(m)
    data.table(draw = k, err = FALSE, n_msgs = length(ms), grad = g,
               lme4_pass = length(ms) == 0,
               newton_pass = is.finite(g) && g <= GRAD_TOL,
               opt0 = isTRUE(m@optinfo$conv$opt == 0),
               msgs = paste(ms, collapse = " || "))
  }, error = function(e) data.table(draw = k, err = TRUE, n_msgs = NA_integer_,
       grad = NA_real_, lme4_pass = FALSE, newton_pass = FALSE, opt0 = FALSE,
       msgs = conditionMessage(e)))
  out
}))
print(res[, .(draw, n_msgs, grad = sprintf("%.3e", grad), lme4_pass, newton_pass, opt0)])
cat(sprintf("\nlme4 criterion (zero conv messages)  passed: %d/40\n", sum(res$lme4_pass)))
cat(sprintf("Newton criterion (scaled_grad<=%.3f) passed: %d/40\n", GRAD_TOL, sum(res$newton_pass)))
cat(sprintf("joint gate (both + opt_code==0)      passed: %d/40\n",
            sum(res$lme4_pass & res$newton_pass & res$opt0)))
cat(sprintf("hard errors: %d/40\n", sum(res$err)))
cat("\nmessage texts seen:\n")
print(res[nchar(msgs) > 0, .N, by = msgs])
cat("\nsessionInfo:\n"); print(sessionInfo())
