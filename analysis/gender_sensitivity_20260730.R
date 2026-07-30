# Post-hoc gender-adjusted sensitivity (2026-07-30, DA C3 部分采纳)
# (a) registered GLMM + gender; (b) CAS mean difference adjusted for gender
# 输入 = 冻结管道 merged_cells_screened.csv + questionnaire_items.csv 的 F2
suppressMessages({library(data.table); library(lme4)})
dt <- fread("data/real-20260727/merged_cells_screened.csv", na.strings = "")
qi <- fread("data/real-20260727/questionnaire_items.csv", na.strings = "")
g  <- qi[phase == "pretask" & grepl("^F2", question) & answer %in% c("Woman", "Man"), .(participant_id, gender = answer)]
g  <- unique(g)
stopifnot(all(g$gender %in% c("Woman", "Man")))
dt <- merge(dt, g, by = "participant_id", all.x = TRUE)
stopifnot(!anyNA(dt$gender), uniqueN(dt$participant_id) == 36, nrow(dt) == 387)
dt[, is_man := as.integer(gender == "Man")]
dt[, condition := relevel(factor(condition), ref = "llm-control")]

m <- glmer(followed ~ condition + literacy + is_HCI + is_man +
             factor(month) + factor(round) + (1 | participant_id),
           data = dt, family = binomial,
           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
b  <- fixef(m)["conditionllm-reasoning"]
se <- sqrt(diag(vcov(m)))["conditionllm-reasoning"]
cat(sprintf("GLMM+gender: OR=%.3f  90%%CI=[%.3f, %.3f]  beta=%.4f se=%.4f\n",
            exp(b), exp(b - 1.6449 * se), exp(b + 1.6449 * se), b, se))
cat(sprintf("gender coef (man): OR=%.3f\n", exp(fixef(m)["is_man"])))
msgs <- unlist(summary(m)$optinfo$conv$lme4$messages); cat("convergence messages:", if (is.null(msgs)) "none" else msgs, "\n")

# (b) CAS: follow-direction adopted cells, OLS + participant cluster bootstrap
cc <- dt[direction == "follow" & followed == 1 & !is.na(cas)]
cat(sprintf("CAS cells: %d reasoning / %d control\n",
            nrow(cc[condition == "llm-reasoning"]), nrow(cc[condition == "llm-control"])))
adj_diff <- function(d) {
  f <- lm(cas ~ condition + is_man, data = d)
  coef(f)["conditionllm-reasoning"]
}
set.seed(20260727)
ids <- split(seq_len(nrow(cc)), cc$participant_id)
boots <- replicate(2000, {
  samp <- sample(names(ids), length(ids), replace = TRUE)
  adj_diff(cc[unlist(ids[samp])])
})
pt <- adj_diff(cc)
cat(sprintf("CAS gender-adjusted diff: %+.3f  90%%CI=[%.3f, %.3f]\n",
            pt, quantile(boots, .05), quantile(boots, .95)))
