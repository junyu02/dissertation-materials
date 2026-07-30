suppressMessages({library(data.table); library(logistf)})
TERM <- "conditionllm-reasoning"
load_dt <- function(f){dt<-fread(f,na.strings="");dt[,followed:=as.integer(followed)]
 dt[,condition:=relevel(factor(condition),ref="llm-control")][,participant_id:=factor(participant_id)][]}
firth_of <- function(d){
  agg <- d[, .(n_follow=sum(followed), n_total=.N), by=.(participant_id, condition, literacy, is_HCI)]
  long <- rbindlist(lapply(seq_len(nrow(agg)), function(i)
    data.table(y=c(1L,0L), w=c(agg$n_follow[i], agg$n_total[i]-agg$n_follow[i]),
               condition=agg$condition[i], literacy=agg$literacy[i], is_HCI=agg$is_HCI[i])))
  long <- long[w>0]; long[, condition := relevel(factor(condition), ref="llm-control")]
  f <- logistf(y ~ condition + literacy + is_HCI, data=long, weights=w, alpha=0.10)
  coef(f)[[which(names(coef(f))==TERM)]]}
for (nm in c("38","36")) {
  d <- load_dt(sprintf("/tmp/rv2/merged%s.csv", nm))
  set.seed(20260727); ids <- unique(d$participant_id); B <- 2000; out <- numeric(B)
  for (b in 1:B) { s <- sample(ids, length(ids), replace=TRUE)
    dbs <- rbindlist(lapply(seq_along(s), function(j){x<-copy(d[participant_id==s[j]]); x[,participant_id:=paste0("B",j)]; x}))
    out[b] <- tryCatch(firth_of(dbs), error=function(e) NA_real_) }
  out <- out[is.finite(out)]
  cat(sprintf("[N%s] seed=20260727 B=2000 -> 90%% percentile CI = [%.4f, %.4f]  median=%.4f\n",
      nm, exp(quantile(out,.05)), exp(quantile(out,.95)), exp(median(out))))
}
