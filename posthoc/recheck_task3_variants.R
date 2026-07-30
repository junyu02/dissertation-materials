suppressMessages({ library(data.table); library(lme4) })
DIR <- "/Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道/data/real-20260727"
dt <- fread(file.path(DIR,"merged_cells_screened.csv"), na.strings="")
dt[, round := as.integer(round)]
dt[, condition := relevel(factor(condition), ref="llm-control")]
dt[, participant_id := factor(participant_id)]
sg <- function(m){dd<-m@optinfo$derivs; tryCatch(max(abs(solve(dd$Hessian,dd$gradient))),error=function(e)NA_real_)}
run <- function(ids, tag){
  set.seed(20260728)
  r <- rbindlist(lapply(1:40, function(k){
    keep <- sample(ids, 34)
    d <- droplevels(dt[participant_id %in% keep])
    m <- glmer(followed ~ condition + literacy + is_HCI + factor(month) + factor(round) +
               (1|participant_id), data=d, family=binomial, nAGQ=1)
    ms <- m@optinfo$conv$lme4$messages; if(is.null(ms)) ms <- character(0)
    g <- sg(m)
    data.table(k=k, lme4=length(ms)==0, newton=is.finite(g)&&g<=0.002)
  }))
  cat(sprintf("%-34s lme4 %2d/40   newton %2d/40\n", tag, sum(r$lme4), sum(r$newton)))
}
suppressWarnings({
run(levels(dt$participant_id),                 "ids = levels() [sorted]")
run(unique(as.character(dt$participant_id)),   "ids = unique() [first-appearance]")
run(sort(unique(as.character(dt$participant_id))), "ids = sort(unique())")
})
