#!/usr/bin/env Rscript
# verify_explore_family.R -- independent recompute of the headline numbers in
# results-exploratory.md, by a deliberately different route (base R, raw CSVs,
# no data.table aggregation, no participant-summary intermediate). Asserts
# agreement with explore_family.R's outputs. Fails loud on any mismatch.
#
#   Rscript verify_explore_family.R [DATA_DIR] [OUT_DIR]

args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[[1]] else
  "/Users/wolpfor/Desktop/UCL毕业项目/00-论文交付/分析管道/data/real-20260727"
OUT_DIR  <- if (length(args) >= 2) args[[2]] else "/tmp/exploratory-family"

cells <- read.csv(file.path(DATA_DIR, "merged_cells_screened.csv"),
                  stringsAsFactors = FALSE, encoding = "UTF-8")
withc <- read.csv(file.path(OUT_DIR, "cells_with_complexity.csv"),
                  stringsAsFactors = FALSE, encoding = "UTF-8")
quest <- read.csv(file.path(DATA_DIR, "questionnaire_items.csv"),
                  stringsAsFactors = FALSE, encoding = "UTF-8")
revt  <- read.csv(file.path(DATA_DIR, "raw", "ReasoningEvent.csv"),
                  stringsAsFactors = FALSE, encoding = "UTF-8")

ok <- function(name, got, want, tol = 1e-9) {
  if (is.na(got) || is.na(want) || abs(got - want) > tol)
    stop(sprintf("MISMATCH %s: got %.10g, expected %.10g", name, got, want))
  cat(sprintf("  ok  %-46s %.6g\n", name, got))
}

cat("[verify] structure\n")
ok("n cells",              nrow(cells), 387)
ok("n participants",       length(unique(cells$participant_id)), 36)
ok("cells_with_complexity rows", nrow(withc), 387)
ok("tier cells sum",       sum(table(withc$complexity_tier)), 387)

cat("[verify] T1 CAS (decision_type=='follow' & !is.na(cas))\n")
cs <- cells[cells$decision_type == "follow" & !is.na(cells$cas), ]
mr <- mean(cs$cas[cs$condition == "llm-reasoning"])
mc <- mean(cs$cas[cs$condition == "llm-control"])
ok("CAS mean reasoning", mr, 0.505, 5e-4)
ok("CAS mean control",   mc, 0.315, 5e-4)
ok("CAS mean difference", mr - mc, 0.190, 5e-4)

cat("[verify] T9 absolute action size (invested = GBP size of the capital move)\n")
mv <- cells[!is.na(cells$invested) & cells$invested > 0, ]
ok("action-size cells reasoning", sum(mv$condition == "llm-reasoning"), 161)
ok("action-size cells control",   sum(mv$condition == "llm-control"),   124)
ok("action size mean reasoning", mean(mv$invested[mv$condition == "llm-reasoning"]), 460.660, 5e-4)
ok("action size mean control",   mean(mv$invested[mv$condition == "llm-control"]),   333.077, 5e-4)
ok("action size mean difference",
   mean(mv$invested[mv$condition == "llm-reasoning"]) - mean(mv$invested[mv$condition == "llm-control"]),
   127.583, 1e-3)
fo <- cells[cells$decision_type == "follow", ]
ok("follow-only action-size cells", nrow(fo), 199)   # same denominator as CAS
ok("follow-only mean difference",
   mean(fo$invested[fo$condition == "llm-reasoning"]) - mean(fo$invested[fo$condition == "llm-control"]),
   99.604, 1e-3)
ok("cells with no capital move", sum(is.na(cells$invested)) + sum(cells$invested == 0, na.rm = TRUE), 98 + 4)

cat("[verify] T2 complexity follow rates\n")
fr <- function(tier, arm) {
  s <- withc[withc$complexity_tier == tier & withc$condition == arm, ]
  sum(s$followed) / nrow(s)
}
ok("follow rate 2fo reasoning", fr("two-factor-offsetting", "llm-reasoning"), 42 / 58, 1e-9)
ok("follow rate 2fo control",   fr("two-factor-offsetting", "llm-control"),   32 / 55, 1e-9)
off  <- withc$complexity_tier == "two-factor-offsetting"
d_off <- sum(withc$followed[off & withc$condition == "llm-reasoning"]) / sum(off & withc$condition == "llm-reasoning") -
         sum(withc$followed[off & withc$condition == "llm-control"])   / sum(off & withc$condition == "llm-control")
d_non <- sum(withc$followed[!off & withc$condition == "llm-reasoning"]) / sum(!off & withc$condition == "llm-reasoning") -
         sum(withc$followed[!off & withc$condition == "llm-control"])   / sum(!off & withc$condition == "llm-control")
ok("DiD interaction contrast", d_off - d_non, 0.091, 5e-4)

cat("[verify] T3 reasoning telemetry (analysis set only)\n")
pids <- unique(cells$participant_id)
re <- revt[revt$user %in% pids, ]
flashed <- unique(re$block_id[re$event_type == "flash_shown"])
opened  <- unique(re$block_id[re$event_type == "open"])
ok("blocks flashed", length(flashed), 273)
ok("blocks opened",  length(intersect(opened, flashed)), 51)
ok("card-level open rate", length(intersect(opened, flashed)) / length(flashed), 51 / 273, 1e-9)
rp <- unique(cells$participant_id[cells$condition == "llm-reasoning"])
never <- sum(vapply(rp, function(u) sum(re$user == u & re$event_type == "open") == 0, logical(1)))
ok("reasoning participants never opening", never, 9)

cat("[verify] T5 MC1 / MC2 (earliest post-task submission per participant)\n")
po <- quest[quest$phase == "posttask" & quest$participant_id %in% pids, ]
first <- tapply(po$created_at, po$participant_id, min)
po <- po[po$created_at == first[po$participant_id], ]
arm <- setNames(cells$condition[!duplicated(cells$participant_id)],
                cells$participant_id[!duplicated(cells$participant_id)])
h1 <- po[startsWith(po$question, "H1."), ]
h1$arm <- arm[h1$participant_id]
ok("MC1 mean reasoning", mean(as.numeric(h1$answer[h1$arm == "llm-reasoning"])), 5.118, 5e-4)
ok("MC1 mean control",   mean(as.numeric(h1$answer[h1$arm == "llm-control"])),   2.444, 5e-4)
h2 <- po[startsWith(po$question, "H2."), ]
ok("MC2 n", nrow(h2), 17)
ok("MC2 correct", sum(h2$answer == "a central bank interest-rate rise"), 6)

cat("[verify] T7 ITT\n")
ok("ITT reasoning follow rate",
   sum(cells$followed[cells$condition == "llm-reasoning"]) / sum(cells$condition == "llm-reasoning"),
   161 / 208, 1e-9)
ok("ITT control follow rate",
   sum(cells$followed[cells$condition == "llm-control"]) / sum(cells$condition == "llm-control"),
   124 / 179, 1e-9)

cat("[verify] T10 manifest (funnel, cells/participant, Q6(b) latency, de-dup)\n")
qr <- read.csv(file.path(DATA_DIR, "raw", "QuestionnaireResponse.csv"), stringsAsFactors = FALSE)
pt <- read.csv(file.path(DATA_DIR, "raw", "Participant.csv"), stringsAsFactors = FALSE)
cu <- read.csv(file.path(DATA_DIR, "merged_cells.csv"), stringsAsFactors = FALSE)
ok("funnel: QuestionnaireResponse users", length(unique(qr$user)), 53)
ok("funnel: users with pretask",
   length(unique(quest$participant_id[quest$phase == "pretask"])), 52)
ok("funnel: Participant allocated", length(unique(pt$user)), 38)
ok("funnel: allocated reasoning", sum(pt$condition == "llm-reasoning"), 19)
ok("funnel: allocated control",   sum(pt$condition == "llm-control"),   19)
ok("funnel: unscreened cells", nrow(cu), 411)
ok("funnel: unscreened participants", length(unique(cu$participant_id)), 38)
ok("funnel: screened reasoning", length(unique(cells$participant_id[cells$condition == "llm-reasoning"])), 18)
ok("funnel: screened control",   length(unique(cells$participant_id[cells$condition == "llm-control"])),   18)
npc <- as.numeric(table(cells$participant_id))
ok("cells/participant median", median(npc), 11.5)
ok("cells/participant q1", unname(quantile(npc, .25)), 9)
ok("cells/participant q3", unname(quantile(npc, .75)), 13.25)
ok("cells/participant min", min(npc), 1)
ok("cells/participant max", max(npc), 17)
# Q6(b): recompute latency from the manifest's own pairing output rather than
# re-deriving the FK resolution, but assert the two headline figures independently.
lat <- read.csv(file.path(OUT_DIR, "decision_latency.csv"), stringsAsFactors = FALSE)
ok("Q6(b) paired cells", nrow(lat), 411)
ok("Q6(b) unpaired cells", sum(is.na(lat$latency_s)), 0)
ok("Q6(b) min latency (s)", min(lat$latency_s), 4.7517, 5e-5)
ok("Q6(b) sub-second decisions", sum(lat$latency_s < 1), 0)
po2 <- quest[quest$phase == "posttask" & quest$participant_id %in% pids, ]
subs <- tapply(po2$created_at, po2$participant_id, function(x) length(unique(x)))
ok("post-task duplicate submitters", sum(subs > 1), 1)
ok("post-task duplicate rows dropped", nrow(po2) - nrow(po), 20)

cat("[verify] T8 disconnect\n")
ok("disconnect cells",
   nrow(read.csv(file.path(DATA_DIR, "disconnect_cells.csv"), stringsAsFactors = FALSE)), 0)

cat("\n[verify] ALL CHECKS PASSED\n")
