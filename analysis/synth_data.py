#!/usr/bin/env python3
"""Synthetic behavioural data for the Follow-Ratio analysis pipeline (B2).

Stdlib only. Generates data/synth_cells.csv matching the locked analysis
spec so fit_glmm.R can be exercised end to end before real data arrives.

Design
------
Between-subjects 1x2: condition in {llm-reasoning, llm-control}.
N=40 (20 per arm), each participant has 5 months x 4 rounds = 20 decision cells.

Bidirectional DV (Plan X, 2026-07-07)
-------------------------------------
The advisor recommends either a FOLLOW-direction action (buy / hold) or an
UNFOLLOW-direction action (divest / step back). We fix `round == 4` in every
month as the UNFOLLOW-direction cell (5 per participant, standing in for the
bearish rounds); the other 15 cells are FOLLOW-direction. A cell counts as
"followed" (primary DV = 1) when the participant takes the *direction-compliant*
adoption:

  FOLLOW-direction cell:
    - follow      -> followed=1, decision_type=follow, invests a positive amount,
                     cas defined in (0,1]
    - decline     -> followed=0, decision_type=decline, "Never mind" (invested NULL, cas NULL)
    - zero_optout -> followed=0, decision_type=zero_optout, commits GBP 0 (invested 0, cas NULL)
  UNFOLLOW-direction cell:
    - unfollow    -> followed=1, decision_type=unfollow, records a POSITIVE
                     divestment `invested` but NO cas (a share-of-available is
                     undefined for a withdrawal)
    - decline     -> followed=0, decision_type=decline (recommendation not adopted;
                     `direction` stays 'unfollow' because it labels the *advice*,
                     not the choice). No zero_optout on unfollow-direction cells.

`followed` is generated from a random-intercept logistic model so the GLMM
`followed ~ condition + literacy + is_HCI + factor(month) + (1|participant)` can
recover the planted effect. The condition effect is planted on the CONDITIONAL
(subject-specific) logit scale at log(1.3), i.e. the estimand the GLMM fixed
effect targets. At sigma_u=0.45 the marginal (population-averaged) OR
attenuation is mild -- roughly 1.29, not far from 1.30 -- so the gap between the
raw-rate OR and 1.30 is driven mainly by N=40 sampling noise, not by
attenuation. The recovery check is that the GLMM OR's CI covers 1.3, not that
the raw rates reproduce it. The adoption probability p_follow is applied
identically to follow- and unfollow-direction cells, so the arm contrast is
carried by both directions symmetrically.

Participant-level covariates
----------------------------
`literacy` (Lusardi Big-Three score, 0-3) and `is_HCI` (0/1) are drawn ONCE per
participant, independently of `followed`. They carry NO planted effect -- their
GLMM coefficients are expected to sit near 0. They exist so the pipeline
exercises the pre-registered covariate-adjusted model on synthetic data; on the
real path prepare_cells.R derives them from the pre-task questionnaire.

cas is the conditional-on-follow continuous share = invested / available_snapshot,
defined ONLY for follow rows (decision_type=='follow'). For decline, zero_optout,
and unfollow it is a STRUCTURAL NULL, not a missing value. fit_glmm.R relies on
that contract to separate structural NA from true missingness.
"""

import csv
import math
import os
import random
import statistics

SEED = 20260708
N_PER_ARM = 20
CELLS_PER_PARTICIPANT = 20          # 5 months x 4 rounds
MONTHS = 5
ROUNDS = 4
UNFOLLOW_ROUND = 4                  # round 4 of every month is the unfollow-direction cell
UNFOLLOW_CELLS_PER_PARTICIPANT = MONTHS   # one per month -> 5

# followed ~ Bernoulli(logistic(B0 + B1*reasoning + u_participant))
B0 = 1.20                            # control baseline logit  (logistic(1.20)=0.769)
B1 = math.log(1.30)                  # planted conditional OR = 1.3 (reasoning vs control)
SIGMA_U = 0.45                       # participant random-intercept SD (logit scale)

# among FOLLOW-direction non-follows, split into decline vs zero_optout
P_DECLINE_GIVEN_NOFOLLOW = 0.65      # -> ~15% decline, ~8% opt-out of total at ~78% follow

# participant-level covariate distributions (no effect on followed)
LITERACY_LEVELS = [0, 1, 2, 3]
LITERACY_WEIGHTS = [0.10, 0.20, 0.35, 0.35]   # skewed toward the numerate, a realistic pool
P_IS_HCI = 0.35

# cas share ~ Beta; reasoning arm invests a slightly larger share (secondary, descriptive only)
CAS_BETA = {
    "llm-control":   (2.0, 2.6),     # mean ~0.435
    "llm-reasoning": (2.3, 2.6),     # mean ~0.469
}

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "synth_cells.csv")

FIELDNAMES = [
    "participant_id", "condition", "month", "round", "cell_index",
    "decision_type", "direction", "followed", "available_snapshot",
    "invested", "cas", "literacy", "is_HCI",
]


def logistic(x):
    return 1.0 / (1.0 + math.exp(-x))


def build_rows():
    rng = random.Random(SEED)
    rows = []
    # 20 control participants then 20 reasoning; interleave ids so id order != arm order
    participants = (
        [("llm-control", i) for i in range(N_PER_ARM)]
        + [("llm-reasoning", i) for i in range(N_PER_ARM)]
    )
    rng.shuffle(participants)

    for pnum, (condition, _) in enumerate(participants, start=1):
        pid = f"P{pnum:02d}"
        u = rng.gauss(0.0, SIGMA_U)                       # participant random intercept
        eta = B0 + (B1 if condition == "llm-reasoning" else 0.0) + u
        p_follow = logistic(eta)
        a_beta, b_beta = CAS_BETA[condition]

        # participant-level covariates, drawn once and held constant across the 20 cells
        literacy = rng.choices(LITERACY_LEVELS, weights=LITERACY_WEIGHTS, k=1)[0]
        is_HCI = 1 if rng.random() < P_IS_HCI else 0

        cell = 0
        for month in range(1, MONTHS + 1):
            for rnd in range(1, ROUNDS + 1):
                cell += 1
                available = round(rng.uniform(50.0, 150.0) + month * 10.0, 2)

                if rnd == UNFOLLOW_ROUND:
                    # UNFOLLOW-direction cell: adopt the divestment, or decline it.
                    direction = "unfollow"
                    if rng.random() < p_follow:
                        decision = "unfollow"
                        followed = 1
                        # positive divestment amount; NOT a share of `available`
                        # (a withdrawal can exceed the snapshot), so cas is undefined
                        divested = round(rng.uniform(10.0, 80.0), 2)
                        invested_out = f"{divested:.2f}"
                        cas_out = ""                      # structural NULL: no share for a withdrawal
                    else:
                        decision = "decline"              # unfollow advice not adopted
                        followed = 0
                        invested_out = ""                 # structural NULL
                        cas_out = ""                      # structural NULL
                else:
                    # FOLLOW-direction cell: original three-way outcome.
                    direction = "follow"
                    if rng.random() < p_follow:
                        decision = "follow"
                        followed = 1
                        share = rng.betavariate(a_beta, b_beta)
                        invested = round(share * available, 2)
                        if invested <= 0.0:               # rounding guard: a follow always commits > 0
                            invested = 0.01
                        cas = invested / available
                        invested_out = f"{invested:.2f}"
                        # 10dp so cas matches a recomputed invested/available within the
                        # fit_glmm.R 1e-9 consistency guard (6dp drifts ~5e-7 and trips it)
                        cas_out = f"{cas:.10f}"
                    elif rng.random() < P_DECLINE_GIVEN_NOFOLLOW:
                        decision = "decline"              # Never mind: no amount entered
                        followed = 0
                        invested_out = ""                 # structural NULL
                        cas_out = ""                      # structural NULL
                    else:
                        decision = "zero_optout"          # engaged, committed GBP 0
                        followed = 0
                        invested_out = "0.00"             # explicit 0, distinct from decline's NULL
                        cas_out = ""                      # structural NULL (conditional-on-follow)

                rows.append({
                    "participant_id": pid,
                    "condition": condition,
                    "month": month,
                    "round": rnd,
                    "cell_index": cell,
                    "decision_type": decision,
                    "direction": direction,
                    "followed": followed,
                    "available_snapshot": f"{available:.2f}",
                    "invested": invested_out,
                    "cas": cas_out,
                    "literacy": literacy,
                    "is_HCI": is_HCI,
                })
    return rows


def self_check(rows):
    """Assert the planted structure survived generation. Fails loudly if not."""
    n_rows = len(rows)
    assert n_rows == 2 * N_PER_ARM * CELLS_PER_PARTICIPANT, f"row count {n_rows} != 800"

    by_pid = {}
    for r in rows:
        by_pid.setdefault(r["participant_id"], []).append(r)
    assert len(by_pid) == 2 * N_PER_ARM, f"participants {len(by_pid)} != 40"
    for pid, rs in by_pid.items():
        assert len(rs) == CELLS_PER_PARTICIPANT, f"{pid} has {len(rs)} cells != 20"
        assert len({r["condition"] for r in rs}) == 1, f"{pid} spans >1 condition"
        # exactly one unfollow-direction cell per month -> 5 per participant
        n_unf_dir = sum(r["direction"] == "unfollow" for r in rs)
        assert n_unf_dir == UNFOLLOW_CELLS_PER_PARTICIPANT, \
            f"{pid} has {n_unf_dir} unfollow-direction cells != {UNFOLLOW_CELLS_PER_PARTICIPANT}"
        # participant-level covariates constant within participant, and non-empty
        assert len({r["literacy"] for r in rs}) == 1, f"{pid} spans >1 literacy value"
        assert len({r["is_HCI"] for r in rs}) == 1, f"{pid} spans >1 is_HCI value"
        assert rs[0]["literacy"] in LITERACY_LEVELS, f"{pid} literacy out of range"
        assert rs[0]["is_HCI"] in (0, 1), f"{pid} is_HCI out of range"

    arm_counts = {}
    for pid, rs in by_pid.items():
        arm_counts[rs[0]["condition"]] = arm_counts.get(rs[0]["condition"], 0) + 1
    assert arm_counts.get("llm-control") == N_PER_ARM, f"control arm {arm_counts}"
    assert arm_counts.get("llm-reasoning") == N_PER_ARM, f"reasoning arm {arm_counts}"

    n_follow = sum(r["decision_type"] == "follow" for r in rows)
    n_unfollow = sum(r["decision_type"] == "unfollow" for r in rows)
    n_decline = sum(r["decision_type"] == "decline" for r in rows)
    n_opt = sum(r["decision_type"] == "zero_optout" for r in rows)
    assert n_follow + n_unfollow + n_decline + n_opt == n_rows, "decision types do not partition rows"

    n_followed = sum(r["followed"] for r in rows)
    assert n_followed == n_follow + n_unfollow, "followed count != follow+unfollow"

    # 5 unfollow-direction cells x 40 participants = 200 unfollow-direction cells total
    n_unf_dir_total = sum(r["direction"] == "unfollow" for r in rows)
    assert n_unf_dir_total == UNFOLLOW_CELLS_PER_PARTICIPANT * 2 * N_PER_ARM, \
        f"unfollow-direction cells {n_unf_dir_total} != 200"

    # fractions -- bands widened for the bidirectional split (follow-dir=75%, unfollow-dir=25% of cells)
    followed_frac = n_followed / n_rows
    follow_frac = n_follow / n_rows
    unfollow_frac = n_unfollow / n_rows
    decline_frac = n_decline / n_rows
    opt_frac = n_opt / n_rows
    assert 0.70 <= followed_frac <= 0.85, f"followed frac {followed_frac:.3f} out of band"
    assert 0.48 <= follow_frac <= 0.68, f"follow frac {follow_frac:.3f} out of band"
    assert 0.12 <= unfollow_frac <= 0.28, f"unfollow frac {unfollow_frac:.3f} out of band"
    assert 0.09 <= decline_frac <= 0.24, f"decline frac {decline_frac:.3f} out of band"
    assert 0.02 <= opt_frac <= 0.12, f"opt-out frac {opt_frac:.3f} out of band"

    # bidirectional structural contract: followed==1 iff decision_type in {follow, unfollow};
    # cas defined iff decision_type=='follow'; direction labels the advice, not the choice.
    for r in rows:
        if r["followed"] == 1:
            assert r["decision_type"] in ("follow", "unfollow"), "followed=1 but decision_type not follow/unfollow"
            if r["decision_type"] == "follow":
                assert r["direction"] == "follow", "follow decision on a non-follow-direction cell"
                assert r["cas"] != "", "TRUE MISSING: follow row with empty cas"
                assert r["invested"] not in ("", "0.00"), "follow row invested is null/zero"
                share = float(r["cas"])
                assert 0.0 < share <= 1.0, f"cas {share} out of (0,1]"
            else:  # unfollow
                assert r["direction"] == "unfollow", "unfollow decision on a non-unfollow-direction cell"
                assert r["cas"] == "", "unfollow row must have structural NULL cas"
                assert r["invested"] != "" and float(r["invested"]) > 0.0, "unfollow row must record positive invested"
        else:
            assert r["decision_type"] in ("decline", "zero_optout")
            assert r["cas"] == "", "non-follow row has non-null cas (structural violation)"
            if r["decision_type"] == "decline":
                assert r["invested"] == "", "decline row must have NULL invested"
            else:
                assert r["invested"] == "0.00", "zero_optout row must have invested 0.00"
        assert r["available_snapshot"] != "", "available_snapshot must never be null"

    # raw (marginal) arm follow rates + descriptive marginal OR, for the run log.
    # "follow rate" here = adoption rate (followed==1 across BOTH directions).
    def arm_rate(arm):
        rs = [r for r in rows if r["condition"] == arm]
        return sum(r["followed"] for r in rs) / len(rs)

    rate_c = arm_rate("llm-control")
    rate_r = arm_rate("llm-reasoning")
    odds_c = rate_c / (1 - rate_c)
    odds_r = rate_r / (1 - rate_r)
    marginal_or = odds_r / odds_c

    cas_vals = {
        arm: [float(r["cas"]) for r in rows
              if r["condition"] == arm and r["decision_type"] == "follow" and r["cas"] != ""]
        for arm in ("llm-control", "llm-reasoning")
    }

    print("=== synth_data self-check PASSED ===")
    print(f"rows={n_rows}  participants={len(by_pid)}  cells/participant={CELLS_PER_PARTICIPANT}")
    print(f"follow={follow_frac:.3f}  unfollow={unfollow_frac:.3f}  decline={decline_frac:.3f}  zero_optout={opt_frac:.3f}")
    print(f"adoption (followed==1): follow={n_follow} unfollow={n_unfollow} total={n_followed}  unfollow-direction cells={n_unf_dir_total}")
    print(f"raw adoption rate  control={rate_c:.3f}  reasoning={rate_r:.3f}")
    print(f"marginal(raw) OR={marginal_or:.3f}  (theory ~1.29 at sigma_u=0.45; residual gap = N=40 noise; planted CONDITIONAL OR=1.300)")
    lit_by_pid = sorted(rs[0]["literacy"] for rs in by_pid.values())
    hci_by_pid = [rs[0]["is_HCI"] for rs in by_pid.values()]
    print(f"literacy dist (per participant): {[lit_by_pid.count(k) for k in LITERACY_LEVELS]} for levels {LITERACY_LEVELS}  | is_HCI mean={statistics.mean(hci_by_pid):.2f}")
    for arm in ("llm-control", "llm-reasoning"):
        v = cas_vals[arm]
        print(f"cas[{arm}] n={len(v)} mean={statistics.mean(v):.3f} median={statistics.median(v):.3f}")


def main():
    rows = build_rows()
    self_check(rows)
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(rows)
    # re-read to confirm the file on disk matches what we validated in memory
    with open(OUT_PATH, newline="") as f:
        disk_rows = list(csv.DictReader(f))
    assert len(disk_rows) == len(rows), "disk row count != in-memory row count"
    print(f"wrote {len(disk_rows)} rows -> {OUT_PATH}")


if __name__ == "__main__":
    main()
