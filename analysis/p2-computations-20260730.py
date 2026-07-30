# P2 批准项计算（2026-07-30）：month×arm、方向标准化 bound、臂内相关、风险集分布、MDE/successor 算术
# 数据源 = exploratory-注册族-20260729/ 的 canonical CSV（冻结）
import csv, math, statistics
base = "exploratory-注册族-20260729/"
cells = list(csv.DictReader(open(base + "cells_with_complexity.csv")))
ps    = list(csv.DictReader(open(base + "participant_summary.csv")))
out = []
def w(s): out.append(s); print(s)

# 1) month×arm follow rate
w("## month×arm follow rate")
for arm in ["llm-reasoning", "llm-control"]:
    row = []
    for m in range(1, 6):
        sub = [r for r in cells if r["condition"] == arm and int(r["month"]) == m]
        f = sum(int(r["followed"]) for r in sub)
        row.append(f"M{m} {f}/{len(sub)}={f/len(sub):.3f}")
    w(arm + ": " + " | ".join(row))

# 2) direction-standardised gap (pooled-mix standardisation)
w("## direction-standardised arm gap")
def rate(arm, d):
    sub = [r for r in cells if r["condition"] == arm and r["direction"] == d]
    return sum(int(r["followed"]) for r in sub) / len(sub), len(sub)
rf, nf_r = rate("llm-reasoning", "follow"); rs, ns_r = rate("llm-reasoning", "unfollow")
cf, nf_c = rate("llm-control", "follow");   cs, ns_c = rate("llm-control", "unfollow")
wf = (nf_r + nf_c) / len(cells); ws = (ns_r + ns_c) / len(cells)
std_r = wf * rf + ws * rs; std_c = wf * cf + ws * cs
raw_r = sum(int(r["followed"]) for r in cells if r["condition"]=="llm-reasoning") / (nf_r + ns_r)
raw_c = sum(int(r["followed"]) for r in cells if r["condition"]=="llm-control") / (nf_c + ns_c)
w(f"raw gap = {raw_r:.4f}-{raw_c:.4f} = {raw_r-raw_c:.4f}")
w(f"pooled weights follow={wf:.4f} stepback={ws:.4f}")
w(f"standardised: reasoning {std_r:.4f}, control {std_c:.4f}, gap {std_r-std_c:.4f}")

# 3) within-arm correlations (scale vs follow_ratio)
w("## within-arm correlations (r, n)")
def pearson(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    sx = math.sqrt(sum((a-mx)**2 for a in x)); sy = math.sqrt(sum((b-my)**2 for b in y))
    return sum((a-mx)*(b-my) for a, b in zip(x, y)) / (sx * sy)
for arm in ["llm-reasoning", "llm-control"]:
    line = []
    for k in ["U", "R", "TR", "SP", "MC1"]:
        pairs = [(float(r[k]), float(r["follow_ratio"])) for r in ps if r["condition"]==arm and r[k] not in ("", "NA")]
        line.append(f"{k} r={pearson([a for a,_ in pairs],[b for _,b in pairs]):+.3f}(n={len(pairs)})")
    w(arm + ": " + " | ".join(line))

# 4) per-arm cells-per-participant distribution
w("## per-arm cells per participant")
for arm in ["llm-reasoning", "llm-control"]:
    ns = sorted(int(r["n_cells"]) for r in ps if r["condition"]==arm)
    w(f"{arm}: total={sum(ns)} median={statistics.median(ns)} range=[{ns[0]},{ns[-1]}]")
w(f"decision-bearing share of scheduled rounds: 387/720 = {387/720:.3f}")

# 5) MDE + successor N arithmetic (from frozen SE=0.432, two-sided alpha=.10, power .80)
w("## MDE / successor arithmetic")
se = 0.432; z_a = 1.6449; z_b = 0.8416
mde = (z_a + z_b) * se
w(f"post-hoc MDE = ({z_a}+{z_b})*{se} = {mde:.4f} log-odds -> OR {math.exp(mde):.3f}")
hw = z_a * se
w(f"realised 90% half-width = {hw:.4f} log-odds (OR factor {math.exp(hw):.3f})")
w(f"halve half-width -> 4x participants ~= {36*4}")
se_t = 0.30 / z_a
w(f"target half-width 0.30 -> SE {se_t:.4f} -> ({se}/{se_t:.4f})^2 = {(se/se_t)**2:.2f}x -> N ~= {36*(se/se_t)**2:.0f}")
open("p2-computations-输出-20260730.txt", "w").write("\n".join(out))
