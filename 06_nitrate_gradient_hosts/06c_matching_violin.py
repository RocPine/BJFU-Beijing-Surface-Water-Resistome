# -*- coding: utf-8 -*-
"""Fig 6c: 丰度匹配后 ARG hosts vs matched non-hosts 的 NAR/NAP/NAS 基因数小提琴图"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "08_core_statistics"))
from narnas_stats import compute_matching

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_HOST = "#E64B35"; C_MATCH = "#00A087"

r = compute_matching(cfg)
h_cnt = r["host_counts"]
m_cnt = r["matched_counts"]
mean_diff = r["mean_diff"]
w_p = r["p_wilcoxon"]

fig, ax = plt.subplots(figsize=(230/72, 200/72), dpi=300)

vp = ax.violinplot([h_cnt, m_cnt], positions=[0, 1], widths=0.55,
                   showmeans=False, showmedians=False, showextrema=False, bw_method=0.35)
for body, c in zip(vp["bodies"], [C_HOST, C_MATCH]):
    body.set_facecolor(c); body.set_alpha(0.35); body.set_edgecolor(c); body.set_linewidth(0.6)

bp = ax.boxplot([h_cnt, m_cnt], positions=[0, 1], widths=0.16, patch_artist=True, showfliers=False,
                medianprops=dict(color="white", linewidth=1.2),
                boxprops=dict(facecolor=C_HOST, edgecolor="black", linewidth=0.6, alpha=0.9),
                whiskerprops=dict(color="black", linewidth=0.5),
                capprops=dict(color="black", linewidth=0.5))
bp["boxes"][1].set_facecolor(C_MATCH)

rng = np.random.default_rng(42)
x_host = rng.uniform(-0.12, 0.12, len(h_cnt))
x_match = 1 + rng.uniform(-0.12, 0.12, len(m_cnt))
for x1, x2, h, m in zip(x_host, x_match, h_cnt, m_cnt):
    ax.plot([x1, x2], [h, m], color="#808080", alpha=0.12, linewidth=0.35, zorder=1)
ax.scatter(x_host, h_cnt, s=5, color="#808080", alpha=0.45,
           edgecolor="none", zorder=3)
ax.scatter(x_match, m_cnt, s=5, color="#808080", alpha=0.45,
           edgecolor="none", zorder=3)

y_max_c = max(h_cnt.max(), m_cnt.max())
ax.text(0.5, y_max_c + 2.5, f"Mean diff = +{mean_diff:.2f}\nP = {w_p:.1e}", ha="center", fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels([f"ARG hosts\n(n = {len(h_cnt)})", f"Matched\nnon-hosts\n(n = {len(m_cnt)})"],
                   fontsize=7, linespacing=1.3)
ax.set_ylabel("NAR/NAP/NAS genes per SGB", fontsize=8)
ax.set_ylim(0, y_max_c + 6)
ax.set_yticks([5, 10, 15, 20, 25, 30, 35])
ax.text(-0.12, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False)
ax.tick_params(axis="y", right=False, labelright=False)

fig.savefig(OUT / "06c_matching_violin.png", dpi=300, facecolor="white")
fig.savefig(OUT / "06c_matching_violin.pdf", facecolor="white")
print("saved 06c_matching_violin.png/pdf")
plt.close()
