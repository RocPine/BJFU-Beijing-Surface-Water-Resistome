# -*- coding: utf-8 -*-
"""Fig 4e: Donor/Recipient/Other ARG host NCG violin plot (N_total)"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fig4_data_loader import load_fig4_data

import yaml, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

data = load_fig4_data(cfg)

donor_n_e = data["donor_n_e"]
recipient_n_e = data["recipient_n_e"]
non_donor_n_e = data["non_donor_n_e"]

DONOR_COLOR = "#6A3D9A"
RECIPIENT_COLOR = "#FF7F00"
NONDONOR_COLOR = "#7F7F7F"
EDGE_COLOR = "#333333"

groups_e = [
    (donor_n_e, "Donor\n(n=32)", DONOR_COLOR),
    (recipient_n_e, "Recipient\n(n=25)", RECIPIENT_COLOR),
    (non_donor_n_e, "Other ARG host\n(n=622)", NONDONOR_COLOR),
]

positions = [0.2, 0.5, 0.8]

# 尺寸/字体复刻原总图 (plot_fig4_v50.py): e 面板 0.38x0.28 于 183x220mm, 全框无 top/right 刻度
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})
fig, ax = plt.subplots(figsize=(69.5/25.4, 61.6/25.4), dpi=300)

for i, (data_arr, label, color) in enumerate(groups_e):
    vp = ax.violinplot(
        data_arr, positions=[positions[i]], widths=0.35,
        showmeans=False, showmedians=False, showextrema=False,
        bw_method=0.3
    )
    for body in vp["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor(EDGE_COLOR)
        body.set_linewidth(0.5)
        body.set_alpha(0.6)

    bp = ax.boxplot(
        [data_arr], positions=[positions[i]], widths=0.08,
        patch_artist=True, showfliers=False,
        medianprops=dict(color="white", linewidth=1.5),
        boxprops=dict(facecolor="white", edgecolor=EDGE_COLOR, linewidth=0.8),
        whiskerprops=dict(color=EDGE_COLOR, linewidth=0.6),
        capprops=dict(color=EDGE_COLOR, linewidth=0.6)
    )

    jitter = np.random.uniform(-0.03, 0.03, len(data_arr))
    ax.scatter(
        np.full(len(data_arr), positions[i]) + jitter,
        data_arr, s=3, color=color, alpha=0.3, edgecolor="none", zorder=3
    )

    ax.text(positions[i], np.median(data_arr)*1.02, f"{np.median(data_arr):.0f}",
            ha="center", va="bottom", fontsize=6, color=color, fontweight="bold")

# Stats
stat_d_e, p_d_e = stats.mannwhitneyu(donor_n_e, non_donor_n_e, alternative="greater")
stat_r_e, p_r_e = stats.mannwhitneyu(recipient_n_e, non_donor_n_e, alternative="greater")

y_max_e = max(np.max(donor_n_e), np.max(recipient_n_e), np.max(non_donor_n_e))

ax.plot([0.2, 0.2, 0.8, 0.8], [y_max_e*1.08, y_max_e*1.11, y_max_e*1.11, y_max_e*1.08],
        color="black", linewidth=0.6)
ax.text(0.35, y_max_e*1.12, f"MW P = {p_d_e:.2e}", ha="center", va="bottom", fontsize=6)

ax.plot([0.5, 0.5, 0.8, 0.8], [y_max_e*1.18, y_max_e*1.21, y_max_e*1.21, y_max_e*1.18],
        color="black", linewidth=0.6)
ax.text(0.65, y_max_e*1.22, f"MW P = {p_r_e:.2e}", ha="center", va="bottom", fontsize=6)

ax.set_xticks([0.2, 0.5, 0.8])
ax.set_xticklabels([
    "Donor\n(n=32)",
    "Recipient\n(n=25)",
    "Other ARG host\n(n=622)"
], fontsize=7, linespacing=1.3)

ax.set_ylabel("NCGs per SGB", fontsize=8)
ax.set_ylim(0, max(np.max(donor_n_e), np.max(recipient_n_e), np.max(non_donor_n_e)) * 1.4)
ax.set_xlim(0, 1.0)
ax.set_title("e", fontsize=11, fontweight="bold", loc="left")

# 全框 (与原总图一致: 4 边可见, 无 top/right 刻度线)
for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.6)
ax.tick_params(axis="x", top=False, labeltop=False)
ax.tick_params(axis="y", right=False, labelright=False)

fig.savefig(OUT / "04e_hgt_violin.png", dpi=300, facecolor="white")
fig.savefig(OUT / "04e_hgt_violin.pdf", facecolor="white")
print("saved 04e_hgt_violin.png/pdf")
plt.close()