# -*- coding: utf-8 -*-
"""Fig 4d: HGT events vs phylum abundance scatter plot"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
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

df = data["df"]
EPS = 1e-6

# ---- 数据准备 ----
d_sub = df[df["donor_events"] > 0].copy()
r_sub = df[df["recipient_events"] > 0].copy()

d_x = np.log10(d_sub["donor_abund"] + EPS)
d_y = np.log10(d_sub["donor_events"] + 1)
r_x = np.log10(r_sub["recipient_abund"] + EPS)
r_y = np.log10(r_sub["recipient_events"] + 1)

DONOR_COLOR = "#6A3D9A"
RECIPIENT_COLOR = "#FF7F00"

# ---- 绘图 ----
# 尺寸/字体复刻原总图 (plot_fig4_v50.py): d 面板 0.38x0.28 于 183x220mm, 全框
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

ax.scatter(d_x, d_y, s=20, color=DONOR_COLOR, alpha=0.55, edgecolor="none", label="Donor", zorder=3)
ax.scatter(r_x, r_y, s=20, color=RECIPIENT_COLOR, alpha=0.55, edgecolor="none", label="Recipient", zorder=3)

z_d = np.polyfit(d_x, d_y, 1)
z_r = np.polyfit(r_x, r_y, 1)
xs_line = np.linspace(min(d_x.min(), r_x.min()), max(d_x.max(), r_x.max()), 100)
ax.plot(xs_line, np.polyval(z_d, xs_line), color=DONOR_COLOR, lw=1.5, ls="-", label="Donor fit")
ax.plot(xs_line, np.polyval(z_r, xs_line), color=RECIPIENT_COLOR, lw=1.5, ls="--", label="Recipient fit")

rho_d, p_d = stats.spearmanr(d_sub["donor_abund"], d_sub["donor_events"])
rho_r, p_r = stats.spearmanr(r_sub["recipient_abund"], r_sub["recipient_events"])

ax.text(0.97, 0.97,
        f"Donor: ρ = {rho_d:.2f}, P = {p_d:.3g} (n = {len(d_sub)})\n"
        f"Recipient: ρ = {rho_r:.2f}, P = {p_r:.3g} (n = {len(r_sub)})",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, linespacing=1.5)

ax.set_xlabel("log$_{10}$ phylum abundance (%)")
ax.set_ylabel("log$_{10}$ (HGT events + 1)")
all_x = np.concatenate([d_x, r_x])
ax.set_xlim(min(all_x.min(), 0) - 0.2, np.ceil(all_x.max()))
ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", labelsize=9)
ax.set_title("d", fontsize=11, fontweight="bold", loc="left")
ax.legend(fontsize=8, frameon=True, loc="lower right", edgecolor="#CCCCCC")

# 全框 (与原总图一致)

fig.savefig(OUT / "04d_hgt_scatter.png", dpi=300, facecolor="white")
fig.savefig(OUT / "04d_hgt_scatter.pdf", facecolor="white")
print("saved 04d_hgt_scatter.png/pdf")
plt.close()