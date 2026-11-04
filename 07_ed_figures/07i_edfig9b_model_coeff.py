# -*- coding: utf-8 -*-
"""ED Fig 9b: NO3-N 系数哑铃图 (m0/m1/m2 全站点 + m0 在 high-DO/low-DO 子集, β ± 95% CI)

ED Fig 9 = Fig 6 梯度关联的稳健性与特异性:
a=bootstrap ρ KDE, b=本图 (模型系数 + DO 分层), c=氮形态特异性
"""
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
from model_coefficients import compute_models

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_MAIN = "#E64B35"
C_SUBSET = "#4DBBD5"
C_NS = "#9E9E9E"

r = compute_models(cfg)
models = r["models"]
betas = [m[1] for m in models]
ses = [m[2] for m in models]
pvals = [m[3] for m in models]
do_models = r["do_models"]

xpos = [0.5, 1.55, 2.6, 3.75, 4.8]
xlabels = ["m0\nall sites", "m1\n+WQ", "m2\n+MEM",
           "m0\nhigh-DO", "m0\nlow-DO"]

fig, ax = plt.subplots(figsize=(230/72, 200/72), dpi=300)
ax.set_box_aspect(1)

# 全站点三模型
for i, (b, s, p) in enumerate(zip(betas, ses, pvals)):
    lo, hi = b - 1.96*s, b + 1.96*s
    c = C_MAIN if p < 0.05 else C_NS
    ax.vlines(xpos[i], lo, hi, color=c, lw=2.4, alpha=0.55, zorder=2)
    ax.scatter([xpos[i], xpos[i]], [lo, hi], s=18, color=c,
               edgecolor="black", linewidth=0.3, zorder=4)
    ax.scatter([xpos[i]], [b], s=22, marker="D", color="white",
               edgecolor=c, linewidth=0.7, zorder=5)
    ax.text(xpos[i], hi + 0.003, f"\u03b2 = {b:.3f}\nP = {p:.4g}",
            ha="center", va="bottom", fontsize=5, linespacing=1.2)

# m0 型模型在 DO 子集
for j, (name, b, s, p, n) in enumerate(do_models):
    i = 3 + j
    lo, hi = b - 1.96*s, b + 1.96*s
    c = C_SUBSET if p < 0.05 else C_NS
    ax.vlines(xpos[i], lo, hi, color=c, lw=2.4, alpha=0.55, zorder=2)
    ax.scatter([xpos[i], xpos[i]], [lo, hi], s=18, color=c,
               edgecolor="black", linewidth=0.3, zorder=4)
    ax.scatter([xpos[i]], [b], s=22, marker="D", color="white",
               edgecolor=c, linewidth=0.7, zorder=5)
    sig = "" if p < 0.05 else " (n.s.)"
    ax.text(xpos[i], hi + 0.003, f"\u03b2 = {b:.3f}{sig}\nP = {p:.4g}",
            ha="center", va="bottom", fontsize=5, linespacing=1.2)

# 分组分隔线
ax.axvline(3.15, color="#CCCCCC", lw=0.6, ls=":")
ax.text(1.825, 0.152, "All sites (n = %d)" % r["n_sites"], ha="center", fontsize=4.8, color="#888888")
ax.text(4.275, 0.152, "DO strata (n per subset)", ha="center", fontsize=4.8, color="#888888")

ax.set_xticks(xpos)
ax.set_xticklabels(xlabels, fontsize=5.5, linespacing=1.2)
ax.set_ylabel("\u03b2 (NO$_3^-$-N coefficient) (95% CI)", fontsize=7.5)
ax.set_xlim(-0.15, 5.45)
ax.set_ylim(-0.07, 0.17)
ax.text(-0.12, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False, width=0.5, length=2.5)
ax.tick_params(axis="y", right=False, labelright=False, width=0.5, length=2.5)

fig.savefig(OUT / "07i_edfig9b_model_coeff.png", dpi=300, facecolor="white")
fig.savefig(OUT / "07i_edfig9b_model_coeff.pdf", facecolor="white")
print("saved 07i_edfig9b_model_coeff.png/pdf (ED Fig 9b, 5-point dumbbell)")
plt.close()
