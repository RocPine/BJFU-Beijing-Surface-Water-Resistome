# -*- coding: utf-8 -*-
"""Fig 6d: NAR/NAP/NAS gene plot using matched-pair McNemar results."""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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

r = compute_matching(cfg)
gene_results = r["gene_results"]
bonf = r["bonf"]

# 基因 -> 类型映射 (NCycDB gene_type)
GENE_GTYPE = {
    "nar": "NAR", "narb": "NAR", "narc": "NAR", "narg": "NAR", "narh": "NAR",
    "nari": "NAR", "narj": "NAR", "nary": "NAR", "narz": "NAR",
    "napa": "NAP", "napb": "NAP", "napc": "NAP", "ansb": "NAP", "gdh": "NAP",
    "glsa": "NAP", "gs": "NAP", "nmo": "NAP",
    "nasa": "NAS", "nasb": "NAS",
}
GTYPE_COLORS = {"NAR": "#D62728", "NAP": "#FF7F0E", "NAS": "#1F77B4"}
GTYPE_SHORT = {"NAR": "NAR (dissimilatory)", "NAP": "NAP (dissimilatory)", "NAS": "NAS (assimilatory)"}
GENE_DISPLAY = {
    "nar": "nar", "narb": "narB", "narc": "narC", "narg": "narG", "narh": "narH",
    "nari": "narI", "narj": "narJ", "nary": "narY", "narz": "narZ",
    "napa": "napA", "napb": "napB", "napc": "napC", "ansb": "ansB",
    "gdh": "gdh", "glsa": "glsA", "gs": "gs", "nmo": "nmo",
    "nasa": "nasA", "nasb": "nasB",
}

genes = [r["gene"] for r in gene_results]
ors = [r["or"] for r in gene_results]
ps = [r["p"] for r in gene_results]
x_v = [np.log2(o) for o in ors]
p_adj = [min(p * len(gene_results), 1.0) for p in ps]
y_v = [-np.log10(max(p, np.finfo(float).tiny)) for p in p_adj]

C_GREY = "#B0B0B0"

fig, ax = plt.subplots(figsize=(230/72, 200/72), dpi=300)

for i, (g, x, y) in enumerate(zip(genes, x_v, y_v)):
    gtype = GENE_GTYPE.get(g, "NAR")
    c = GTYPE_COLORS[gtype] if ps[i] < bonf else C_GREY
    ax.scatter(x, y, s=32, color=c, alpha=0.85 if ps[i] < bonf else 0.6,
               edgecolor="none", linewidth=0, zorder=3)

ax.axhline(-np.log10(bonf), color="#999999", lw=0.7, ls="--")
ax.text(ax.get_xlim()[0] + 0.1, -np.log10(bonf) + 0.8,
        f"Bonferroni adjusted P = 0.05", fontsize=6, color="#666666")
ax.axvline(0, color="black", lw=0.7)
for g, x, y in zip(genes, x_v, y_v):
    ax.annotate(GENE_DISPLAY.get(g, g), (x, y), xytext=(x + 0.08, y + 1.0),
                fontsize=6, ha="left", va="center", fontstyle="italic")

ax.set_xlabel("log2(OR) (ARG hosts vs matched non-hosts)", fontsize=8)
ax.set_ylabel("-log10(adj. P)", fontsize=8)
ax.set_xlim(-3.5, 3.5)
handles_d = [Line2D([0], [0], marker="o", color="w", markerfacecolor=GTYPE_COLORS[t],
                    markersize=6, label=GTYPE_SHORT[t]) for t in ["NAR", "NAP", "NAS"]]
ax.legend(handles=handles_d, loc="upper left", fontsize=6.5, frameon=False)
ax.text(-0.12, 1.02, "d", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False)
ax.tick_params(axis="y", right=False, labelright=False)

fig.savefig(OUT / "06d_narnas_volcano.png", dpi=300, facecolor="white")
fig.savefig(OUT / "06d_narnas_volcano.pdf", facecolor="white")
print("saved 06d_narnas_volcano.png/pdf")
plt.close()
