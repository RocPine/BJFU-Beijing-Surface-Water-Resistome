# -*- coding: utf-8 -*-
"""Fig 4c: HGT subtype × donor phylum bubble plot"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import pandas as pd
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

hgt = data["hgt"]
subtype_class = data["subtype_class"]
short_name = data["short_name"]
CLASS_COLORS = data["CLASS_COLORS"]
donor_phys = data["donor_phys"]

# ---- 数据准备 ----
cross_sub_phy = pd.crosstab(hgt["ARG_Subtype"], hgt["donor_p"])
top_subs = hgt["ARG_Subtype"].value_counts()
maj_subs = top_subs[top_subs >= 2].index.tolist()
maj_phys_d = donor_phys.sort_values(ascending=False).head(6).index.tolist()

ypos_d = np.arange(len(maj_subs))[::-1]
xpos_d = np.arange(len(maj_phys_d))
S_SQRT = 60

# ---- 绘图 ----
# 尺寸/字体复刻原总图 (plot_fig4_v50.py): c 面板 0.38x0.28 于 183x220mm, 全框
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

for i, sub in enumerate(maj_subs):
    for j, phy in enumerate(maj_phys_d):
        cnt = cross_sub_phy.loc[sub, phy] if sub in cross_sub_phy.index and phy in cross_sub_phy.columns else 0
        if cnt > 0:
            c = CLASS_COLORS.get(subtype_class.get(sub, ""), "#8491B4")
            ax.scatter(j, ypos_d[i], s=S_SQRT * np.sqrt(cnt), color=c, alpha=0.75,
                       edgecolor="white", linewidth=0.5)
            if cnt >= 5:
                ax.text(j, ypos_d[i], str(int(cnt)), ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")

ax.set_xticks(xpos_d)
ax.set_xticklabels([p.replace("p__", "") for p in maj_phys_d], rotation=45, ha="right", fontsize=8)
ax.set_yticks(ypos_d)
ax.set_yticklabels([short_name(s) for s in maj_subs], fontsize=7)
ax.set_xlim(-0.8, len(maj_phys_d) - 0.2)
ax.set_ylim(-0.8, len(maj_subs) - 0.2)
ax.set_xlabel("Donor phylum")
ax.set_ylabel("ARG subtype")
ax.set_title("c", fontsize=11, fontweight="bold", loc="left")

# 全框 (与原总图一致)

fig.savefig(OUT / "04c_hgt_bubble.png", dpi=300, facecolor="white")
fig.savefig(OUT / "04c_hgt_bubble.pdf", facecolor="white")
print("saved 04c_hgt_bubble.png/pdf")
plt.close()