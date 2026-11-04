# -*- coding: utf-8 -*-
"""Fig 4b: HGT 频率 lollipop plot by ARG subtype
数据来源: SI Table 13 (HGT frequency: SampleID, ARG_subtype, HGT_abundance_standardized)"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

# 亚型 -> 类别映射 (Table 12)
wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)
ws12 = wb["Table 12"]
it = ws12.iter_rows(values_only=True)
next(it)
hdr12 = [str(h) if h else "" for h in next(it)]
i_sub = hdr12.index("ARG_Subtype")
i_type = hdr12.index("ARG_Type")
subtype_class = {}
for r in it:
    if r[0] is None: continue
    subtype_class[str(r[i_sub]).strip()] = str(r[i_type]).strip() if r[i_type] else "unknown"
# 频率表 (Table 13) 用短名 vanH-A/vanH-O, Table 12 用长名; 补特例 (与原总图一致)
subtype_class["vanH-A"] = "glycopeptide antibiotic"
subtype_class["vanH-O"] = "glycopeptide antibiotic"

# Table 13: HGT 频率
ws13 = wb["Table 13"]
it13 = ws13.iter_rows(values_only=True)
next(it13)
hdr13 = [str(h) if h else "" for h in next(it13)]
i_ab = hdr13.index("HGT_abundance_standardized")
rows = []
for r in it13:
    if r[0] is None: continue
    rows.append((str(r[1]).strip(), r[i_ab] if isinstance(r[i_ab], (int, float)) else 0.0))
wb.close()

freq = pd.DataFrame(rows, columns=["ARG_subtype", "HGT_abundance_standardized"])
sub_stats = freq.groupby("ARG_subtype")["HGT_abundance_standardized"].agg(total="sum")
sub_stats["total_pct"] = sub_stats["total"] * 100

type_order = {}
for s in sub_stats.index:
    t = subtype_class.get(s, "unknown")
    type_order.setdefault(t, []).append(s)
for t in type_order:
    type_order[t].sort(key=lambda s: -sub_stats["total_pct"][s])
type_tot = {t: sum(sub_stats["total_pct"][s] for s in ss) for t, ss in type_order.items()}
type_sorted = sorted(type_order.keys(), key=lambda t: -type_tot[t])

labels = []
for t in type_sorted:
    for s in type_order[t]:
        labels.append((t, s))

CLASS_COLORS = {
    "macrolide-lincosamide-streptogramin": "#E64B35", "aminoglycoside antibiotic": "#4DBBD5",
    "glycopeptide antibiotic": "#00A087", "multidrug": "#3C5488",
    "polymyxin": "#F39B7F", "mupirocin": "#8491B4", "pleuromutilin_tiamulin": "#7E6148",
}

def short_name(s):
    return str(s).replace("vanH gene in ", "").replace(" antibiotic", "").replace("_tiamulin", "")

# 尺寸/字体复刻原总图 (plot_fig4_v50.py): b 面板 0.38x0.28 于 183x220mm, 全框
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

for i, (t, s) in enumerate(labels):
    c = CLASS_COLORS.get(t, "#8491B4")
    val_pct = sub_stats["total_pct"][s]
    ax.hlines(i, 0, val_pct, color=c, lw=1.0)
    ax.scatter(val_pct, i, s=25, color=c, zorder=3)

ax.axvline(0, color="gray", lw=0.5)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels([short_name(s) for _, s in labels], fontsize=7)
ax.set_xlabel("Total HGT frequency (%)")
ax.set_xlim(0, sub_stats["total_pct"].max() * 1.2)
ax.set_title("b", fontsize=11, fontweight="bold", loc="left")

# 全框 (与原总图一致: 默认 4 边可见, 无 top/right 刻度线)

fig.savefig(OUT / "04b_hgt_lollipop.png", dpi=300, facecolor="white")
fig.savefig(OUT / "04b_hgt_lollipop.pdf", facecolor="white")
print("saved 04b_hgt_lollipop.png/pdf (frequency)")
top = sub_stats.sort_values("total_pct", ascending=False).head(3)
for s, v in top.iterrows():
    print(f"  {s}: {v['total_pct']:.3f}%")
plt.close()