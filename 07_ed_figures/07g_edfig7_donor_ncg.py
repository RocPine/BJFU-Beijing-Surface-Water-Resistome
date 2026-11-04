# -*- coding: utf-8 -*-
"""ED Fig 7: 12 transferred ARG subtypes donor Mean NCGs (N_total)
数据来源: 仅发布版 SI 表 (Table 12 HGT 事件, Table 9 MAG->SGB, Table 15 SGB 注释)"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

CLASS_COLORS = {
    "macrolide-lincosamide-streptogramin": "#E64B35",
    "aminoglycoside antibiotic": "#4DBBD5",
    "glycopeptide antibiotic": "#00A087",
    "multidrug": "#3C5488",
    "polymyxin": "#F39B7F",
    "mupirocin": "#8491B4",
    "pleuromutilin_tiamulin": "#7E6148",
}
TYPE_SHORT = {
    "macrolide-lincosamide-streptogramin": "MLS",
    "aminoglycoside antibiotic": "Aminoglycoside",
    "glycopeptide antibiotic": "Glycopeptide",
    "multidrug": "Multidrug",
    "polymyxin": "Polymyxin",
    "mupirocin": "Mupirocin",
    "pleuromutilin_tiamulin": "Pleuromutilin",
}

wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)

# ---------- Table 15: SGB 注释 (ARG_host, N_total_NCGs) ----------
ws15 = wb["Table 15"]
it15 = ws15.iter_rows(values_only=True)
next(it15)
hdr15 = [str(h) if h else "" for h in next(it15)]
i_host = hdr15.index("ARG_host")
i_nt = hdr15.index("N_total_NCGs")
sgb_anno = {}
for row in it15:
    if row[0] is None: continue
    sgb_anno[str(row[0]).strip()] = {
        "ARG_host": row[i_host],
        "N_total": row[i_nt] if isinstance(row[i_nt], (int, float)) else 0,
    }

# ---------- Table 9: MAG -> SGB ----------
ws9 = wb["Table 9"]
it9 = ws9.iter_rows(values_only=True)
next(it9)
hdr9 = [str(h) if h else "" for h in next(it9)]
i_m = hdr9.index("Member_MAG")
mag2sgb = {}
for row in it9:
    if row[0] is None: continue
    sgb = str(row[0]).strip()
    if row[i_m]:
        for mag in str(row[i_m]).split(";"):
            mag = mag.strip()
            if mag: mag2sgb[mag] = sgb
sgb_names_set = set(sgb_anno)

# ---------- Table 12: HGT 事件 ----------
ws12 = wb["Table 12"]
it12 = ws12.iter_rows(values_only=True)
next(it12)
hdr12 = [str(h) if h else "" for h in next(it12)]
i_dmag = hdr12.index("donor_mag")
i_sub = hdr12.index("ARG_Subtype")
i_type = hdr12.index("ARG_Type")
wb.close()

donor_subtypes = defaultdict(list)
for row in it12:
    if row[0] is None: continue
    d_mag = str(row[i_dmag]).strip()
    subtype = row[i_sub]
    arg_type = row[i_type]
    # MAG -> SGB (与 Fig 4 loader 相同规则)
    d_sgb = None
    if d_mag in sgb_names_set:
        d_sgb = d_mag
    elif d_mag in mag2sgb:
        d_sgb = mag2sgb[d_mag]
    if d_sgb and d_sgb in sgb_anno and sgb_anno[d_sgb]["ARG_host"] == "Yes":
        donor_subtypes[subtype].append({
            "N_total": sgb_anno[d_sgb]["N_total"],
            "ARG_Type": arg_type,
        })

# 短亚型名
def short_name(s):
    return str(s).replace("vanH gene in ", "vanH-").replace(" antibiotic", "").replace("_tiamulin", "")

subtype_stats = []
for st, lst in donor_subtypes.items():
    vals = [x["N_total"] for x in lst]
    mean = np.mean(vals)
    boots = [np.mean(np.random.choice(vals, len(vals), replace=True)) for _ in range(5000)] if len(vals) > 1 else [mean]
    subtype_stats.append({
        "subtype": st,
        "short": short_name(st),
        "n": len(lst),
        "mean": mean,
        "ci_low": np.percentile(boots, 2.5),
        "ci_high": np.percentile(boots, 97.5),
        "type": lst[0]["ARG_Type"],
    })

subtype_stats.sort(key=lambda x: -x["mean"])

# ============ Figure ============
fig, ax = plt.subplots(figsize=(183/25.4, 75/25.4), dpi=300)
x_pos = np.arange(len(subtype_stats))
BLACK = "#000000"

for i, s in enumerate(subtype_stats):
    c = CLASS_COLORS.get(s["type"], "#8491B4")
    if s["n"] >= 3:
        ax.plot([i, i], [s["ci_low"], s["ci_high"]], color=BLACK, lw=1.5, alpha=0.8, zorder=3)
        ax.scatter([i, i], [s["ci_low"], s["ci_high"]], color=BLACK, s=15,
                   edgecolor="white", linewidth=0.4, zorder=5)
    ax.scatter([i], [s["mean"]], s=40, color=c, edgecolor="black", linewidth=0.5, zorder=5)
    ax.text(i, s["mean"] + 0.5, f"n={s['n']}", ha="center", fontsize=6.5, color="#555555")

ax.set_xticks(x_pos)
ax.set_xticklabels([s["short"] for s in subtype_stats], rotation=45, ha="right", fontsize=6.5)
ax.set_ylabel("Number of NCGs in donor ARG hosts", fontsize=8)
ax.set_ylim(0, 27)
ax.set_xlim(-0.5, len(subtype_stats) - 0.5)

# ED Fig 7 为独立单面板图, 无 panel letter

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.6)
ax.tick_params(axis="x", width=0.5, length=3, direction="out")
ax.tick_params(axis="y", width=0.5, length=3, direction="out")

handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=CLASS_COLORS.get(t, "#8491B4"),
                  markersize=6, label=TYPE_SHORT.get(t, t)) for t in sorted(set(s["type"] for s in subtype_stats))]
ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right", ncol=2)

fig.savefig(OUT / "07g_edfig7_donor_ncg.png", dpi=300, facecolor="white")
fig.savefig(OUT / "07g_edfig7_donor_ncg.pdf", facecolor="white")
print("Saved 07g_edfig7_donor_ncg.png/pdf")
for s in subtype_stats:
    print(f"  {s['subtype']} ({s['short']}): mean={s['mean']:.1f} [{s['ci_low']:.1f}, {s['ci_high']:.1f}], n={s['n']}, {s['type']}")
plt.close()