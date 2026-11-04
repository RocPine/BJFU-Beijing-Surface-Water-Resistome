# -*- coding: utf-8 -*-
"""Fig 4a: HGT Sankey diagram (donor vs recipient phyla)"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from fig4_data_loader import load_fig4_data

# ---- 加载配置 ----
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

# ---- 解包数据 ----
hgt = data["hgt"]
cross_pp = data["cross_pp"]
node_order = data["node_order"]
donor_phys = data["donor_phys"]
recv_phys = data["recv_phys"]
CLASS_COLORS = data["CLASS_COLORS"]
PHYLA_COLORS = data["PHYLA_COLORS"]
TYPE_SHORT = data["TYPE_SHORT"]

# ---- Sankey 计算 ----
donor_phys_s = donor_phys.reindex(node_order).fillna(0)
recv_phys_s = recv_phys.reindex(node_order).fillna(0)

MAX_BLOCK = 0.35
max_val = max(donor_phys_s.max(), recv_phys_s.max())
scale = MAX_BLOCK / max_val if max_val > 0 else 1
d_h_dict = {p: donor_phys_s.get(p, 0) * scale for p in node_order}
r_h_dict = {p: recv_phys_s.get(p, 0) * scale for p in node_order}

d_heights = [d_h_dict[p] for p in node_order]
d_bottoms = np.concatenate([[0.0], np.cumsum(d_heights)[:-1]])
r_heights = [r_h_dict[p] for p in node_order]
r_bottoms = np.concatenate([[0.0], np.cumsum(r_heights)[:-1]])

flows = []
for d in node_order:
    if d not in cross_pp.index: continue
    for r in node_order:
        if r not in cross_pp.columns: continue
        cnt = cross_pp.loc[d, r]
        if cnt > 0:
            sub = hgt[(hgt["donor_phy"] == d) & (hgt["recipient_phy"] == r)]
            main_type = sub["ARG_Type"].mode()[0] if len(sub) > 0 else "unknown"
            flows.append((d, r, int(cnt), main_type))

STACK_PADDING = 0.9
def stack_flows(node_flows, total_h):
    total_cnt = sum(c for c, _ in node_flows)
    stack_h = total_h * STACK_PADDING
    node_flows = sorted(node_flows, key=lambda x: -x[0])
    centers = {}
    y_start = (total_h - stack_h) / 2
    y_cursor = y_start
    for cnt, fidx in node_flows:
        h = cnt / total_cnt * stack_h
        centers[fidx] = y_cursor + h / 2
        y_cursor += h
    return centers

d_flow_centers = {d: {} for d in node_order}
r_flow_centers = {r: {} for r in node_order}
for fidx, (d, r, cnt, mtype) in enumerate(flows):
    d_flow_centers[d][fidx] = None
    r_flow_centers[r][fidx] = None

for d in node_order:
    d_fids = list(d_flow_centers[d].keys())
    d_cnts = [flows[f][2] for f in d_fids]
    if d_fids:
        centers = stack_flows(list(zip(d_cnts, d_fids)), d_h_dict[d])
        for fidx, yc in centers.items():
            d_flow_centers[d][fidx] = yc

for r in node_order:
    r_fids = list(r_flow_centers[r].keys())
    r_cnts = [flows[f][2] for f in r_fids]
    if r_fids:
        centers = stack_flows(list(zip(r_cnts, r_fids)), r_h_dict[r])
        for fidx, yc in centers.items():
            r_flow_centers[r][fidx] = yc

d_total_h = sum(d_heights)
r_total_h = sum(r_heights)
y_max = max(d_total_h, r_total_h)

# ---- 绘图 ----
# 尺寸复刻原总图 (plot_fig4_v50.py): a 面板 0.90x0.28 于 183x220mm
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
fig = plt.figure(figsize=(164.7/25.4, 61.6/25.4), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

x_d = 0.05
x_r = 0.62
sankey_left = 0.02
sankey_right = 0.65

for fidx, (d, r, cnt, mtype) in enumerate(flows):
    y1_center = d_flow_centers[d][fidx]
    y2_center = r_flow_centers[r][fidx]
    w_half = cnt / (donor_phys_s[d] if donor_phys_s[d] > 0 else 1) * d_h_dict[d] / 2
    color = CLASS_COLORS.get(mtype, "#A6CEE3")
    alpha = 0.5
    y1 = (d_bottoms[node_order.index(d)] + y1_center) / y_max
    y2 = (r_bottoms[node_order.index(r)] + y2_center) / y_max
    h_half = w_half / y_max
    t = np.linspace(0, 1, 40)
    x0, y0 = x_d, y1
    x3, y3 = x_r, y2
    ctrl_x1 = x0 + (x3 - x0) * 0.45
    ctrl_x2 = x0 + (x3 - x0) * 0.55
    bx = (1-t)**3*x0 + 3*(1-t)**2*t*ctrl_x1 + 3*(1-t)*t**2*ctrl_x2 + t**3*x3
    by = (1-t)**3*y0 + 3*(1-t)**2*t*y1 + 3*(1-t)*t**2*y2 + t**3*y3
    w_half_d = cnt / donor_phys_s[d] * d_h_dict[d] * STACK_PADDING / 2 / y_max
    w_half_r = cnt / recv_phys_s[r] * r_h_dict[r] * STACK_PADDING / 2 / y_max
    w_half_t = w_half_d + (w_half_r - w_half_d) * t
    by_top = by + w_half_t
    by_bot = by - w_half_t
    poly_pts = list(zip(bx, by_top)) + list(zip(bx[::-1], by_bot[::-1]))
    ax.add_patch(Polygon(poly_pts, closed=True, facecolor=color, alpha=alpha,
                       edgecolor="none", zorder=2 if cnt >= 5 else 1))

for i, p in enumerate(node_order):
    c = PHYLA_COLORS.get(p, "#666666")
    h_d = d_h_dict[p] / y_max
    h_r = r_h_dict[p] / y_max
    if h_d > 0.003:
        yb_d = d_bottoms[i] / y_max
        ax.add_patch(Rectangle((sankey_left, yb_d), x_d - sankey_left, h_d,
                       facecolor=c, alpha=0.95, edgecolor="none", zorder=3))
    if h_r > 0.003:
        yb_r = r_bottoms[i] / y_max
        ax.add_patch(Rectangle((x_r, yb_r), sankey_right - x_r, h_r,
                       facecolor=c, alpha=0.85, edgecolor="none", zorder=3))

diag_self = sum(cross_pp.loc[d, d] if (d in cross_pp.index and d in cross_pp.columns) else 0 for d in node_order)
ax.text(x_d/2, 0.97, "Donor", ha="center", fontsize=10, fontweight="bold")
ax.text((x_r + sankey_right)/2, 0.97, "Recipient", ha="center", fontsize=10, fontweight="bold")
ax.set_title("a", fontsize=11, fontweight="bold", loc="left")
ax.text(0.02, 0.02, f"Intra-phylum: {diag_self}/{len(hgt)} events", transform=ax.transAxes, fontsize=7, color="gray")

# ---- 图例 (单独保存或内嵌) ----
# Phylum legend
phylum_items = [(p, PHYLA_COLORS.get(p, "#666666")) for p in node_order if (donor_phys.get(p,0)+recv_phys.get(p,0)) > 0]
legend_left = 0.67
legend_y_min = 0.05
legend_y_max = 0.95
legend_height = legend_y_max - legend_y_min
item_height = legend_height / max(len(phylum_items), 1)

for i, (phy, color) in enumerate(phylum_items):
    y = legend_y_max - (i + 0.5) * item_height
    block_height = item_height * 0.8
    ax.add_patch(Rectangle((legend_left, y - block_height/2), 0.025, block_height,
                   facecolor=color, edgecolor="none", zorder=5, transform=ax.transAxes))
    ax.text(legend_left + 0.03, y, phy.replace("p__", ""),
            fontsize=7, va="center", ha="left", transform=ax.transAxes, zorder=5)

# ARG class legend
class_legend_left = legend_left + 0.11
class_items = [(t, CLASS_COLORS[t]) for t in CLASS_COLORS]
for i, (cls, color) in enumerate(class_items):
    y = legend_y_max - (i + 0.5) * item_height
    block_height = item_height * 0.8
    ax.add_patch(Rectangle((legend_left + 0.10, y - block_height/2), 0.025, block_height,
                   facecolor=color, edgecolor="none", zorder=5, transform=ax.transAxes))
    ax.text(legend_left + 0.13, y, TYPE_SHORT.get(cls, cls),
            fontsize=7, va="center", ha="left", transform=ax.transAxes, zorder=5)

ax.set_title("a", fontsize=11, fontweight="bold", loc="left")
ax.text(0.02, 0.02, f"Intra-phylum: {diag_self}/{len(hgt)} events", transform=ax.transAxes, fontsize=7, color="gray")

# ---- 保存 ----
fig.savefig(OUT / "04a_hgt_sankey.png", dpi=300, facecolor="white")
fig.savefig(OUT / "04a_hgt_sankey.pdf", facecolor="white")
print("saved 04a_hgt_sankey.png/pdf")

# 单独保存 ARG class 图例
from matplotlib.lines import Line2D
fig_leg = plt.figure(figsize=(180/25.4, 20/25.4), dpi=300)
legend_handles_arg = [Line2D([0],[0], marker="o", color="w", markerfacecolor=c, markersize=10, label=TYPE_SHORT.get(lab, lab))
                      for lab, c in CLASS_COLORS.items()]
fig_leg.legend(handles=legend_handles_arg, loc="center",
               fontsize=11, frameon=True, ncol=4, title="ARG class", title_fontsize=12,
               edgecolor="#CCCCCC", columnspacing=1.5, handletextpad=0.5)
fig_leg.savefig(OUT / "Fig4_ARGclass_legend.png", dpi=300, facecolor="white")
fig_leg.savefig(OUT / "Fig4_ARGclass_legend.pdf", facecolor="white")
print("saved ARG class legend")
plt.close('all')