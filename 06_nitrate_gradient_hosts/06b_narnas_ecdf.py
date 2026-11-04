# -*- coding: utf-8 -*-
"""Fig 6b (原 Fig 5e): NAR/NAP/NAS 基因数 per SGB 的 ECDF (ARG hosts vs non-hosts)

尺寸/样式与 06 系列面板 (06a/06c/06d) 一致: figsize 230x200 pt, 面板字母 b
数据来源: SI Table 15
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
from narnas_stats import compute_fig5e

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_HOST = "#E64B35"
C_NONHOST = "#7F7F7F"
C_EDGE = "#333333"

r = compute_fig5e(cfg)
h_counts = r["host_counts"]
nh_counts = r["nonhost_counts"]
n_hosts = r["n_hosts"]
n_nonhosts = r["n_nonhosts"]
h_pct = r["carriage_host_pct"]
nh_pct = r["carriage_nonhost_pct"]

fig, ax = plt.subplots(figsize=(230/72, 200/72), dpi=300)

def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(data)+1) / len(data) * 100
    return x, y

x_h, y_h = ecdf(h_counts)
x_n, y_n = ecdf(nh_counts)

ax.step(x_h, y_h, where="post", color=C_HOST, lw=1.8, label=f"ARG hosts (n = {n_hosts})")
ax.step(x_n, y_n, where="post", color=C_NONHOST, lw=1.8, ls="--", label=f"Non-hosts (n = {n_nonhosts})")

ax.axvline(1, color=C_EDGE, lw=0.7, ls=":", alpha=0.6)
ax.scatter([1], [h_pct], s=32, color=C_HOST, zorder=6, edgecolor="white", linewidth=0.6)
ax.scatter([1], [nh_pct], s=32, color=C_NONHOST, zorder=6, edgecolor="white", linewidth=0.6)
ax.text(1.5, h_pct, f"{h_pct:.1f}%", fontsize=8, va="center", color=C_HOST, fontweight="bold")
ax.text(1.5, nh_pct, f"{nh_pct:.1f}%", fontsize=8, va="center", color=C_NONHOST)
ax.text(1.5, 2.5, "\u22651 gene", fontsize=6.5, va="bottom", color=C_EDGE, style="italic")

ax.text(0.97, 0.045,
        f"Fisher exact OR = {r['or']:.2f}\nP = {r['p_fisher']:.1e}\nMann-Whitney P = {r['p_mw']:.1e}",
        transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.9, edgecolor="#CCCCCC", linewidth=0.5))

ax.set_xlabel("Nitrate reduction genes (NAR/NAP/NAS) per SGB", fontsize=9)
ax.set_ylabel("Cumulative proportion of SGBs (%)", fontsize=9)
ax.set_xlim(-0.4, 16)
ax.set_ylim(0, 103)
ax.tick_params(axis="x", top=False, labeltop=False)
ax.tick_params(axis="y", right=False, labelright=False)
ax.legend(fontsize=8, frameon=False, loc="lower right")
ax.text(-0.12, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)

fig.savefig(OUT / "06b_narnas_ecdf.png", dpi=300, facecolor="white")
fig.savefig(OUT / "06b_narnas_ecdf.pdf", facecolor="white")
print("saved 06b_narnas_ecdf.png/pdf (Fig 6b panel)")
plt.close()

