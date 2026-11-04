# -*- coding: utf-8 -*-
"""Fig 6a: 三组梯度散点图 (NAR/NAP/NAS-carrying hosts vs hosts without vs non-hosts)"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "08_core_statistics"))
from narnas_stats import compute_fig6a

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_HIGH = "#E64B35"; C_LOW = "#4DBBD5"; C_NONHOST = "#00A087"

r = compute_fig6a(cfg)

fig, ax = plt.subplots(figsize=(230/72, 176/72), dpi=300)

def plot_group(x, y, c, label, ls="-"):
    ax.scatter(x, y, s=28, color=c, alpha=0.55, edgecolor="none", zorder=3, label=label)
    z = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    fit = np.polyval(z, xs)
    n = len(x)
    residuals = y - np.polyval(z, x)
    s_err = np.sqrt(np.sum(residuals**2) / (n - 2))
    x_bar = np.mean(x)
    sxx = np.sum((x - x_bar)**2)
    t_crit = stats.t.ppf(0.975, n - 2)
    band = t_crit * s_err * np.sqrt(1 / n + (xs - x_bar)**2 / sxx)
    ax.fill_between(xs, fit - band, fit + band, color=c, alpha=0.14, linewidth=0, zorder=1)
    ax.plot(xs, fit, color=c, lw=1.6, ls=ls, zorder=4)

plot_group(r["x_h"], r["y_h"], C_HIGH, "ARG hosts with NAR/NAP/NAS")
plot_group(r["x_l"], r["y_l"], C_LOW, "ARG hosts without NAR/NAP/NAS", ls="--")
plot_group(r["x_n"], r["y_n"], C_NONHOST, "Non-hosts")

ax.set_xlabel("NO3-N (mg/L)")
ax.set_ylabel("Relative abundance (%)")
ax.set_xlim(0, 10)
ax.set_ylim(0, 100)

ax.text(0.03, 0.94, f"ρ = {r['rho_h']:.2f}, P = {r['p_h']:.1e}", transform=ax.transAxes,
        ha="left", va="top", fontsize=7.5, color=C_HIGH)
ax.text(0.03, 0.83, f"ρ = {r['rho_l']:.2f}, P = {r['p_l']:.2g} (n.s.)", transform=ax.transAxes,
        ha="left", va="top", fontsize=7.5, color=C_LOW)
ax.text(0.03, 0.72, f"ρ = {r['rho_n']:.2f}, P = {r['p_n']:.1e}", transform=ax.transAxes,
        ha="left", va="top", fontsize=7.5, color=C_NONHOST)
ax.text(-0.12, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

ax.legend(fontsize=7, frameon=False, loc="upper left")
for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False)
ax.tick_params(axis="y", right=False, labelright=False)

fig.savefig(OUT / "06a_three_group_gradient.png", dpi=300, facecolor="white")
fig.savefig(OUT / "06a_three_group_gradient.pdf", facecolor="white")
print("saved 06a_three_group_gradient.png/pdf")
plt.close()
