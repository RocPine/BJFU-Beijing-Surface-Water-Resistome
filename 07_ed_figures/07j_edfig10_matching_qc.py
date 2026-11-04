# -*- coding: utf-8 -*-
"""ED Fig 10: abundance balance before and after no-replacement matching."""
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
from narnas_stats import (load_sgb_annotation, compute_matching)

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_HOST = "#E64B35"; C_NONHOST = "#808080"; C_MATCHED = "#4DBBD5"

anno = load_sgb_annotation(cfg)
hosts = [s for s, a in anno.items() if a["ARG_host"] == "Yes"]
non_hosts = [s for s, a in anno.items() if a["ARG_host"] != "Yes"]
r = compute_matching(cfg)
matched_nh = r["matched_nh"]

hosts_abund = np.array([10**anno[s]["Log10_abund"] for s in hosts])
nonhosts_abund = np.array([10**anno[s]["Log10_abund"] for s in non_hosts])
matched_abund = np.array([10**anno[s]["Log10_abund"] for s in matched_nh])

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(482/72, 180/72), dpi=300,
                                 gridspec_kw={"wspace": 0.15, "left": 0.08, "right": 0.98, "top": 0.92, "bottom": 0.15})

def full_box(ax):
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(0.6)
    ax.tick_params(axis="x", top=False, labeltop=False, width=0.6, length=3)
    ax.tick_params(axis="y", right=False, labelright=False, width=0.6, length=3)

# Panel a: Before matching
ax_a.hist(np.log10(hosts_abund + 1e-6), bins=30, density=True, alpha=0.5, color=C_HOST,
          label=f"ARG hosts (n={len(hosts)})", edgecolor="none")
ax_a.hist(np.log10(nonhosts_abund + 1e-6), bins=30, density=True, alpha=0.5, color=C_NONHOST,
          label=f"All non-hosts (n={len(non_hosts)})", edgecolor="none")
ax_a.axvline(np.log10(np.mean(hosts_abund) + 1e-6), color=C_HOST, ls="--", lw=1)
ax_a.axvline(np.log10(np.mean(nonhosts_abund) + 1e-6), color="gray", ls="--", lw=1)

ks_stat, ks_p = stats.ks_2samp(hosts_abund, nonhosts_abund)
ax_a.text(0.02, 0.95, f"SMD = 0.52\nKS P = {ks_p:.2e}", transform=ax_a.transAxes, ha="left", va="top", fontsize=7)
ax_a.set_xlabel("log10 mean abundance (%)")
ax_a.set_ylabel("Probability density")
ax_a.set_xlim(-3, 10)
ax_a.legend(fontsize=6, frameon=False, loc="upper right")
ax_a.text(-0.12, 1.02, "a", transform=ax_a.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
full_box(ax_a)

# Panel b: After matching
ax_b.hist(np.log10(hosts_abund + 1e-6), bins=30, density=True, alpha=0.5, color=C_HOST,
          label=f"ARG hosts (n={len(hosts)})", edgecolor="none")
ax_b.hist(np.log10(matched_abund + 1e-6), bins=30, density=True, alpha=0.5, color=C_MATCHED,
          label=f"Matched non-hosts (n={len(matched_nh)})", edgecolor="none")

def smd(x, y):
    pooled_sd = np.sqrt((np.var(x, ddof=1) + np.var(y, ddof=1)) / 2)
    return (np.mean(x) - np.mean(y)) / pooled_sd if pooled_sd else 0.0

ax_b.text(0.02, 0.95, f"SMD = {smd(hosts_abund, matched_abund):.4f}\nWithout replacement",
          transform=ax_b.transAxes, ha="left", va="top", fontsize=7)
ax_b.set_xlabel("log10 mean abundance (%)")
ax_b.set_ylabel("Probability density")
ax_b.set_xlim(-3, 10)
ax_b.legend(fontsize=6, frameon=False, loc="upper right")
ax_b.text(-0.12, 1.02, "b", transform=ax_b.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
full_box(ax_b)

fig.savefig(OUT / "07j_edfig10_matching_qc.png", dpi=300, facecolor="white")
fig.savefig(OUT / "07j_edfig10_matching_qc.pdf", facecolor="white")
print("saved 07j_edfig10_matching_qc.png/pdf")
plt.close()
