# -*- coding: utf-8 -*-
"""ED Fig 9c: 氮形态特异性 (x=氮形态, y=Spearman rho vs 组合型宿主丰度, bootstrap CI)"""
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
from narnas_stats import load_sgb_annotation, load_site_data, load_site_env, site_abund, _split_hosts

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_FOCUS = "#E64B35"
C_OTHER_SIG = "#4DBBD5"
C_NS = "#9E9E9E"

anno = load_sgb_annotation(cfg)
hosts, _ = _split_hosts(anno)
combined = [s for s in hosts if anno[s]["N_narnas"] >= 1]

site_data = load_site_data(cfg)
sa = site_abund(site_data, combined)
env = load_site_env(cfg, ["TN", "NO2_N", "NH3_N"])
env["NO3_N"] = site_data["site_no3"]

FORMS = [("NO3_N", "NO$_3^-$-N"), ("TN", "TN"), ("NO2_N", "NO$_2^-$-N"), ("NH3_N", "NH$_3$-N")]

rng = np.random.default_rng(42)
rows = []
for key, label in FORMS:
    common = sorted(set(env[key]) & set(sa))
    x = np.array([env[key][s] for s in common])
    y = np.array([sa[s] for s in common])
    rho, p = stats.spearmanr(x, y)
    boot = []
    for _ in range(1000):
        idx = rng.integers(0, len(common), len(common))
        boot.append(stats.spearmanr(x[idx], y[idx])[0])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    rows.append((label, rho, lo, hi, p))

fig, ax = plt.subplots(figsize=(230/72, 200/72), dpi=300)
ax.set_box_aspect(1)

xpos = list(range(len(rows)))
for (label, rho, lo, hi, p), xp in zip(rows, xpos):
    if "NO$_3^-$" in label:
        c = C_FOCUS
    elif p < 0.05:
        c = C_OTHER_SIG
    else:
        c = C_NS
    ax.vlines(xp, lo, hi, color=c, lw=2.4, alpha=0.55, zorder=2)
    ax.scatter([xp, xp], [lo, hi], s=18, color=c, edgecolor="black", linewidth=0.3, zorder=4)
    ax.scatter([xp], [rho], s=22, marker="D", color="white", edgecolor=c, linewidth=0.7, zorder=5)
    sig = "" if p < 0.05 else " (n.s.)"
    ax.text(xp, hi + 0.02, f"\u03c1 = {rho:.3f}{sig}\nP = {p:.2g}",
            ha="center", va="bottom", fontsize=5, linespacing=1.2)

ax.set_xticks(xpos)
ax.set_xticklabels([r[0] for r in rows], fontsize=7)
ax.set_ylabel("Spearman \u03c1 (N form vs combined-host abundance)", fontsize=7.5)
ax.set_xlim(-0.5, len(rows) - 0.5)
ax.set_ylim(-0.25, 0.75)
ax.text(-0.12, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False, width=0.5, length=2.5)
ax.tick_params(axis="y", right=False, labelright=False, width=0.5, length=2.5)

fig.savefig(OUT / "07i_edfig9c_nform_specificity.png", dpi=300, facecolor="white")
fig.savefig(OUT / "07i_edfig9c_nform_specificity.pdf", facecolor="white")
print("saved 07i_edfig9c_nform_specificity.png/pdf (ED Fig 9c)")
plt.close()
