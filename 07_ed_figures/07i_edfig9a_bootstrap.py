# -*- coding: utf-8 -*-
"""ED Fig 9a (Fig 6a 梯度关联的 bootstrap 验证): 三组 bootstrap ρ 的 KDE

ED Fig 9 = Fig 6 梯度关联的稳健性与特异性: a=bootstrap ρ KDE (本脚本),
b=模型系数 + DO 分层 (07i_edfig9b_model_coeff.py),
c=氮形态特异性 (07i_edfig9c_nform_specificity.py), 脚本名含 9a 前缀"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "08_core_statistics"))
from narnas_stats import (load_sgb_annotation, load_site_data, site_abund,
                          site_spearman, _split_hosts)

with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
for k, v in cfg.items():
    if isinstance(v, str) and not pathlib.Path(v).is_absolute():
        cfg[k] = str(ROOT / v)

C_HIGH = "#E64B35"; C_LOW = "#4DBBD5"; C_NONHOST = "#00A087"

anno = load_sgb_annotation(cfg)
hosts, non_hosts = _split_hosts(anno)
high = [s for s in hosts if anno[s]["N_narnas"] >= 1]
low = [s for s in hosts if anno[s]["N_narnas"] < 1]

site_data = load_site_data(cfg)
sa_high = site_abund(site_data, high)
sa_low = site_abund(site_data, low)
sa_nonhost = site_abund(site_data, non_hosts)

common_sites = sorted(set(site_data["site_no3"]) & set(sa_high) & set(sa_low) & set(sa_nonhost))
x_no3 = np.array([site_data["site_no3"][s] for s in common_sites])
y_h = np.array([sa_high[s] for s in common_sites])
y_l = np.array([sa_low[s] for s in common_sites])
y_n = np.array([sa_nonhost[s] for s in common_sites])

rng = np.random.default_rng(42)
n_boot = 1000
boot_h = np.zeros(n_boot); boot_l = np.zeros(n_boot); boot_n = np.zeros(n_boot)
n_hl = 0; n_hn = 0
for b in range(n_boot):
    idx = rng.integers(0, len(common_sites), len(common_sites))
    r_h = stats.spearmanr(x_no3[idx], y_h[idx])[0]
    r_l = stats.spearmanr(x_no3[idx], y_l[idx])[0]
    r_n = stats.spearmanr(x_no3[idx], y_n[idx])[0]
    boot_h[b] = r_h; boot_l[b] = r_l; boot_n[b] = r_n
    if r_h > r_l: n_hl += 1
    if r_h > r_n: n_hn += 1

fig, ax = plt.subplots(figsize=(326/72, 255/72), dpi=300)
x_grid = np.linspace(-0.8, 0.8, 400)
obs_h, obs_l, obs_n = 0.437, 0.146, -0.429

groups = [
    (boot_h, C_HIGH, f"ARG hosts with NAR/NAP/NAS (n={len(high)})", obs_h),
    (boot_l, C_LOW, f"ARG hosts without NAR/NAP/NAS (n={len(low)})", obs_l),
    (boot_n, C_NONHOST, f"Non-hosts (n={len(non_hosts)})", obs_n),
]

for data, c, label, obs in groups:
    kde = gaussian_kde(data, bw_method=0.25)
    dens = kde(x_grid)
    ax.plot(x_grid, dens, color=c, lw=1.5, label=label)
    ax.axvline(obs, color=c, lw=0.8, ls="--", alpha=0.7)

ax.set_xlabel("Spearman ρ")
ax.set_ylabel("Probability density")
ax.set_xlim(-0.8, 0.8)
ax.set_ylim(0, None)

ax.text(0.02, 0.02,
        f"Bootstrap proportion (High > Non-host) = {n_hn/n_boot:.3f}\n"
        f"Bootstrap proportion (High > Low) = {n_hl/n_boot:.3f}",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=7, linespacing=1.5)
ax.legend(fontsize=7, frameon=False, loc="upper right")
ax.text(-0.12, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

for spine in ["top", "right", "left", "bottom"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(0.7)
ax.tick_params(axis="x", top=False, labeltop=False, width=0.5, length=3)
ax.tick_params(axis="y", right=False, labelright=False, width=0.5, length=3)

fig.savefig(OUT / "07i_edfig9a_bootstrap.png", dpi=300, facecolor="white")
fig.savefig(OUT / "07i_edfig9a_bootstrap.pdf", facecolor="white")
print("saved 07i_edfig9a_bootstrap.png/pdf")
plt.close()
