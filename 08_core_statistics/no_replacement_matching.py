"""No-replacement abundance matching sensitivity analysis for Fig. 6c/6d.

Uses globally optimal one-to-one matching on log10 mean SGB abundance.
This script writes new results only and does not modify the existing analysis.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest, fisher_exact, mannwhitneyu, wilcoxon

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "no_replacement_matching"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "08_core_statistics"))
from narnas_stats import NARNAS_GENES, _parse_genes, _split_hosts, load_config, load_sgb_annotation  # noqa: E402


def smd(x: np.ndarray, y: np.ndarray) -> float:
    pooled_sd = np.sqrt((np.var(x, ddof=1) + np.var(y, ddof=1)) / 2)
    return float((np.mean(x) - np.mean(y)) / pooled_sd) if pooled_sd else 0.0


def write_csv(path: pathlib.Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    anno = load_sgb_annotation(load_config())
    hosts, non_hosts = _split_hosts(anno)
    host_values = np.array([anno[s]["Log10_abund"] for s in hosts], dtype=float)
    nonhost_values = np.array([anno[s]["Log10_abund"] for s in non_hosts], dtype=float)
    if len(hosts) > len(non_hosts):
        raise ValueError("No-replacement matching requires at least as many non-hosts as hosts")

    # Global assignment prevents arbitrary dependence on host ordering.
    cost = np.abs(host_values[:, None] - nonhost_values[None, :])
    host_idx, nonhost_idx = linear_sum_assignment(cost)
    if len(host_idx) != len(hosts):
        raise RuntimeError("Assignment did not produce one control for every host")
    matched_hosts = [hosts[i] for i in host_idx]
    matched_nonhosts = [non_hosts[i] for i in nonhost_idx]
    distances = cost[host_idx, nonhost_idx]

    host_counts = np.array([anno[s]["N_narnas"] for s in matched_hosts], dtype=float)
    nonhost_counts = np.array([anno[s]["N_narnas"] for s in matched_nonhosts], dtype=float)
    w_stat, w_p = wilcoxon(host_counts, nonhost_counts, alternative="two-sided")
    mw_stat, mw_p = mannwhitneyu(host_counts, nonhost_counts, alternative="two-sided")

    genes = sorted(
        set().union(*[_parse_genes(anno[s]["Specific_genes"]) for s in matched_hosts])
        & NARNAS_GENES
    )
    gene_rows = []
    for gene in genes:
        host_has = np.array([gene in _parse_genes(anno[s]["Specific_genes"]) for s in matched_hosts])
        nonhost_has = np.array([gene in _parse_genes(anno[s]["Specific_genes"]) for s in matched_nonhosts])
        table = [[int(host_has.sum()), int((~host_has).sum())],
                 [int(nonhost_has.sum()), int((~nonhost_has).sum())]]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        discordant_host_only = int((host_has & ~nonhost_has).sum())
        discordant_nonhost_only = int((~host_has & nonhost_has).sum())
        discordant_total = discordant_host_only + discordant_nonhost_only
        if discordant_total:
            mcnemar_p = float(2 * min(
                binomtest(discordant_host_only, discordant_total, 0.5, alternative="less").pvalue,
                binomtest(discordant_host_only, discordant_total, 0.5, alternative="greater").pvalue,
            ))
            mcnemar_p = min(mcnemar_p, 1.0)
        else:
            mcnemar_p = 1.0
        gene_rows.append({
            "gene": gene,
            "host_carriers": table[0][0],
            "host_noncarriers": table[0][1],
            "matched_nonhost_carriers": table[1][0],
            "matched_nonhost_noncarriers": table[1][1],
            "odds_ratio": float(odds_ratio),
            "fisher_p_unpaired_reference": float(p_value),
            "host_only_discordant": discordant_host_only,
            "nonhost_only_discordant": discordant_nonhost_only,
            "mcnemar_exact_p": mcnemar_p,
            "mcnemar_bonferroni_p": min(mcnemar_p * len(genes), 1.0),
        })

    pair_rows = []
    for h, n, d in zip(matched_hosts, matched_nonhosts, distances):
        pair_rows.append({
            "host_sgb": h,
            "nonhost_sgb": n,
            "host_log10_mean_abundance": anno[h]["Log10_abund"],
            "nonhost_log10_mean_abundance": anno[n]["Log10_abund"],
            "absolute_log10_distance": float(d),
            "host_narnas_count": anno[h]["N_narnas"],
            "nonhost_narnas_count": anno[n]["N_narnas"],
        })
    write_csv(OUT / "matched_pairs.csv", pair_rows)
    write_csv(OUT / "gene_tests.csv", gene_rows)

    summary = {
        "matching": "global minimum-cost one-to-one assignment without replacement",
        "matching_variable": "Log10_mean_abundance",
        "n_hosts_available": len(hosts),
        "n_nonhosts_available": len(non_hosts),
        "n_pairs": len(pair_rows),
        "n_unique_nonhosts_used": len(set(matched_nonhosts)),
        "mean_absolute_log10_distance": float(np.mean(distances)),
        "median_absolute_log10_distance": float(np.median(distances)),
        "max_absolute_log10_distance": float(np.max(distances)),
        "smd_before_matching": smd(host_values, nonhost_values),
        "smd_after_matching": smd(host_values[host_idx], nonhost_values[nonhost_idx]),
        "paired_wilcoxon_statistic": float(w_stat),
        "paired_wilcoxon_p": float(w_p),
        "two_sample_mann_whitney_statistic": float(mw_stat),
        "two_sample_mann_whitney_p": float(mw_p),
        "n_genes_tested": len(gene_rows),
        "n_mcnemar_bonferroni_significant": sum(r["mcnemar_bonferroni_p"] < 0.05 for r in gene_rows),
        "host_narnas_mean": float(np.mean(host_counts)),
        "matched_nonhost_narnas_mean": float(np.mean(nonhost_counts)),
        "mean_difference_host_minus_nonhost": float(np.mean(host_counts - nonhost_counts)),
        "replacement_control_reference": {
            "current_method": "nearest neighbour with replacement",
            "current_script": "08_core_statistics/narnas_stats.py:252-256",
            "warning": "This output is a sensitivity analysis, not a replacement of current results.",
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
