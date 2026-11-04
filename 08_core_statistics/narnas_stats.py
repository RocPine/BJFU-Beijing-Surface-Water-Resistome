# -*- coding: utf-8 -*-
"""核心统计: NAR/NAP/NAS 口径 (Fig 6a 梯度 / 6b ECDF / 6c 匹配 / 6d 火山及 ED Fig 9/10 共用)

数据来源: SI Tables (Supplementary_Tables.xlsx)
- Table 7:  环境变量 (NO3_N)
- Table 9:  SGB 分类 (MAG -> SGB 映射)
- Table 10: SGB 相对丰度矩阵
- Table 15: 氮循环基因携带 (ARG_host, NAR_NAP_NAS_total, Specific_genes,
             Log10_mean_abundance)
"""
import sys, re, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import numpy as np
from scipy.stats import binomtest, fisher_exact, mannwhitneyu, spearmanr, wilcoxon
from scipy.optimize import linear_sum_assignment


def load_config():
    ROOT = pathlib.Path(__file__).resolve().parents[1]
    with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    # 相对路径统一解析为基于仓库根目录的绝对路径
    for k, v in cfg.items():
        if isinstance(v, str) and not pathlib.Path(v).is_absolute():
            cfg[k] = str(ROOT / v)
    return cfg


def load_sgb_annotation(cfg):
    """Table 15: SGB 注释 (ARG_host, NAR_NAP_NAS_total, Specific_genes, Log10_mean_abundance)"""
    wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)
    ws = wb["Table 15"]
    it = ws.iter_rows(values_only=True)
    next(it)
    hdr = [str(h) if h else "" for h in next(it)]
    idx = {h: i for i, h in enumerate(hdr)}
    anno = {}
    for row in ws.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        sgb = str(row[0]).strip()
        anno[sgb] = {
            "ARG_host": row[idx["ARG_host"]],
            "N_narnas": row[idx["NAR_NAP_NAS_total"]] if isinstance(row[idx["NAR_NAP_NAS_total"]], (int, float)) else 0,
            "Specific_genes": str(row[idx["Specific_genes"]]) if row[idx["Specific_genes"]] else "",
            "Log10_abund": row[idx["Log10_mean_abundance"]] if isinstance(row[idx["Log10_mean_abundance"]], (int, float)) else -9,
        }
    wb.close()
    return anno


def sample_to_site(s):
    """站点名推导: 样本名去掉尾部复制编号 (ASB1->ASB, BAIH1_1->BAIH1)"""
    return re.sub(r'(_\d+|\d+)$', '', str(s).strip())


def load_site_data(cfg):
    """Table 10 (SGB 丰度) + Table 7 (环境): 站点级 NO3-N 均值与三组 SGB 丰度"""
    wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)

    # Table 10: SGB 丰度矩阵
    ws11 = wb["Table 10"]
    it = ws11.iter_rows(values_only=True)
    next(it)
    header = [str(h) if h else "" for h in next(it)]
    sample_cols = [s for s in header[1:] if s and s != 'None']
    sites = [sample_to_site(s) for s in sample_cols]
    sample_idx = {s: i for i, s in enumerate(sample_cols)}
    sgb_abund = {}
    for row in ws11.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        sgb = str(row[0]).strip()
        vals = [float(x) if isinstance(x, (int, float)) else 0.0 for x in row[1:1 + len(sample_cols)]]
        sgb_abund[sgb] = vals

    # Table 7: 环境变量
    ws7 = wb["Table 7"]
    it7 = ws7.iter_rows(values_only=True)
    next(it7)
    hdr7 = [str(h) if h else "" for h in next(it7)]
    idx_no3 = hdr7.index("NO3_N")
    env_by_sample = {}
    for row in ws7.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        env_by_sample[str(row[0]).strip()] = row[idx_no3]
    wb.close()

    site_samples = {}
    for s, site in zip(sample_cols, sites):
        site_samples.setdefault(site, []).append(s)

    site_no3 = {}
    for site, ss in site_samples.items():
        vals = [env_by_sample[s] for s in ss if s in env_by_sample and env_by_sample[s] is not None]
        if vals:
            site_no3[site] = np.mean([float(v) for v in vals])

    return {
        "sample_idx": sample_idx,
        "sgb_abund": sgb_abund,
        "site_samples": site_samples,
        "site_no3": site_no3,
    }


def load_site_env(cfg, keys):
    """Table 7 站点级环境变量均值 {key: {site: mean}} (供 ED Fig 9c/9d 使用)"""
    wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)
    ws7 = wb["Table 7"]
    it7 = ws7.iter_rows(values_only=True)
    next(it7)
    hdr7 = [str(h) if h else "" for h in next(it7)]
    idx = {h: i for i, h in enumerate(hdr7)}
    per_sample = {k: {} for k in keys}
    for row in ws7.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        site = sample_to_site(str(row[0]).strip())
        for k in keys:
            v = row[idx[k]]
            if isinstance(v, (int, float)):
                per_sample[k].setdefault(site, []).append(v)
    wb.close()
    return {k: {s: np.mean(v) for s, v in d.items()} for k, d in per_sample.items()}


def site_abund(site_data, group_sgbs):
    """站点级相对丰度: 组内 SGB 丰度求和, 复制样本取均值"""
    out = {}
    for site, ss in site_data["site_samples"].items():
        sums = []
        for s in ss:
            if s not in site_data["sample_idx"]:
                continue
            si = site_data["sample_idx"][s]
            total = sum(site_data["sgb_abund"][sgb][si] for sgb in group_sgbs if sgb in site_data["sgb_abund"])
            sums.append(total)
        if sums:
            out[site] = np.mean(sums)
    return out


def site_spearman(site_no3, sa):
    common = sorted(set(site_no3) & set(sa))
    x = np.array([site_no3[s] for s in common])
    y = np.array([sa[s] for s in common])
    rho, p = spearmanr(x, y)
    return rho, p, x, y


def _split_hosts(anno):
    hosts = [s for s, a in anno.items() if a["ARG_host"] == "Yes"]
    non_hosts = [s for s, a in anno.items() if a["ARG_host"] != "Yes"]
    return hosts, non_hosts


def compute_fig5e(cfg):
    """Fig 6b (原 Fig 5e) 携带率: Fisher exact + 基因数 Mann-Whitney
    (函数名保留旧编号, 被 06b_narnas_ecdf.py 调用)"""
    anno = load_sgb_annotation(cfg)
    hosts, non_hosts = _split_hosts(anno)
    h_counts = np.array([anno[s]["N_narnas"] for s in hosts])
    n_counts = np.array([anno[s]["N_narnas"] for s in non_hosts])
    h_car = int((h_counts >= 1).sum())
    n_car = int((n_counts >= 1).sum())
    table = [[h_car, len(hosts) - h_car], [n_car, len(non_hosts) - n_car]]
    # 样本 OR (cross-product), P 用 Fisher exact (单侧, 与原文口径一致)
    or_sample = table[0][0] * table[1][1] / (table[0][1] * table[1][0])
    _, p_fisher = fisher_exact(table, alternative="greater")
    _, p_mw = mannwhitneyu(h_counts, n_counts, alternative="greater")
    return {
        "n_hosts": len(hosts), "n_nonhosts": len(non_hosts), "n_high": h_car,
        "carriage_host_pct": h_car / len(hosts) * 100,
        "carriage_nonhost_pct": n_car / len(non_hosts) * 100,
        "or": or_sample, "p_fisher": p_fisher,
        "host_mean": h_counts.mean(), "nonhost_mean": n_counts.mean(),
        "p_mw": p_mw, "host_counts": h_counts, "nonhost_counts": n_counts,
        "hosts": hosts, "non_hosts": non_hosts,
    }


def compute_fig6a(cfg, anno=None):
    """Fig 6a 三组梯度: Spearman + bootstrap
    (函数名保留旧编号, 被 06a_three_group_gradient.py / 07i_edfig9a_bootstrap.py 调用)"""
    if anno is None:
        anno = load_sgb_annotation(cfg)
    hosts, non_hosts = _split_hosts(anno)
    high = [s for s in hosts if anno[s]["N_narnas"] >= 1]
    low = [s for s in hosts if anno[s]["N_narnas"] < 1]

    site_data = load_site_data(cfg)
    sa_high = site_abund(site_data, high)
    sa_low = site_abund(site_data, low)
    sa_nonhost = site_abund(site_data, non_hosts)

    rho_h, p_h, x_h, y_h = site_spearman(site_data["site_no3"], sa_high)
    rho_l, p_l, x_l, y_l = site_spearman(site_data["site_no3"], sa_low)
    rho_n, p_n, x_n, y_n = site_spearman(site_data["site_no3"], sa_nonhost)

    # Bootstrap
    rng = np.random.default_rng(42)
    common_sites = sorted(set(site_data["site_no3"]) & set(sa_high) & set(sa_low) & set(sa_nonhost))
    xn = np.array([site_data["site_no3"][s] for s in common_sites])
    yh = np.array([sa_high[s] for s in common_sites])
    yl = np.array([sa_low[s] for s in common_sites])
    yn = np.array([sa_nonhost[s] for s in common_sites])
    n_hl = 0
    n_hn = 0
    for _ in range(1000):
        idx = rng.integers(0, len(common_sites), len(common_sites))
        r_h = spearmanr(xn[idx], yh[idx])[0]
        r_l = spearmanr(xn[idx], yl[idx])[0]
        r_n = spearmanr(xn[idx], yn[idx])[0]
        if r_h > r_l:
            n_hl += 1
        if r_h > r_n:
            n_hn += 1

    return {
        "high": high, "low": low, "non_hosts": non_hosts,
        "n_high": len(high), "n_low": len(low),
        "rho_h": rho_h, "p_h": p_h, "x_h": x_h, "y_h": y_h,
        "rho_l": rho_l, "p_l": p_l, "x_l": x_l, "y_l": y_l,
        "rho_n": rho_n, "p_n": p_n, "x_n": x_n, "y_n": y_n,
        "boot_high_gt_low": n_hl / 1000, "boot_high_gt_nonhost": n_hn / 1000,
        "site_no3": site_data["site_no3"],
    }


def _parse_genes(specific_genes):
    """Table 15 Specific_genes 列 -> 基因名集合"""
    return set(g.strip().lower() for g in str(specific_genes).split(";") if g.strip())


# Fig 6d 只统计 NAR/NAP/NAS 家族基因 (与火山图 GENE_GTYPE 白名单一致)
NARNAS_GENES = {
    "nar", "narb", "narc", "narg", "narh", "nari", "narj", "nary", "narz",
    "napa", "napb", "napc", "ansb", "gdh", "glsa", "gs", "nmo",
    "nasa", "nasb",
}


def compute_matching(cfg, anno=None):
    """Fig 6c/6d: global one-to-one abundance matching without replacement."""
    if anno is None:
        anno = load_sgb_annotation(cfg)
    hosts, non_hosts = _split_hosts(anno)

    h_ab = np.array([anno[s]["Log10_abund"] for s in hosts], dtype=float)
    nh_ab = np.array([anno[s]["Log10_abund"] for s in non_hosts], dtype=float)
    if len(hosts) > len(non_hosts):
        raise ValueError("No-replacement matching requires at least as many non-hosts as hosts")
    cost = np.abs(h_ab[:, None] - nh_ab[None, :])
    host_idx, nonhost_idx = linear_sum_assignment(cost)
    matched_nh = [non_hosts[i] for i in nonhost_idx]

    h_cnt = np.array([anno[s]["N_narnas"] for s in hosts])
    m_cnt = np.array([anno[s]["N_narnas"] for s in matched_nh])
    _, w_p = wilcoxon(h_cnt, m_cnt)
    mean_diff = h_cnt.mean() - m_cnt.mean()

    # Paired exact McNemar tests for each detected NAR/NAP/NAS gene.
    gene_set_all = sorted(
        set().union(*[_parse_genes(anno[s]["Specific_genes"]) for s in hosts]) & NARNAS_GENES)
    gene_results = []
    for gene in gene_set_all:
        h_has_vec = np.array([gene in _parse_genes(anno[s]["Specific_genes"]) for s in hosts])
        m_has_vec = np.array([gene in _parse_genes(anno[s]["Specific_genes"]) for s in matched_nh])
        h_has = int(h_has_vec.sum())
        m_has = int(m_has_vec.sum())
        host_only = int((h_has_vec & ~m_has_vec).sum())
        nonhost_only = int((~h_has_vec & m_has_vec).sum())
        discordant = host_only + nonhost_only
        if discordant:
            pv = min(2 * min(
                binomtest(host_only, discordant, 0.5, alternative="less").pvalue,
                binomtest(host_only, discordant, 0.5, alternative="greater").pvalue,
            ), 1.0)
        else:
            pv = 1.0
        orv = (host_only + 0.5) / (nonhost_only + 0.5)
        gene_results.append({"gene": gene, "or": orv, "p": pv,
                             "host_only": host_only, "nonhost_only": nonhost_only,
                             "host_carriers": h_has, "matched_carriers": m_has})
    bonf = 0.05 / len(gene_results)
    n_sig = sum(r["p"] < bonf for r in gene_results)

    return {
        "matched_nh": matched_nh, "host_counts": h_cnt, "matched_counts": m_cnt,
        "mean_diff": mean_diff, "p_wilcoxon": w_p,
        "or_results": [(r["gene"], r["or"], r["p"]) for r in gene_results],
        "gene_results": gene_results, "bonf": bonf, "n_sig": n_sig,
        "hosts": hosts, "matched_nh": matched_nh,
    }


def main():
    cfg = load_config()
    anno = load_sgb_annotation(cfg)
    r5 = compute_fig5e(cfg)
    print(f"Fig 6b (ECDF): hosts {r5['carriage_host_pct']:.1f}% vs non-hosts {r5['carriage_nonhost_pct']:.1f}% "
          f"(OR={r5['or']:.2f}, P={r5['p_fisher']:.2e})")
    print(f"Fig 6b (ECDF): hosts mean={r5['host_mean']:.2f}, non-hosts mean={r5['nonhost_mean']:.2f}, MW P={r5['p_mw']:.2e}")
    r6 = compute_fig6a(cfg, anno)
    print(f"Fig 6a: High(n={r6['n_high']}) rho={r6['rho_h']:.3f} P={r6['p_h']:.2e} | "
          f"Low(n={r6['n_low']}) rho={r6['rho_l']:.3f} P={r6['p_l']:.2e} | "
          f"Nonhost rho={r6['rho_n']:.3f} P={r6['p_n']:.2e}")
    print(f"Fig 6a bootstrap (ED Fig 9a): High>Low={r6['boot_high_gt_low']:.3f}, High>Non-host={r6['boot_high_gt_nonhost']:.3f}")
    rm = compute_matching(cfg, anno)
    print(f"Fig 6c: host mean={rm['host_counts'].mean():.2f}, matched={rm['matched_counts'].mean():.2f}, "
          f"diff=+{rm['mean_diff']:.2f}, paired Wilcoxon P={rm['p_wilcoxon']:.2e}")
    print(f"Fig 6d: {len(rm['or_results'])} genes, {rm['n_sig']} significant after Bonferroni (P < {rm['bonf']:.4f})")
    for g, o, p in rm["or_results"]:
        print(f"  {g}: OR={o:.2f}, P={p:.2e}")


if __name__ == "__main__":
    main()
