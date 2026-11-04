# -*- coding: utf-8 -*-
"""Fig 4 共享数据加载模块 (仅依赖随论文发布的 Supplementary Tables.xlsx)

数据来源:
- Table 12: ARG-associated HGT events (Sample, donor_mag, recipient_mag,
  ARG_Subtype, ARG_Type, donor_p, recipient_p)
- Table 9:  SGB classification (SGB, Member_MAG, MAG_taxonomy) -> MAG->SGB 映射
- Table 10: SGB relative abundance matrix (SGB x sample)
- Table 15: SGB annotation (ARG_host, N_total_NCGs, Log10_mean_abundance)
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import numpy as np
import pandas as pd


def _read_rows(ws, skip=1):
    """读取工作表: 跳过标题行, 返回 [(表头), [数据行...]]"""
    it = ws.iter_rows(values_only=True)
    for _ in range(skip):
        hdr = next(it)
    rows = []
    for r in it:
        if r[0] is None:
            continue
        rows.append(r)
    return [str(h) if h else "" for h in hdr], rows


def load_fig4_data(cfg):
    """
    仅从发布版 SI Tables 加载 Fig 4 所需数据。
    cfg: config.yaml 解析结果 (需含 si_tables_xlsx)
    """
    wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)

    # ---------- Table 12: HGT 事件 ----------
    ws = wb["Table 12"]
    it = ws.iter_rows(values_only=True)
    next(it)
    hdr = [str(h) if h else "" for h in next(it)]
    hgt_rows = []
    for r in it:
        if r[0] is None:
            continue
        hgt_rows.append(dict(zip(hdr, r)))
    hgt = pd.DataFrame(hgt_rows)
    hgt["sample"] = hgt["Sample"].astype(str)
    # donor/recipient 门 (无 p__ 前缀则补上, 保证与分类学命名一致)
    hgt["donor_phy"] = hgt["donor_p"].fillna("unknown").astype(str)
    hgt["recipient_phy"] = hgt["recipient_p"].fillna("unknown").astype(str)

    # ---------- Table 9: MAG -> SGB 映射 ----------
    ws9 = wb["Table 9"]
    it9 = ws9.iter_rows(values_only=True)
    next(it9)
    hdr9 = [str(h) if h else "" for h in next(it9)]
    i_member = hdr9.index("Member_MAG")
    mag2sgb = {}
    for r in it9:
        if r[0] is None:
            continue
        sgb = str(r[0]).strip()
        member_mags = r[i_member]
        if member_mags:
            for mag in str(member_mags).split(";"):
                mag = mag.strip()
                if mag:
                    mag2sgb[mag] = sgb

    # ---------- Table 10: SGB 丰度矩阵 ----------
    ws11 = wb["Table 10"]
    it11 = ws11.iter_rows(values_only=True)
    next(it11)
    hdr11 = [str(h) if h else "" for h in next(it11)]
    samples = [s for s in hdr11[1:] if s and s != "None"]
    sample_idx = {s: i for i, s in enumerate(samples)}
    sgb_abund = {}
    for r in it11:
        if r[0] is None:
            continue
        sgb = str(r[0]).strip()
        vals = [float(x) if isinstance(x, (int, float)) else 0.0 for x in r[1:1 + len(samples)]]
        sgb_abund[sgb] = vals
    sgb_set = set(sgb_abund)

    # ---------- Table 15: SGB 注释 (氮代谢基因携带表) ----------
    ws15 = wb["Table 15"]
    it15 = ws15.iter_rows(values_only=True)
    next(it15)
    hdr15 = [str(h) if h else "" for h in next(it15)]
    i_host = hdr15.index("ARG_host")
    i_nt = hdr15.index("N_total_NCGs")
    i_log = hdr15.index("Log10_mean_abundance")
    sgb_anno = {}
    for r in it15:
        if r[0] is None:
            continue
        sgb = str(r[0]).strip()
        sgb_anno[sgb] = {
            "ARG_host": r[i_host],
            "N_total": r[i_nt] if isinstance(r[i_nt], (int, float)) else 0,
            "Log10_abund": r[i_log] if isinstance(r[i_log], (int, float)) else -9,
        }
    wb.close()

    # ---------- 解析 MAG -> SGB ----------
    # 规则: MAG 名若本身就是 SGB 名则直接使用, 否则经 Member_MAG 映射
    sgb_names_set = set(sgb_abund)

    def resolve(mag):
        mag = str(mag).strip()
        if mag in sgb_names_set:
            return mag
        if mag in mag2sgb:
            return mag2sgb[mag]
        return None

    hgt["donor_sgb"] = hgt["donor_mag"].map(resolve)
    hgt["recipient_sgb"] = hgt["recipient_mag"].map(resolve)

    # ---------- 样本级门丰度与 HGT 事件汇总 (Fig 4d) ----------
    all_phys = set(hgt["donor_phy"]) | set(hgt["recipient_phy"])
    records = []
    for s in samples:
        si = sample_idx[s]
        for phy in all_phys:
            d_ev = len(hgt[(hgt["sample"] == s) & (hgt["donor_phy"] == phy)])
            r_ev = len(hgt[(hgt["sample"] == s) & (hgt["recipient_phy"] == phy)])
            if d_ev == 0 and r_ev == 0:
                continue
            d_sgbs = hgt[(hgt["sample"] == s) & (hgt["donor_phy"] == phy)]["donor_sgb"].dropna().unique()
            r_sgbs = hgt[(hgt["sample"] == s) & (hgt["recipient_phy"] == phy)]["recipient_sgb"].dropna().unique()
            d_abund = sum(sgb_abund[g][si] for g in d_sgbs if g in sgb_set)
            r_abund = sum(sgb_abund[g][si] for g in r_sgbs if g in sgb_set)
            records.append({"sample": s, "phylum": phy,
                            "donor_abund": d_abund, "recipient_abund": r_abund,
                            "donor_events": d_ev, "recipient_events": r_ev})
    df = pd.DataFrame(records)

    # ---------- 频率: 事件计数 (Table 12 每行一事件; SI 无标准化丰度) ----------
    subtype_class = dict(zip(hgt["ARG_Subtype"], hgt["ARG_Type"]))
    subtype_class["vanH-A"] = "glycopeptide antibiotic"
    subtype_class["vanH-O"] = "glycopeptide antibiotic"

    def short_name(s):
        return str(s).replace("vanH gene in ", "").replace(" antibiotic", "").replace("_tiamulin", "")

    # ---------- 门级事件汇总 ----------
    donor_phys = hgt.groupby("donor_phy").size()
    recv_phys = hgt.groupby("recipient_phy").size()
    all_phy = set(donor_phys.index) | set(recv_phys.index)
    total_by_phy = {p: donor_phys.get(p, 0) + recv_phys.get(p, 0) for p in all_phy}
    node_order = sorted(all_phy, key=lambda p: -total_by_phy.get(p, 0))
    cross_pp = pd.crosstab(hgt["donor_phy"], hgt["recipient_phy"])

    # ---------- Panel e 数据: donor/recipient/other hosts 的 N_total ----------
    # 与原文一致: 事件行级计数 (每个 HGT 事件行计一次), 且仅统计 ARG-host 的 donor/recipient
    donor_n_e = np.array([
        sgb_anno[d]["N_total"]
        for d in hgt.loc[hgt["donor_sgb"].notna(), "donor_sgb"]
        if d in sgb_anno and sgb_anno[d]["ARG_host"] == "Yes"
    ], float)
    recipient_n_e = np.array([
        sgb_anno[r]["N_total"]
        for r in hgt.loc[hgt["recipient_sgb"].notna(), "recipient_sgb"]
        if r in sgb_anno and sgb_anno[r]["ARG_host"] == "Yes"
    ], float)
    donor_sgbs = {s for s in hgt["donor_sgb"].dropna() if s in sgb_anno and sgb_anno[s]["ARG_host"] == "Yes"}
    recipient_sgbs = {s for s in hgt["recipient_sgb"].dropna() if s in sgb_anno and sgb_anno[s]["ARG_host"] == "Yes"}
    hgt_sgbs_e = donor_sgbs | recipient_sgbs
    non_donor_n_e = np.array(
        [a["N_total"] for s, a in sgb_anno.items()
         if a["ARG_host"] == "Yes" and s not in hgt_sgbs_e], float)

    return {
        "hgt": hgt,
        "samples": samples,
        "sample_idx": sample_idx,
        "sgb_set": sgb_set,
        "df": df,
        "subtype_class": subtype_class,
        "short_name": short_name,
        "donor_phys": donor_phys,
        "recv_phys": recv_phys,
        "all_phy": all_phy,
        "total_by_phy": total_by_phy,
        "node_order": node_order,
        "cross_pp": cross_pp,
        "donor_n_e": donor_n_e,
        "recipient_n_e": recipient_n_e,
        "non_donor_n_e": non_donor_n_e,
        "donor_sgbs": donor_sgbs,
        "recipient_sgbs": recipient_sgbs,
        "hgt_sgbs": hgt_sgbs_e,
        "sgb_anno": sgb_anno,
        "CLASS_COLORS": {
            "macrolide-lincosamide-streptogramin": "#E64B35", "aminoglycoside antibiotic": "#4DBBD5",
            "glycopeptide antibiotic": "#00A087", "multidrug": "#3C5488",
            "polymyxin": "#F39B7F", "mupirocin": "#8491B4", "pleuromutilin_tiamulin": "#7E6148",
        },
        "TYPE_SHORT": {
            "macrolide-lincosamide-streptogramin": "MLS",
            "aminoglycoside antibiotic": "Aminoglycoside",
            "glycopeptide antibiotic": "Glycopeptide",
            "multidrug": "Multidrug",
            "polymyxin": "Polymyxin",
            "mupirocin": "Mupirocin",
            "pleuromutilin_tiamulin": "Pleuromutilin",
        },
        "PHYLA_COLORS": {
            "p__Actinomycetota": "#D62728", "p__Bacteroidota": "#1F77B4",
            "p__Pseudomonadota": "#2CA02C", "p__Verrucomicrobiota": "#9467BD",
            "p__Planctomycetota": "#FF7F0E", "p__Acidobacteriota": "#17BECF",
            "p__Gemmatimonadota": "#8C564B", "p__Cyanobacteriota": "#E377C2",
            "p__Patescibacteria": "#BCBD22",
        }
    }