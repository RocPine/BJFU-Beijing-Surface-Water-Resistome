# -*- coding: utf-8 -*-
"""ED Fig 9b (原 Fig 6b): NAR/NAP/NAS 口径下 High 组丰度的 OLS 三模型 (m0/m1/m2)

数据来源: SI Tables (Supplementary_Tables.xlsx)
- Table 1:  站点坐标 (Sample-id, Longitude, Latitude)
- Table 7:  环境变量 (NO3_N, TOC, COD, TP)
- Table 10: SGB 相对丰度矩阵
- Table 15: SGB 注释 (ARG_host, NAR_NAP_NAS_total)
"""
import sys, yaml, pathlib
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import numpy as np
import pandas as pd
from scipy import stats
from pyproj import Transformer

EPS = 1e-6


def load_config():
    ROOT = pathlib.Path(__file__).resolve().parents[1]
    with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    # 相对路径统一解析为基于仓库根目录的绝对路径
    for k, v in cfg.items():
        if isinstance(v, str) and not pathlib.Path(v).is_absolute():
            cfg[k] = str(ROOT / v)
    return cfg


def dbmem(coords):
    """Distance-based Moran eigenvector maps"""
    n = len(coords)
    D = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    trunc = np.max(D)
    D_tr = D.copy()
    D_tr[D_tr > trunc] = 4 * trunc
    A = -0.5 * D_tr ** 2
    I_ = np.eye(n) - np.ones((n, n)) / n
    G = I_ @ A @ I_
    evals, evecs = np.linalg.eigh(G)
    pos = evals > 1e-8
    V = evecs[:, pos] * np.sqrt(np.abs(evals[pos]))
    order = np.argsort(evals[pos])[::-1]
    return V[:, order], evals[pos][order]


def ols(y_col, x_cols, data):
    n = len(data)
    X = np.column_stack([np.ones(n)] + [data[c].values for c in x_cols])
    yv = data[y_col].values
    beta = np.linalg.solve(X.T @ X, X.T @ yv)
    resid = yv - X @ beta
    dof = n - X.shape[1]
    s2 = resid @ resid / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = beta / se
    p_val = 2 * (1 - stats.t.cdf(np.abs(t), dof))
    return beta, se, p_val


def compute_models(cfg):
    """返回站点级数据框 df 与三个模型结果 [(beta_no3, se_no3, p_no3), ...]"""
    wb = openpyxl.load_workbook(cfg["si_tables_xlsx"], read_only=True)

    # Table 15: ARG_host + NAR_NAP_NAS_total
    ws15 = wb["Table 15"]
    it = ws15.iter_rows(values_only=True)
    next(it)
    hdr = [str(h) if h else "" for h in next(it)]
    i_host = hdr.index("ARG_host")
    i_ntt = hdr.index("NAR_NAP_NAS_total")
    sgb_anno = {}
    for row in ws15.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        sgb_anno[str(row[0]).strip()] = {
            "ARG_host": row[i_host],
            "N_narnas": row[i_ntt] if isinstance(row[i_ntt], (int, float)) else 0,
        }
    hosts = [s for s, a in sgb_anno.items() if a["ARG_host"] == "Yes"]
    high = [s for s in hosts if sgb_anno[s]["N_narnas"] >= 1]

    # Table 10: SGB 丰度
    ws11 = wb["Table 10"]
    it = ws11.iter_rows(values_only=True)
    next(it)
    hdr11 = [str(h) if h else "" for h in next(it)]
    sample_cols = [s for s in hdr11[1:] if s and s != 'None']
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
    env_by_sample = {}
    for row in ws7.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        env_by_sample[str(row[0]).strip()] = {h: row[i] for i, h in enumerate(hdr7) if h}

    # Table 1: 站点坐标
    ws1 = wb["Table 1"]
    it1 = ws1.iter_rows(values_only=True)
    next(it1)
    hdr1 = [str(h) if h else "" for h in next(it1)]
    i_lon = hdr1.index("Longitude")
    i_lat = hdr1.index("Latitude")
    site_coord = {}
    for row in ws1.iter_rows(min_row=3, values_only=True):
        if row[0] is None:
            continue
        try:
            site_coord[str(row[0]).strip()] = (float(row[i_lon]), float(row[i_lat]))
        except (TypeError, ValueError):
            continue
    wb.close()

    # 站点聚合 (与 narnas_stats 相同的 sample_to_site 逻辑; QSH -> QSH2 对齐 Table 1)
    import re
    sites = [re.sub(r'(_\d+|\d+)$', '', s.strip()) for s in sample_cols]
    sites = ["QSH2" if s == "QSH" else s for s in sites]
    site_samples = {}
    for s, site in zip(sample_cols, sites):
        site_samples.setdefault(site, []).append(s)

    def site_abund(group_sgbs):
        out = {}
        for site, ss in site_samples.items():
            sums = []
            for s in ss:
                if s not in sample_idx:
                    continue
                si = sample_idx[s]
                sums.append(sum(sgb_abund[sgb][si] for sgb in group_sgbs if sgb in sgb_abund))
            if sums:
                out[site] = np.mean(sums)
        return out

    sa_high = site_abund(high)

    def site_env_mean(site, key):
        ss = site_samples[site]
        vals = [env_by_sample[s][key] for s in ss if s in env_by_sample and env_by_sample[s].get(key) is not None]
        return np.mean([float(v) for v in vals]) if vals else np.nan

    def site_no3_mean(site):
        ss = site_samples[site]
        vals = [env_by_sample[s]["NO3_N"] for s in ss if s in env_by_sample and env_by_sample[s].get("NO3_N") is not None]
        return np.mean([float(v) for v in vals]) if vals else np.nan

    rows = []
    for site in sorted(site_samples):
        if site not in sa_high or site not in site_coord:
            continue
        no3 = site_no3_mean(site)
        toc = site_env_mean(site, "TOC")
        cod = site_env_mean(site, "COD")
        tp = site_env_mean(site, "TP")
        do = site_env_mean(site, "DO")
        if np.isnan(no3) or np.isnan(toc) or np.isnan(cod) or np.isnan(tp) or np.isnan(do):
            continue
        lon, lat = site_coord[site]
        tr = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
        e, n = tr.transform(lon, lat)
        rows.append({"site": site, "high": sa_high[site], "NO3": no3,
                     "TOC": toc, "COD": cod, "TP": tp, "DO": do, "east": e, "north": n})
    df = pd.DataFrame(rows)
    df["high_log"] = np.log10(df["high"] + EPS)
    df["NO3_log"] = np.log10(df["NO3"] + EPS)
    df["TOC_log"] = np.log10(df["TOC"] + EPS)
    df["COD_log"] = np.log10(df["COD"] + EPS)
    df["TP_log"] = np.log10(df["TP"] + EPS)

    b0, se0, p0 = ols("high_log", ["NO3_log"], df)
    b1, se1, p1 = ols("high_log", ["NO3_log", "TOC_log", "COD_log", "TP_log"], df)

    # MEM 前向选择
    coords2 = np.column_stack([df.east, df.north])
    V, _ = dbmem(coords2)
    y = df.high_log.values
    X_cur = np.ones((len(df), 1))
    selected = []
    remain = list(range(V.shape[1]))
    for _ in range(min(15, len(remain))):
        best_p, best_j = 1.0, None
        for j in remain:
            Xt = np.column_stack([X_cur, V[:, j]])
            try:
                beta = np.linalg.lstsq(Xt, y, rcond=None)[0]
                resid = y - Xt @ beta
                dof = len(y) - Xt.shape[1]
                se = np.sqrt(np.diag(resid @ resid / dof * np.linalg.inv(Xt.T @ Xt)))[-1]
                t = beta[-1] / se
                p = 2 * (1 - stats.t.cdf(abs(t), dof))
            except Exception:
                continue
            if p < best_p:
                best_p, best_j = p, j
        if best_p < 0.05 and best_j is not None:
            selected.append(best_j)
            X_cur = np.column_stack([X_cur, V[:, best_j]])
            remain.remove(best_j)
        else:
            break

    df_mem = pd.DataFrame({"high_log": df.high_log.values, "NO3_log": df.NO3_log.values})
    for k, j in enumerate(selected):
        df_mem[f"mem{k}"] = V[:, j]
    cols_mem = [f"mem{k}" for k in range(len(selected))]
    b2, se2, p2 = ols("high_log", ["NO3_log"] + cols_mem, df_mem)

    # DO 分层: m0 型 NO3 系数 (ED Fig 9b 后半段)
    med_do = df["DO"].median()
    dh = df[df["DO"] >= med_do]
    dl = df[df["DO"] < med_do]
    bh, seh, ph = ols("high_log", ["NO3_log"], dh)
    bl, sel, pl = ols("high_log", ["NO3_log"], dl)

    return {
        "df": df, "n_high": len(high), "n_sites": len(df), "n_mem_selected": len(selected),
        "models": [("m0", b0[1], se0[1], p0[1]), ("m1", b1[1], se1[1], p1[1]), ("m2", b2[1], se2[1], p2[1])],
        "do_models": [("high-DO", bh[1], seh[1], ph[1], len(dh)), ("low-DO", bl[1], sel[1], pl[1], len(dl))],
    }


def main():
    cfg = load_config()
    r = compute_models(cfg)
    print(f"High (NAR/NAP/NAS >=1): {r['n_high']} SGBs; {r['n_sites']} sites; MEM selected: {r['n_mem_selected']}")
    for name, b, se, p in r["models"]:
        print(f"  {name}: beta={b:.4f} se={se:.4f} P={p:.4g}")


if __name__ == "__main__":
    main()
