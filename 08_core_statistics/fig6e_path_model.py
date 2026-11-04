# -*- coding: utf-8 -*-
"""Exploratory site-level path model for Fig. 6e.

Uses original SI workbook. Site is independent unit (three subsamples/site).
This script reports associational paths and host-adjustment attenuation; it
does not estimate or label causal mediation.
"""
from __future__ import annotations

import pathlib
import re
import sys

import numpy as np
import pandas as pd
import yaml
from openpyxl import load_workbook
from scipy import stats
from scipy.stats import chi2
from pyproj import Transformer

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "fig6e_path_model"
EPS = 1e-6
N_BOOT = 5000
N_MORAN_PERM = 9999
SEED = 20260906


def read_table(wb, name):
    ws = wb[name]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(x) if x is not None else "" for x in rows[1]]
    return pd.DataFrame([r for r in rows[2:] if r and r[0] is not None], columns=header)


def site_name(sample):
    value = re.sub(r"(_\d+|\d+)$", "", str(sample).strip())
    return "QSH2" if value == "QSH" else value


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def load_site_frame():
    with open(ROOT / "config" / "config.yaml", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)
    si_path = ROOT / cfg["si_tables_xlsx"]
    wb = load_workbook(si_path, read_only=True, data_only=True)
    t4, t7, t10, t15, t1 = (read_table(wb, name) for name in ("Table 4", "Table 7", "Table 10", "Table 15", "Table 1"))
    wb.close()

    # Table 7 supplies sample-level environmental and human-activity values.
    t7["site"] = t7["Sample"].map(site_name)
    env_cols = ["NO3_N", "TOC", "COD", "TP", "HSCI"]
    env = t7[["site"] + env_cols].copy()
    for col in env_cols:
        env[col] = numeric(env[col])
    env = env.groupby("site", as_index=False)[env_cols].mean()

    # Table 10 supplies SGB abundance; sum nitrate-reducing ARG-host SGBs.
    anno = t15.set_index(t15.columns[0])
    host_col = "ARG_host"
    narnas_col = "NAR_NAP_NAS_total"
    high_sgbs = set(anno.index[(anno[host_col] == "Yes") & (numeric(anno[narnas_col]).fillna(0) >= 1)])
    samples = [str(x) for x in t10.columns[1:] if x not in ("", "None")]
    sgb = t10.set_index(t10.columns[0])
    abundance = sgb.apply(pd.to_numeric, errors="coerce").fillna(0)
    high_sgbs &= set(abundance.index)
    all_arg_hosts = set(anno.index[anno[host_col] == "Yes"]) & set(abundance.index)
    low_sgbs = all_arg_hosts - high_sgbs
    high_by_sample = abundance.loc[sorted(high_sgbs)].sum(axis=0)
    low_by_sample = abundance.loc[sorted(low_sgbs)].sum(axis=0)
    host = pd.DataFrame({
        "sample": samples,
        "host_abundance": [high_by_sample.get(s, 0.0) for s in samples],
        "noncarrying_host_abundance": [low_by_sample.get(s, 0.0) for s in samples],
    })
    host["site"] = host["sample"].map(site_name)
    host = host.groupby("site", as_index=False)[["host_abundance", "noncarrying_host_abundance"]].mean()

    # Table 4 is sample x ARG subtype; sum each sample row directly.
    arg_values = t4.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    arg = pd.DataFrame({"sample": t4.iloc[:, 0].astype(str), "arg_abundance": arg_values.to_numpy()})
    arg["site"] = arg["sample"].map(site_name)
    arg = arg.groupby("site", as_index=False)["arg_abundance"].mean()

    # Table 1 coordinates enable a prespecified spatial residual diagnostic.
    coords = t1[["Sample-id", "Longitude", "Latitude"]].copy()
    coords = coords.rename(columns={"Sample-id": "site"})
    coords["Longitude"] = numeric(coords["Longitude"])
    coords["Latitude"] = numeric(coords["Latitude"])
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:32650", always_xy=True)
    coords["east"], coords["north"] = transformer.transform(coords["Longitude"].to_numpy(), coords["Latitude"].to_numpy())
    frame = env.merge(host, on="site", how="inner").merge(arg, on="site", how="inner").merge(coords, on="site", how="inner")
    frame = frame.dropna(subset=env_cols + ["host_abundance", "noncarrying_host_abundance", "arg_abundance"]).copy()
    if len(frame) < 50:
        raise ValueError(f"Only {len(frame)} complete sites; expected approximately 89")
    for col in ["NO3_N", "host_abundance", "noncarrying_host_abundance", "arg_abundance"]:
        frame[f"log_{col}"] = np.log10(frame[col].clip(lower=0) + EPS)
    frame["host_composition_log_ratio"] = frame["log_host_abundance"] - frame["log_noncarrying_host_abundance"]
    return frame.sort_values("site").reset_index(drop=True), si_path, len(high_sgbs)


def fit_standardized(data, outcome, predictors):
    cols = [outcome] + predictors
    z = data[cols].apply(lambda s: (s - s.mean()) / s.std(ddof=1))
    y = z[outcome].to_numpy()
    x = np.column_stack([np.ones(len(z)), z[predictors].to_numpy()])
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    residual = y - x @ beta
    dof = len(y) - x.shape[1]
    mse = (residual @ residual) / dof
    cov = mse * np.linalg.pinv(x.T @ x)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    r2 = 1 - (residual @ residual) / ((y - y.mean()) @ (y - y.mean()))
    leverage = np.sum(x * (x @ np.linalg.pinv(x.T @ x)), axis=1)
    hc3 = (residual / np.maximum(1 - leverage, 1e-8)) ** 2
    robust_cov = np.linalg.pinv(x.T @ x) @ (x.T @ np.diag(hc3) @ x) @ np.linalg.pinv(x.T @ x)
    robust_se = np.sqrt(np.diag(robust_cov))
    robust_t = beta / robust_se
    robust_p = 2 * stats.t.sf(np.abs(robust_t), dof)
    return {"beta": beta[1:], "se": se[1:], "p": p[1:], "robust_se": robust_se[1:], "robust_p": robust_p[1:], "r2": r2, "residual": residual}


def bootstrap_paths(data, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(data)
    formulas = {
        "host_model": ("log_host_abundance", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
        "arg_total_model": ("log_arg_abundance", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
        "arg_host_adjusted_model": ("log_arg_abundance", ["log_NO3_N", "log_host_abundance", "HSCI", "TOC", "COD", "TP"]),
        "host_composition_model": ("host_composition_log_ratio", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
    }
    draws = {name: [] for name in formulas}
    for _ in range(N_BOOT):
        sample = data.iloc[rng.integers(0, n, n)].reset_index(drop=True)
        fitted = {name: fit_standardized(sample, outcome, predictors) for name, (outcome, predictors) in formulas.items()}
        no3_total = fitted["arg_total_model"]["beta"][0]
        no3_adjusted = fitted["arg_host_adjusted_model"]["beta"][0]
        host_to_arg = fitted["arg_host_adjusted_model"]["beta"][1]
        no3_to_host = fitted["host_model"]["beta"][0]
        draws["host_model"].append(no3_to_host)
        draws["arg_total_model"].append(no3_total)
        draws["arg_host_adjusted_model"].append(no3_adjusted)
        draws["host_composition_model"].append(fitted["host_composition_model"]["beta"][0])
        draws.setdefault("indirect_association", []).append(no3_to_host * host_to_arg)
        draws.setdefault("attenuation", []).append((no3_total - no3_adjusted) / no3_total if abs(no3_total) > EPS else np.nan)
    return {key: np.asarray(values, dtype=float) for key, values in draws.items()}


def moran_weights(coords):
    distances = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distances, np.inf)
    scale = np.nanmedian(distances[np.isfinite(distances)])
    weights = np.exp(-distances / scale)
    np.fill_diagonal(weights, 0)
    return weights


def moran_i(values, weights):
    centered = values - values.mean()
    return len(values) / weights.sum() * ((weights * np.outer(centered, centered)).sum() / (centered @ centered))


def moran_permutation(values, weights, n_perm=N_MORAN_PERM, seed=SEED):
    observed = moran_i(values, weights)
    rng = np.random.default_rng(seed)
    permuted = np.empty(n_perm)
    for i in range(n_perm):
        permuted[i] = moran_i(rng.permutation(values), weights)
    p_two_sided = (1 + np.count_nonzero(np.abs(permuted) >= abs(observed))) / (n_perm + 1)
    return observed, p_two_sided


def calculate_vif(data, predictors):
    values = data[predictors].apply(lambda s: (s - s.mean()) / s.std(ddof=1)).to_numpy()
    out = []
    for i, predictor in enumerate(predictors):
        target = values[:, i]
        other = np.delete(values, i, axis=1)
        x = np.column_stack([np.ones(len(values)), other])
        fitted = x @ np.linalg.lstsq(x, target, rcond=None)[0]
        r2 = 1 - np.sum((target - fitted) ** 2) / np.sum((target - target.mean()) ** 2)
        out.append({"predictor": predictor, "vif": 1 / max(1 - r2, 1e-12)})
    return out


def main():
    data, si_path, high_sgb_count = load_site_frame()
    formulas = {
        "host_model": ("log_host_abundance", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
        "arg_total_model": ("log_arg_abundance", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
        "arg_host_adjusted_model": ("log_arg_abundance", ["log_NO3_N", "log_host_abundance", "HSCI", "TOC", "COD", "TP"]),
        "host_composition_model": ("host_composition_log_ratio", ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]),
    }
    results = []
    fitted = {}
    diagnostics = []
    weights = moran_weights(data[["east", "north"]].to_numpy())
    for model, (outcome, predictors) in formulas.items():
        fit = fit_standardized(data, outcome, predictors)
        fitted[model] = fit
        for predictor, beta, se, p in zip(predictors, fit["beta"], fit["se"], fit["p"]):
            index = predictors.index(predictor)
            results.append({"model": model, "outcome": outcome, "predictor": predictor, "beta": beta, "se": se, "p": p, "robust_se": fit["robust_se"][index], "robust_p": fit["robust_p"][index], "r2": fit["r2"]})
        moran_value, moran_p = moran_permutation(fit["residual"], weights)
        diagnostics.append({"model": model, "diagnostic": "Moran_I_residual", "value": moran_value, "p": moran_p})
    boot = bootstrap_paths(data)
    for key, values in boot.items():
        finite = values[np.isfinite(values)]
        results.append({"model": "derived", "outcome": key, "predictor": "derived", "beta": np.nanmean(finite), "se": np.nanstd(finite, ddof=1), "p": np.nan, "r2": np.nan, "ci_low": np.nanpercentile(finite, 2.5), "ci_high": np.nanpercentile(finite, 97.5)})

    # Coordinates are optional for this first spike; use site order only if available in SI.
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(OUT / "path_coefficients.csv", index=False)
    vif_predictors = ["log_NO3_N", "HSCI", "TOC", "COD", "TP"]
    diagnostics.extend({"model": "predictors", "diagnostic": row["predictor"] + "_VIF", "value": row["vif"], "p": np.nan} for row in calculate_vif(data, vif_predictors))
    pd.DataFrame(diagnostics).to_csv(OUT / "diagnostics.csv", index=False)
    data.to_csv(OUT / "site_model_data.csv", index=False)
    summary = {
        "source_workbook": str(si_path),
        "n_sites": len(data),
        "n_high_narnas_arg_host_sgbs": high_sgb_count,
        "bootstrap_replicates": N_BOOT,
        "seed": SEED,
        "no3_beta_arg_total": float(fitted["arg_total_model"]["beta"][0]),
        "no3_beta_arg_host_adjusted": float(fitted["arg_host_adjusted_model"]["beta"][0]),
        "no3_beta_host": float(fitted["host_model"]["beta"][0]),
        "host_beta_arg_adjusted": float(fitted["arg_host_adjusted_model"]["beta"][1]),
        "bootstrap_indirect_association_ci": [float(np.nanpercentile(boot["indirect_association"], q)) for q in (2.5, 97.5)],
        "bootstrap_no3_attenuation_ci": [float(np.nanpercentile(boot["attenuation"], q)) for q in (2.5, 97.5)],
    }
    for model, fit in fitted.items():
        summary[f"moran_i_{model}_residual"] = float(next(row["value"] for row in diagnostics if row["model"] == model and row["diagnostic"] == "Moran_I_residual"))
        summary[f"moran_p_{model}_residual"] = float(next(row["p"] for row in diagnostics if row["model"] == model and row["diagnostic"] == "Moran_I_residual"))
    pd.Series(summary).to_json(OUT / "summary.json", force_ascii=False, indent=2)
    print(pd.DataFrame(results).to_string(index=False))
    print("\nSummary:")
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
