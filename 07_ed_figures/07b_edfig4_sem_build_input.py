# -*- coding: utf-8 -*-
"""从 SI Tables 构建 SEM 输入数据

变量映射:
- 样本标识: Table 5 样本列 (站点内复制用 _N 后缀区分)
- W1-W15 = Table 7 B-P 列 (pH, DO, Temperature, TDS, EC, Turbidity, COD,
  TC, TOC, IC, TN, NO3_N, NO2_N, NH3_N, TP); 模型使用其中 13 个
  (跳过 TC/IC, 因 TOC+IC 冗余)
- H8-H11 = Table 7 X-AA 列 (NTL, PISA, NDVI, HSCI)
- Microbe   = Table 5 每样本物种丰度总和 (四域)
- ARG       = Table 4 每样本 ARG 丰度总和
- Microbe_m = Table 5 中属于 Table 6 人源微生物列表的物种丰度和
- Microbe_p = Table 5 中属于 Table 6 病原微生物列表的物种丰度和
- group     = 站点名 (样本名去 _N 复制编号; QSH -> QSH2 对齐 Table 1 命名)

注意: Table 4/7 样本列序与 Table 5 不同, 必须按样本名 reindex, 不能按位置合并。
"""
import sys, re, pathlib, yaml
sys.stdout.reconfigure(encoding="utf-8")
import openpyxl
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
with open(ROOT / "config" / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)
SI = ROOT / cfg["si_tables_xlsx"]
wb = openpyxl.load_workbook(SI, read_only=True)


def read_table(name, skip=1):
    ws = wb[name]
    it = ws.iter_rows(values_only=True)
    for _ in range(skip):
        next(it)
    hdr = [str(h) if h else "" for h in next(it)]
    rows = []
    for r in it:
        if r[0] is None:
            continue
        rows.append(r)
    return pd.DataFrame(rows, columns=hdr)


t4 = read_table("Table 4")
t5 = read_table("Table 5")
t7 = read_table("Table 7")
t6 = read_table("Table 6")

# 样本标识以 Table 5 列名为准 (与 kraken.xlsx 一致)
samples = [str(x) for x in t5.columns[7:]]
print(f"Table 5 样本数: {len(samples)}")

# --- Microbe / Microbe_m / Microbe_p (Table 5, 四域) ---
t5n = t5.iloc[:, 7:].apply(pd.to_numeric, errors="coerce")
microbe = t5n.sum(axis=0)

t5_species = [str(x) for x in t5.iloc[:, 0]]
human_sp = set(str(x) for x in t6.iloc[:, 0] if x)
patho_sp = set(str(x) for x in t6.iloc[:, 1] if x)
print(f"\nTable 5 物种 {len(t5_species)}, 人源列表 {len(human_sp)}, 病原列表 {len(patho_sp)}")
print(f"人源命中: {len(set(t5_species) & human_sp)}, 病原命中: {len(set(t5_species) & patho_sp)}")

microbe_m = t5n[t5.iloc[:, 0].isin(human_sp)].sum(axis=0)
microbe_p = t5n[t5.iloc[:, 0].isin(patho_sp)].sum(axis=0)

# --- ARG = Table 4 每样本 ARG 丰度和 (按样本名 reindex) ---
t4n = t4.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
arg_s = t4n.sum(axis=1)
arg = pd.Series(arg_s.values, index=[str(x) for x in t4.iloc[:, 0]]).reindex(samples)

# --- W1-W15 水质变量 + H8-H11 人类活动变量 (Table 7, 按样本名 reindex) ---
w_cols = ["pH", "DO", "Temperature", "TDS", "EC", "Turbidity", "COD",
          "TC", "TOC", "IC", "TN", "NO3_N", "NO2_N", "NH3_N", "TP"]
w_names = {c: f"W{i+1}" for i, c in enumerate(w_cols)}
env = t7[["Sample"] + w_cols].copy()
env = env.rename(columns={c: w_names[c] for c in w_cols})
h_cols = ["NTL", "PISA", "NDVI", "HSCI"]
for i, c in enumerate(h_cols):
    env[f"H{8 + i}"] = t7[c]
env = env.set_index("Sample").reindex(samples).reset_index(drop=True)
print("\nW 变量映射:", {w_names[c]: c for c in w_cols})
print("H 变量映射:", {f"H{8+i}": c for i, c in enumerate(h_cols)})
print(f"env NaN: W={env.iloc[:, :15].isna().sum().sum()}, H={env.iloc[:, 15:].isna().sum().sum()}")

# --- group: 样本 -> 站点 (Table 1 站点名, 89 组) ---
def sample_to_type(s):
    s = str(s)
    m = re.match(r"^(.*?)_memb$", s)
    if m:
        site = re.sub(r"(_\d+|\d+)$", "", m.group(1))
    else:
        site = re.sub(r"(_\d+|\d+)$", "", s)
    return "QSH2" if site == "QSH" else site

group = [sample_to_type(s) for s in samples]
print(f"\n分组: {len(set(group))} 组, {pd.Series(group).value_counts().to_dict()}")

df = pd.DataFrame({
    "sample": samples,
    "Microbe": microbe.values,
    "ARG": arg.values,
    "Microbe_m": microbe_m.values,
    "Microbe_p": microbe_p.values,
    "group": group,
})
df = pd.concat([df, env], axis=1)
out = ROOT / "07_ed_figures" / "07b_edfig4_sem_input.csv"
df.to_csv(out, index=False)
print(f"\n07b_edfig4_sem_input.csv 已保存: {df.shape}")
print("列:", list(df.columns))
print(df.head(3).to_string())
