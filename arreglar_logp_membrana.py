from pathlib import Path, PureWindowsPath
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

XML_FILE = Path("datos_sergio/membrana/DMPC_Permeability.txt")
FEATURES_FILE = Path("membrana_features.csv")

def numbers_from_text(text):
    if text is None:
        return []
    nums = re.findall(
        r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?",
        text
    )
    return [float(x) for x in nums]

def clean_name(name):
    if name is None:
        return ""
    name = str(name).strip().strip('"').strip("'")
    name = name.replace("\\\\", "\\")
    name = PureWindowsPath(name).name

    for ext in [".cosmo", ".orcacosmo", ".txt", ".sdf", ".pdbqt"]:
        if name.lower().endswith(ext):
            name = name[:-len(ext)]

    return name

def get_child_text(element, tag):
    for child in element:
        if child.tag.lower() == tag.lower():
            return child.text
    return None

def get_first_number(element, tag):
    vals = numbers_from_text(get_child_text(element, tag))
    return vals[0] if vals else np.nan

def get_vector(element, tag):
    return numbers_from_text(get_child_text(element, tag))

def summarize_vector(prefix, values, depths):
    row = {}

    if not values:
        row[f"{prefix}_mean"] = np.nan
        row[f"{prefix}_min"] = np.nan
        row[f"{prefix}_max"] = np.nan
        row[f"{prefix}_center"] = np.nan
        row[f"{prefix}_argmax_depth"] = np.nan
        row[f"{prefix}_argmin_depth"] = np.nan
        return row

    arr = np.array(values, dtype=float)

    row[f"{prefix}_mean"] = float(np.nanmean(arr))
    row[f"{prefix}_min"] = float(np.nanmin(arr))
    row[f"{prefix}_max"] = float(np.nanmax(arr))

    if depths is not None and len(depths) == len(arr):
        depths_arr = np.array(depths, dtype=float)
        center_idx = int(np.argmin(np.abs(depths_arr - 0.0)))
        row[f"{prefix}_center"] = float(arr[center_idx])
        row[f"{prefix}_argmax_depth"] = float(depths_arr[int(np.nanargmax(arr))])
        row[f"{prefix}_argmin_depth"] = float(depths_arr[int(np.nanargmin(arr))])
    else:
        row[f"{prefix}_center"] = float(arr[0])
        row[f"{prefix}_argmax_depth"] = np.nan
        row[f"{prefix}_argmin_depth"] = np.nan

    return row

def norm_name(s):
    s = str(s).lower()
    s = s.replace(".cosmo", "")
    s = s.replace(".orcacosmo", "")
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s

# Leer XML
raw = XML_FILE.read_bytes()

for enc in ["windows-1252", "utf-8-sig", "utf-8", "latin-1"]:
    try:
        text = raw.decode(enc)
        break
    except UnicodeDecodeError:
        pass

root = ET.fromstring(text)

# Layer positions globales
depths = []
for el in root.iter():
    if el.tag.lower() == "layerposition":
        depths = numbers_from_text(el.text)
        break

rows = []

for solute in root.iter():
    if solute.tag.lower() != "solute":
        continue

    input_name = get_child_text(solute, "input_name")
    name = input_name if input_name else get_child_text(solute, "name")
    ligando = clean_name(name)

    diffusion = get_vector(solute, "diffusion_E9_m2_s")

    row = {
        "ligando": ligando,
        "mem_logP": get_first_number(solute, "logP_micelle_water"),
        "mem_logP_kg_l": get_first_number(solute, "logP_micelle_water_kg_l"),
        "mem_logPerm": get_first_number(solute, "logPerm_membrane_cm_s"),
    }

    row.update(summarize_vector("mem_diffusion_E9", diffusion, depths))

    rows.append(row)

new = pd.DataFrame(rows)
new["key_ligando"] = new["ligando"].apply(norm_name)

old = pd.read_csv(FEATURES_FILE)
old["key_ligando"] = old["ligando"].apply(norm_name)

# Quitar columnas antiguas mal extraídas
cols_to_replace = [
    "mem_logP",
    "mem_logP_kg_l",
    "mem_logPerm",
    "mem_diffusion_E9_mean",
    "mem_diffusion_E9_min",
    "mem_diffusion_E9_max",
    "mem_diffusion_E9_center",
    "mem_diffusion_E9_argmax_depth",
    "mem_diffusion_E9_argmin_depth",
]

old = old.drop(columns=[c for c in cols_to_replace if c in old.columns])

merged = old.merge(
    new.drop(columns=["ligando"]),
    on="key_ligando",
    how="left"
)

merged = merged.drop(columns=["key_ligando"])

merged.to_csv("membrana_features_corregido.csv", index=False)
merged.to_excel("membrana_features_corregido.xlsx", index=False)

print("Creado membrana_features_corregido.csv")
print("Creado membrana_features_corregido.xlsx")

print("\nValores no vacíos:")
print(merged[["mem_logP", "mem_logP_kg_l", "mem_logPerm", "mem_diffusion_E9_mean"]].notna().sum())

print("\nPrimeras filas:")
print(merged[["ligando", "mem_logP", "mem_logP_kg_l", "mem_logPerm", "mem_diffusion_E9_mean"]].head(10).to_string(index=False))