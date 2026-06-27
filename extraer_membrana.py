from pathlib import Path, PureWindowsPath
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd

INPUT = Path("datos_sergio/membrana/DMPC_Permeability.txt")
OUT = Path("membrana_features.csv")

def numbers_from_text(text):
    if text is None:
        return []
    nums = re.findall(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?", text)
    return [float(x) for x in nums]

def clean_name(name):
    if name is None:
        return ""

    name = str(name).strip().strip('"').strip("'")

    # Arreglar dobles barras si aparecen
    name = name.replace("\\\\", "\\")

    # Si viene como ruta de Windows, quedarnos solo con el archivo final
    name = PureWindowsPath(name).name

    # Quitar extensiones
    for ext in [".cosmo", ".orcacosmo", ".txt", ".sdf", ".pdbqt"]:
        if name.lower().endswith(ext):
            name = name[: -len(ext)]

    return name

def find_child_text(element, tag):
    tag = tag.lower()
    for child in element.iter():
        if child.tag.lower() == tag:
            return child.text
    return None

def get_first_number(element, tag):
    vals = numbers_from_text(find_child_text(element, tag))
    return vals[0] if vals else np.nan

def get_vector(element, tag):
    return numbers_from_text(find_child_text(element, tag))

def summarize_vector(prefix, values, depths=None):
    d = {}
    if not values:
        d[f"{prefix}_mean"] = np.nan
        d[f"{prefix}_min"] = np.nan
        d[f"{prefix}_max"] = np.nan
        d[f"{prefix}_center"] = np.nan
        d[f"{prefix}_argmax_depth"] = np.nan
        d[f"{prefix}_argmin_depth"] = np.nan
        return d

    arr = np.array(values, dtype=float)

    d[f"{prefix}_mean"] = float(np.nanmean(arr))
    d[f"{prefix}_min"] = float(np.nanmin(arr))
    d[f"{prefix}_max"] = float(np.nanmax(arr))

    if depths is not None and len(depths) == len(arr):
        depths_arr = np.array(depths, dtype=float)
        center_idx = int(np.argmin(np.abs(depths_arr - 0.0)))
        d[f"{prefix}_center"] = float(arr[center_idx])
        d[f"{prefix}_argmax_depth"] = float(depths_arr[int(np.nanargmax(arr))])
        d[f"{prefix}_argmin_depth"] = float(depths_arr[int(np.nanargmin(arr))])
    else:
        d[f"{prefix}_center"] = float(arr[0])
        d[f"{prefix}_argmax_depth"] = np.nan
        d[f"{prefix}_argmin_depth"] = np.nan

    return d

text = INPUT.read_text(errors="ignore")

# Por si tiene algún carácter raro antes/después del XML
start = text.find("<")
end = text.rfind(">")
if start != -1 and end != -1:
    text_xml = text[start:end+1]
else:
    text_xml = text

root = ET.fromstring(text_xml)

# Profundidades globales, si existen fuera de cada solute
global_depths = []
for el in root.iter():
    if el.tag.lower() == "layerposition":
        global_depths = numbers_from_text(el.text)
        break

rows = []

for solute in root.iter():
    if solute.tag.lower() != "solute":
        continue

    # Nombre del soluto: puede estar como atributo o como texto/hijo
    name = (
        solute.attrib.get("name")
        or solute.attrib.get("id")
        or solute.attrib.get("compound")
    )

    if not name:
        # Si <solute> contiene texto directo
        direct = solute.text.strip() if solute.text else ""
        if direct and not re.search(r"\d", direct):
            name = direct

    if not name:
        # Buscar tags internas típicas
        for possible in ["name", "compound", "filename"]:
            val = find_child_text(solute, possible)
            if val:
                name = val.strip()
                break

    name = clean_name(name)

    depths = get_vector(solute, "layerposition")
    if not depths:
        depths = global_depths

    distribution = get_vector(solute, "distribution")
    free_energy = get_vector(solute, "free_energy")
    meanentropy = get_vector(solute, "meanentropy")
    diffusion = get_vector(solute, "diffusion_E9")

    row = {
        "ligando": name,
        "mem_logP": get_first_number(solute, "logP"),
        "mem_logPerm": get_first_number(solute, "logPerm"),
        "n_layerposition": len(depths),
    }

    row.update(summarize_vector("mem_distribution", distribution, depths))
    row.update(summarize_vector("mem_free_energy", free_energy, depths))
    row.update(summarize_vector("mem_meanentropy", meanentropy, depths))
    row.update(summarize_vector("mem_diffusion_E9", diffusion, depths))

    # Variables más interpretables
    if distribution:
        arr = np.array(distribution, dtype=float)
        row["mem_distribution_sum"] = float(np.nansum(arr))
        row["mem_distribution_nonzero_center"] = int(row["mem_distribution_center"] > 0)
    else:
        row["mem_distribution_sum"] = np.nan
        row["mem_distribution_nonzero_center"] = np.nan

    if free_energy:
        arr = np.array(free_energy, dtype=float)
        row["mem_free_energy_barrier"] = float(np.nanmax(arr) - np.nanmin(arr))
    else:
        row["mem_free_energy_barrier"] = np.nan

    rows.append(row)

df = pd.DataFrame(rows)

# Quitar solutos sin nombre, si los hubiera
df = df[df["ligando"].astype(str).str.len() > 0].copy()

df.to_csv(OUT, index=False)
df.to_excel("membrana_features.xlsx", index=False)

print(f"Creado {OUT}")
print(f"Solutos detectados: {len(df)}")
print(df.head())