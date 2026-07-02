import pandas as pd
import numpy as np

INPUT = "dataset_descriptores_membrana.csv"

df = pd.read_csv(INPUT)

TARGETS = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

print("Dimensiones:", df.shape)
print("\nColumnas objetivo encontradas:")
print([c for c in TARGETS if c in df.columns])

missing_targets = [c for c in TARGETS if c not in df.columns]
if missing_targets:
    raise ValueError(f"Faltan columnas y: {missing_targets}")

numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
feature_cols = [c for c in numeric_cols if c not in TARGETS]

rows = []

for c in feature_cols:
    rows.append({
        "columna": c,
        "missing": int(df[c].isna().sum()),
        "missing_pct": df[c].isna().mean(),
        "n_unique": int(df[c].nunique(dropna=True)),
        "constante": int(df[c].nunique(dropna=True) <= 1),
        "es_membrana": int(c.startswith("mem_") or c == "n_layerposition"),
    })

info = pd.DataFrame(rows)

info.to_csv("diagnostico_columnas.csv", index=False)

print("\nColumnas vacías o constantes:")
print(info[(info["missing_pct"] == 1) | (info["constante"] == 1)][["columna", "missing", "missing_pct", "n_unique"]].to_string(index=False))

# Correlación Spearman media con las 4 salidas
corr_rows = []

for c in feature_cols:
    vals = []

    for y in TARGETS:
        sub = df[[c, y]].dropna()

        if len(sub) > 3 and sub[c].nunique() > 1:
            corr = sub[c].corr(sub[y], method="spearman")
            vals.append(abs(corr))

    if vals:
        corr_rows.append({
            "columna": c,
            "abs_spearman_media_y": float(np.nanmean(vals)),
            "abs_spearman_max_y": float(np.nanmax(vals)),
            "es_membrana": int(c.startswith("mem_") or c == "n_layerposition"),
        })

corr_df = pd.DataFrame(corr_rows).sort_values("abs_spearman_media_y", ascending=False)
corr_df.to_csv("correlaciones_features_targets.csv", index=False)

print("\nTop 25 variables más correlacionadas con las y:")
print(corr_df.head(25).to_string(index=False))

print("\nArchivos creados:")
print("diagnostico_columnas.csv")
print("correlaciones_features_targets.csv")