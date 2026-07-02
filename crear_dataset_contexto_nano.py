import pandas as pd
import numpy as np
from pathlib import Path

# =========================
# ENTRADAS
# =========================

DATASET = "dataset_descriptores_membrana.csv"
NANO = "resultados_docking_nanoparticulas_corregido.csv"

OUT_LONG = "dataset_ligando_nano_proteina_contexto_long.csv"
OUT_RANKING = "ranking_contextual_ligando_nano_proteina.csv"
OUT_NANO_WIDE = "tabla_nano_proteina_contexto_wide.csv"

TARGETS = {
    "1ao6": "y_1ao6",
    "1hzh": "y_1hzh",
    "2hav": "y_2hav",
    "3ghg": "y_3ghg",
}


# =========================
# FUNCIONES
# =========================

def zscore(x):
    x = pd.Series(x, dtype=float)
    sd = x.std(ddof=0)

    if pd.isna(sd) or sd == 0:
        return pd.Series(np.zeros(len(x)), index=x.index)

    return (x - x.mean()) / sd


# =========================
# CARGA
# =========================

df = pd.read_csv(DATASET)
nano = pd.read_csv(NANO)

print("Dataset ligandos:", df.shape)
print("Dataset nano-proteína:", nano.shape)

# Comprobar columnas y
for p, ycol in TARGETS.items():
    if ycol not in df.columns:
        raise ValueError(f"Falta {ycol} en {DATASET}")

if "binding_energy_min" not in nano.columns:
    raise ValueError("Falta binding_energy_min en la tabla de nanopartículas.")

# Limpiar nano
nano = nano.copy()
nano["proteina"] = nano["proteina"].astype(str).str.lower()
nano = nano[nano["estado"].astype(str).str.lower().eq("ok")].copy()

nano = nano[[
    "nanoparticula",
    "proteina",
    "binding_energy_min",
    "n_energias_detectadas",
    "archivo",
]].rename(columns={
    "binding_energy_min": "E_nano_proteina",
    "archivo": "archivo_dlg_nano_proteina",
})

# =========================
# 1. PASAR LIGANDOS A FORMATO LONG
# =========================

feature_cols = [
    c for c in df.columns
    if c not in ["ligando"] + list(TARGETS.values())
]

long_rows = []

for _, row in df.iterrows():
    for proteina, ycol in TARGETS.items():
        base = {
            "ligando": row["ligando"],
            "proteina": proteina,
            "y_ligando_proteina": row[ycol],
        }

        # Añadimos descriptores del ligando y membrana
        for c in feature_cols:
            base[c] = row[c]

        long_rows.append(base)

lig_long = pd.DataFrame(long_rows)

print("Tabla ligando-proteína long:", lig_long.shape)

# =========================
# 2. AÑADIR CONTEXTO NANO-PROTEÍNA
# =========================

# Al hacer merge por proteína, cada ligando-proteína se replica para cada nanopartícula
data = lig_long.merge(
    nano,
    on="proteina",
    how="left"
)

print("Tabla ligando-nano-proteína:", data.shape)

# =========================
# 3. SCORES NORMALIZADOS
# =========================

# Más negativo = mejor.
# Por eso usamos -zscore: valores más favorables quedan más altos.

data["z_ligando_proteina"] = data.groupby("proteina")["y_ligando_proteina"].transform(zscore)
data["score_ligando_proteina"] = -data["z_ligando_proteina"]

# Para nano, calculamos z-score en la tabla única nano-proteína
nano_unique = nano.copy()
nano_unique["z_nano_proteina"] = nano_unique.groupby("proteina")["E_nano_proteina"].transform(zscore)
nano_unique["score_nano_proteina"] = -nano_unique["z_nano_proteina"]

nano_unique = nano_unique[[
    "nanoparticula",
    "proteina",
    "z_nano_proteina",
    "score_nano_proteina",
]]

data = data.merge(
    nano_unique,
    on=["nanoparticula", "proteina"],
    how="left"
)

# Score contextual integrado.
# No suma energías brutas; suma puntuaciones normalizadas.
data["score_contextual"] = (
    data["score_ligando_proteina"] +
    data["score_nano_proteina"]
)

# También dejamos una versión media
data["score_contextual_medio"] = data["score_contextual"] / 2

# =========================
# 4. GUARDAR DATASET LONG
# =========================

data.to_csv(OUT_LONG, index=False)
data.to_excel("dataset_ligando_nano_proteina_contexto_long.xlsx", index=False)

# =========================
# 5. RANKING
# =========================

ranking = data.dropna(subset=["E_nano_proteina"]).copy()

ranking = ranking.sort_values(
    ["proteina", "score_contextual"],
    ascending=[True, False]
)

ranking_cols = [
    "ligando",
    "nanoparticula",
    "proteina",
    "y_ligando_proteina",
    "E_nano_proteina",
    "score_ligando_proteina",
    "score_nano_proteina",
    "score_contextual",
    "score_contextual_medio",
]

ranking[ranking_cols].to_csv(OUT_RANKING, index=False)
ranking[ranking_cols].to_excel("ranking_contextual_ligando_nano_proteina.xlsx", index=False)

# =========================
# 6. TABLA NANO-PROTEÍNA WIDE
# =========================

nano_wide = nano.pivot_table(
    index="nanoparticula",
    columns="proteina",
    values="E_nano_proteina",
    aggfunc="first"
).reset_index()

nano_wide.to_csv(OUT_NANO_WIDE, index=False)
nano_wide.to_excel("tabla_nano_proteina_contexto_wide.xlsx", index=False)

# =========================
# RESUMEN
# =========================

print("\nCreados:")
print(OUT_LONG)
print("dataset_ligando_nano_proteina_contexto_long.xlsx")
print(OUT_RANKING)
print("ranking_contextual_ligando_nano_proteina.xlsx")
print(OUT_NANO_WIDE)
print("tabla_nano_proteina_contexto_wide.xlsx")

print("\nTabla nano-proteína:")
print(nano_wide.to_string(index=False))

print("\nTop 30 ranking contextual:")
print(ranking[ranking_cols].head(30).to_string(index=False))

print("\nAviso:")
print("Este dataset sirve para análisis contextual/ranking, no como sustituto directo del modelo principal.")
print("Para entrenar un modelo real nano-ligando-proteína harían falta y reales de complejos nano-ligando-proteína.")