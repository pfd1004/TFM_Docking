import pandas as pd

dock = pd.read_csv("tabla_resultados_docking_final.csv")
desc = pd.read_csv("descriptores_rdkit_2d3d.csv")

# Solo dockings correctos
dock = dock[dock["estado"] == "ok"].copy()

# Pasar de formato largo a ancho
wide = dock.pivot_table(
    index="ligando",
    columns="proteina",
    values="afinidad_kcal_mol",
    aggfunc="first"
).reset_index()

# Renombrar salidas
wide = wide.rename(columns={
    "1ao6": "y_1ao6",
    "1hzh": "y_1hzh",
    "2hav": "y_2hav",
    "3ghg": "y_3ghg",
})

targets = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

# Quitar ligandos incompletos
wide_complete = wide.dropna(subset=targets).copy()

# Unir con descriptores
df = desc.merge(wide_complete, on="ligando", how="inner")

# Revisar
print("Ligandos con las 4 salidas:", len(df))
print("Columnas:", len(df.columns))

df.to_csv("dataset_multisalida_docking_descriptores.csv", index=False)
df.to_excel("dataset_multisalida_docking_descriptores.xlsx", index=False)

print("Creados:")
print("dataset_multisalida_docking_descriptores.csv")
print("dataset_multisalida_docking_descriptores.xlsx")