import pandas as pd

dock = pd.read_csv("tabla_resultados_docking_final.csv")
desc = pd.read_csv("descriptores_rdkit_2d3d.csv")

# Nos quedamos solo con dockings correctos
dock = dock[dock["estado"] == "ok"].copy()

# Unión por nombre de ligando
df = dock.merge(desc, on="ligando", how="left")

# Comprobar ligandos sin descriptores
sin_desc = df[df["MolWt"].isna()]["ligando"].unique()

if len(sin_desc) > 0:
    print("Ligandos sin descriptores:")
    for lig in sin_desc:
        print(" -", lig)

df.to_csv("dataset_modelo_docking_descriptores.csv", index=False)
df.to_excel("dataset_modelo_docking_descriptores.xlsx", index=False)

print("Creados:")
print("dataset_modelo_docking_descriptores.csv")
print("dataset_modelo_docking_descriptores.xlsx")