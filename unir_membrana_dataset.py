import re
import pandas as pd

DATASET = "dataset_descriptores.csv"
MEMBRANA = "membrana_features_corregido.csv"

def norm_name(s):
    s = str(s).lower()
    s = s.replace(".cosmo", "")
    s = s.replace(".orcacosmo", "")
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s

df = pd.read_csv(DATASET)
mem = pd.read_csv(MEMBRANA)

df["key_ligando"] = df["ligando"].apply(norm_name)
mem["key_ligando"] = mem["ligando"].apply(norm_name)

merged = df.merge(
    mem.drop(columns=["ligando"]),
    on="key_ligando",
    how="left"
)

sin_mem = merged[merged["mem_logP"].isna()]["ligando"].unique()

print("Filas dataset:", len(df))
print("Filas tras unión:", len(merged))
print("Ligandos sin datos de membrana:", len(sin_mem))

if len(sin_mem) > 0:
    print("\nLigandos sin membrana:")
    for x in sin_mem:
        print(" -", x)

merged = merged.drop(columns=["key_ligando"])

merged.to_csv("dataset_descriptores_membrana.csv", index=False)
merged.to_excel("dataset_descriptores_membrana.xlsx", index=False)

print("\nCreados:")
print("dataset_descriptores_membrana.csv")
print("dataset_descriptores_membrana.xlsx")