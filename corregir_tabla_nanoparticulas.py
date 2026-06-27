import pandas as pd
import re

df = pd.read_csv("resultados_docking_nanoparticulas.csv")

proteins = ["1ao6", "1hzh", "2hav", "3ghg"]

nanos = []
prots = []

for archivo in df["archivo"]:
    stem = archivo.replace(".dlg", "")

    found = None
    for p in proteins:
        if stem.lower().endswith("_" + p):
            found = p
            break

    if found is None:
        nanos.append(stem)
        prots.append("")
    else:
        nano = stem[:-(len(found)+1)]
        nanos.append(nano)
        prots.append(found)

df["nanoparticula"] = nanos
df["proteina"] = prots

df.to_csv("resultados_docking_nanoparticulas_corregido.csv", index=False)
df.to_excel("resultados_docking_nanoparticulas_corregido.xlsx", index=False)

print(df)
print("Creado resultados_docking_nanoparticulas_corregido.csv")