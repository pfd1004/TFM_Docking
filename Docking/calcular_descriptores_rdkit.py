from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

sdf_dir = Path("Ligandos_sdf")
rows = []

for sdf in sdf_dir.glob("*.sdf"):
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)

    for mol in suppl:
        if mol is None:
            print(f"No se pudo leer: {sdf.name}")
            continue

        ligando = sdf.stem

        rows.append({
            "ligando": ligando,
            "MolWt": Descriptors.MolWt(mol),
            "ExactMolWt": Descriptors.ExactMolWt(mol),
            "MolLogP": Descriptors.MolLogP(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
            "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
            "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
            "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
            "NumRings": rdMolDescriptors.CalcNumRings(mol),
            "FormalCharge": Chem.GetFormalCharge(mol),
        })

df = pd.DataFrame(rows)
df.to_csv("descriptores_rdkit.csv", index=False)
df.to_excel("descriptores_rdkit.xlsx", index=False)

print("Creados descriptores_rdkit.csv y descriptores_rdkit.xlsx")