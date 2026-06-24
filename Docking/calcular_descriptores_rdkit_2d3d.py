from pathlib import Path
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Descriptors3D

sdf_dir = Path("Ligandos_sdf")
rows = []

for sdf in sdf_dir.glob("*.sdf"):
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)

    for mol in suppl:
        if mol is None:
            print(f"No se pudo leer: {sdf.name}")
            continue

        ligando = sdf.stem

        row = {
            "ligando": ligando,

            # 2D básicos
            "MolWt": Descriptors.MolWt(mol),
            "ExactMolWt": Descriptors.ExactMolWt(mol),
            "MolLogP": Descriptors.MolLogP(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol),
            "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
            "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
            "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "NumHeavyAtoms": mol.GetNumHeavyAtoms(),
            "NumAtoms": mol.GetNumAtoms(),
            "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
            "NumRings": rdMolDescriptors.CalcNumRings(mol),
            "FormalCharge": Chem.GetFormalCharge(mol),

            # Fingerprint-like counts
            "Chi0n": Descriptors.Chi0n(mol),
            "Chi1n": Descriptors.Chi1n(mol),
            "Kappa1": Descriptors.Kappa1(mol),
            "Kappa2": Descriptors.Kappa2(mol),
            "LabuteASA": rdMolDescriptors.CalcLabuteASA(mol),
        }

        # 3D: solo si hay conformación 3D
        if mol.GetNumConformers() > 0:
            try:
                row.update({
                    "Asphericity": Descriptors3D.Asphericity(mol),
                    "Eccentricity": Descriptors3D.Eccentricity(mol),
                    "InertialShapeFactor": Descriptors3D.InertialShapeFactor(mol),
                    "NPR1": Descriptors3D.NPR1(mol),
                    "NPR2": Descriptors3D.NPR2(mol),
                    "PMI1": Descriptors3D.PMI1(mol),
                    "PMI2": Descriptors3D.PMI2(mol),
                    "PMI3": Descriptors3D.PMI3(mol),
                    "RadiusOfGyration": Descriptors3D.RadiusOfGyration(mol),
                    "SpherocityIndex": Descriptors3D.SpherocityIndex(mol),
                })
            except Exception as e:
                print(f"Error 3D en {sdf.name}: {e}")
        else:
            print(f"Sin conformación 3D: {sdf.name}")

        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv("descriptores_rdkit_2d3d.csv", index=False)
df.to_excel("descriptores_rdkit_2d3d.xlsx", index=False)

print("Creados:")
print("descriptores_rdkit_2d3d.csv")
print("descriptores_rdkit_2d3d.xlsx")