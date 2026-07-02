
import pandas as pd
import numpy as np
import joblib
from pathlib import Path
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import ExtraTreesRegressor

INPUT = "dataset_descriptores_membrana.csv"
OUTDIR = Path("modelo_final_rdkit")
OUTDIR.mkdir(exist_ok=True)
TARGETS = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

df = pd.read_csv(INPUT).dropna(subset=TARGETS).copy()
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
feature_cols = [
    c for c in numeric_cols
    if c not in TARGETS and not c.startswith("mem_") and c != "n_layerposition"
]
X = df[feature_cols].replace([np.inf, -np.inf], np.nan)
Y = df[TARGETS]

model = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("variance", VarianceThreshold()),
    ("model", ExtraTreesRegressor(
        n_estimators=300,
        max_depth=None,
        max_features=0.5,
        min_samples_leaf=1,
        random_state=42,
        n_jobs=-1,
    ))
])
model.fit(X, Y)
joblib.dump(model, OUTDIR / "modelo_final_extratrees_rdkit.joblib")
pd.DataFrame({"feature": feature_cols}).to_csv(OUTDIR / "features_modelo_final.csv", index=False)
with open(OUTDIR / "descripcion_modelo_final.txt", "w", encoding="utf-8") as f:
    f.write("Modelo final: ExtraTreesRegressor multisalida\n")
    f.write("Entradas: descriptores RDKit 2D/3D\n")
    f.write("Salidas: y_1ao6, y_1hzh, y_2hav, y_3ghg\n")
    f.write("n_estimators=300; max_depth=None; max_features=0.5; min_samples_leaf=1\n")
print("Modelo final entrenado.")
print("N ligandos:", len(df))
print("N features:", len(feature_cols))
print("Guardado en modelo_final_rdkit/modelo_final_extratrees_rdkit.joblib")
