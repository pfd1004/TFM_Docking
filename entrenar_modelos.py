from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import joblib

from sklearn.model_selection import train_test_split, KFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    make_scorer,
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)
from sklearn.linear_model import Ridge
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor


# =========================
# CONFIGURACIÓN
# =========================

INPUT_FILE = "dataset_descriptores.csv"

TARGETS = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

TEST_SIZE = 0.20
RANDOM_STATE = 42

OUTDIR = Path("resultados_modelos")
OUTDIR.mkdir(exist_ok=True)


# =========================
# MÉTRICAS
# =========================

def rmse_multi(y_true, y_pred):
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


rmse_scorer = make_scorer(rmse_multi, greater_is_better=False)


def calcular_metricas(modelo, y_true, y_pred):
    filas = []

    for i, target in enumerate(TARGETS):
        filas.append({
            "modelo": modelo,
            "salida": target,
            "MAE": mean_absolute_error(y_true[:, i], y_pred[:, i]),
            "RMSE": np.sqrt(mean_squared_error(y_true[:, i], y_pred[:, i])),
            "R2": r2_score(y_true[:, i], y_pred[:, i]),
        })

    filas.append({
        "modelo": modelo,
        "salida": "media",
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "R2": r2_score(y_true, y_pred, multioutput="uniform_average"),
    })

    return filas


# =========================
# CARGA DE DATOS
# =========================

df = pd.read_csv(INPUT_FILE)

print("Columnas del dataset:")
print(list(df.columns))

for target in TARGETS:
    if target not in df.columns:
        raise ValueError(f"Falta la columna objetivo: {target}")

df = df.dropna(subset=TARGETS).copy()

# X = descriptores
# Y = afinidades frente a las 4 proteínas
X = df.drop(columns=["ligando"] + TARGETS)
X = X.select_dtypes(include=[np.number])
X = X.replace([np.inf, -np.inf], np.nan)

# Quitar columnas con más de 30% de NaN
keep_cols = X.columns[X.isna().mean() <= 0.30]
X = X[keep_cols]

Y = df[TARGETS]

print("\nDataset preparado:")
print("Ligandos:", len(df))
print("Descriptores usados:", X.shape[1])
print("Salidas:", TARGETS)

if len(df) < 30:
    print("\nAVISO: dataset pequeño; las métricas pueden ser inestables.")


# =========================
# TRAIN / TEST
# =========================

X_train, X_test, Y_train, Y_test, lig_train, lig_test = train_test_split(
    X,
    Y,
    df["ligando"],
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
)

Y_train_np = Y_train.values
Y_test_np = Y_test.values


# =========================
# MODELOS + GRIDSEARCH
# =========================

max_pls = max(1, min(6, len(X_train) - 1, X_train.shape[1]))

modelos = {
    "Ridge": (
        Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold()),
            ("scaler", StandardScaler()),
            ("model", Ridge()),
        ]),
        {
            "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0],
        }
    ),

    "PLSRegression": (
        Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold()),
            ("scaler", StandardScaler()),
            ("model", PLSRegression()),
        ]),
        {
            "model__n_components": list(range(1, max_pls + 1)),
        }
    ),

    "RandomForest": (
        Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold()),
            ("model", RandomForestRegressor(
                random_state=RANDOM_STATE,
                n_jobs=-1,
            )),
        ]),
        {
            "model__n_estimators": [300, 600],
            "model__max_features": ["sqrt", 0.5, 1.0],
            "model__min_samples_leaf": [1, 2, 4],
            "model__max_depth": [None, 5, 10],
        }
    ),

    "ExtraTrees": (
        Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold()),
            ("model", ExtraTreesRegressor(
                random_state=RANDOM_STATE,
                n_jobs=-1,
            )),
        ]),
        {
            "model__n_estimators": [300, 600],
            "model__max_features": ["sqrt", 0.5, 1.0],
            "model__min_samples_leaf": [1, 2, 4],
            "model__max_depth": [None, 5, 10],
        }
    ),
}


cv = KFold(n_splits=min(5, len(X_train)), shuffle=True, random_state=RANDOM_STATE)

resumen = []
metricas = []
predicciones = []

mejor_modelo = None
mejor_nombre = None
mejor_rmse = np.inf


# =========================
# ENTRENAMIENTO
# =========================

for nombre, (pipeline, param_grid) in modelos.items():
    print(f"\nEntrenando con GridSearchCV: {nombre}")

    grid = GridSearchCV(
        estimator=pipeline,
        param_grid=param_grid,
        scoring=rmse_scorer,
        cv=cv,
        refit=True,
        n_jobs=1,
        verbose=1,
    )

    grid.fit(X_train, Y_train_np)

    pred = grid.predict(X_test)
    pred = np.asarray(pred)

    test_mae = mean_absolute_error(Y_test_np, pred)
    test_rmse = np.sqrt(mean_squared_error(Y_test_np, pred))
    test_r2 = r2_score(Y_test_np, pred, multioutput="uniform_average")

    resumen.append({
        "modelo": nombre,
        "CV_RMSE": -grid.best_score_,
        "TEST_MAE": test_mae,
        "TEST_RMSE": test_rmse,
        "TEST_R2": test_r2,
        "best_params": grid.best_params_,
    })

    metricas.extend(calcular_metricas(nombre, Y_test_np, pred))

    pred_df = pd.DataFrame({
        "modelo": nombre,
        "ligando": list(lig_test),
    })

    for i, target in enumerate(TARGETS):
        pred_df[f"real_{target}"] = Y_test_np[:, i]
        pred_df[f"pred_{target}"] = pred[:, i]
        pred_df[f"error_{target}"] = pred[:, i] - Y_test_np[:, i]

    predicciones.append(pred_df)

    pd.DataFrame(grid.cv_results_).to_csv(
        OUTDIR / f"cv_results_{nombre}.csv",
        index=False
    )

    joblib.dump(grid.best_estimator_, OUTDIR / f"best_{nombre}.joblib")

    print("Mejores parámetros:", grid.best_params_)
    print(f"CV RMSE: {-grid.best_score_:.4f}")
    print(f"TEST RMSE: {test_rmse:.4f}")
    print(f"TEST MAE: {test_mae:.4f}")
    print(f"TEST R2: {test_r2:.4f}")

    if test_rmse < mejor_rmse:
        mejor_rmse = test_rmse
        mejor_nombre = nombre
        mejor_modelo = grid.best_estimator_


# =========================
# GUARDAR RESULTADOS
# =========================

resumen_df = pd.DataFrame(resumen)
metricas_df = pd.DataFrame(metricas)
predicciones_df = pd.concat(predicciones, ignore_index=True)

resumen_df.to_csv(OUTDIR / "resumen_modelos.csv", index=False)
metricas_df.to_csv(OUTDIR / "metricas_por_proteina.csv", index=False)
predicciones_df.to_csv(OUTDIR / "predicciones_test.csv", index=False)

joblib.dump(mejor_modelo, OUTDIR / "mejor_modelo_global.joblib")

with open(OUTDIR / "mejor_modelo_global.txt", "w") as f:
    f.write(f"Mejor modelo por TEST_RMSE: {mejor_nombre}\n")
    f.write(f"TEST_RMSE: {mejor_rmse:.4f}\n")

print("\nTerminado.")
print(f"Mejor modelo: {mejor_nombre}")
print(f"TEST RMSE: {mejor_rmse:.4f}")
print(f"Resultados guardados en: {OUTDIR}")