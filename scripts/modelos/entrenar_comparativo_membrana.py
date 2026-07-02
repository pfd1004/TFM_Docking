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
from sklearn.decomposition import PCA
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

INPUT = "datasets/dataset_descriptores_membrana.csv"

TARGETS = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

OUTDIR = Path("resultados/modelos/resultados_comparacion_membrana")
OUTDIR.mkdir(exist_ok=True)

RANDOM_STATE = 42
TEST_SIZE = 0.20


# =========================
# FUNCIONES
# =========================

def rmse_multi(y_true, y_pred):
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


rmse_scorer = make_scorer(rmse_multi, greater_is_better=False)


def calcular_metricas(nombre_set, nombre_modelo, y_true, y_pred):
    rows = []

    for i, target in enumerate(TARGETS):
        rows.append({
            "feature_set": nombre_set,
            "modelo": nombre_modelo,
            "salida": target,
            "MAE": mean_absolute_error(y_true[:, i], y_pred[:, i]),
            "RMSE": np.sqrt(mean_squared_error(y_true[:, i], y_pred[:, i])),
            "R2": r2_score(y_true[:, i], y_pred[:, i]),
        })

    rows.append({
        "feature_set": nombre_set,
        "modelo": nombre_modelo,
        "salida": "media",
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": np.sqrt(mean_squared_error(y_true, y_pred)),
        "R2": r2_score(y_true, y_pred, multioutput="uniform_average"),
    })

    return rows


def limpiar_columnas(X):
    """
    Quita columnas totalmente vacías o constantes antes de entrenar.
    El resto de NaN se imputan dentro del pipeline.
    """
    X = X.copy()
    X = X.replace([np.inf, -np.inf], np.nan)

    # quitar columnas completamente vacías
    X = X.dropna(axis=1, how="all")

    # quitar columnas constantes
    nunique = X.nunique(dropna=True)
    keep = nunique[nunique > 1].index.tolist()
    X = X[keep]

    return X


# =========================
# CARGA
# =========================

df = pd.read_csv(INPUT)

print("Dataset cargado:", INPUT)
print("Dimensiones iniciales:", df.shape)

for y in TARGETS:
    if y not in df.columns:
        raise ValueError(f"Falta la columna objetivo: {y}")

df = df.dropna(subset=TARGETS).copy()

# Detectar columnas
numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

mem_cols_all = [
    c for c in numeric_cols
    if c.startswith("mem_") or c == "n_layerposition"
]

rdkit_cols_all = [
    c for c in numeric_cols
    if c not in TARGETS and c not in mem_cols_all
]

# Columnas de membrana que NO recomiendo usar porque son vacías,
# constantes o poco informativas en este CSV concreto.
mem_cols_malas = {
    "n_layerposition",
    "mem_meanentropy_max",
    "mem_meanentropy_center",
    "mem_meanentropy_argmax_depth",
    "mem_meanentropy_argmin_depth",
    "mem_distribution_nonzero_center",
    "mem_diffusion_E9_argmax_depth",
    "mem_diffusion_E9_argmin_depth",
}

# Membrana seleccionada: variables interpretables y no constantes
mem_cols_selected = [
    c for c in mem_cols_all
    if c not in mem_cols_malas
]

# Si FormalCharge es constante, tampoco aporta
rdkit_cols_selected = [
    c for c in rdkit_cols_all
    if c != "FormalCharge"
]

print("\nRDKit cols:", len(rdkit_cols_selected))
print("Membrana cols totales:", len(mem_cols_all))
print("Membrana cols seleccionadas:", len(mem_cols_selected))

# Para comparación justa, usamos solo ligandos que tengan membrana real.
# En tu CSV probablemente SDS no tiene membrana.
if "mem_free_energy_barrier" in df.columns:
    sin_mem = df[df["mem_free_energy_barrier"].isna()]["ligando"].tolist()
    if sin_mem:
        print("\nLigandos sin membrana real que se excluirán en esta comparación:")
        for x in sin_mem:
            print(" -", x)

    df_common = df[df["mem_free_energy_barrier"].notna()].copy()
else:
    df_common = df.copy()

print("\nFilas usadas para comparación justa:", len(df_common))

# =========================
# FEATURE SETS
# =========================

feature_sets = {
    "RDKit_solo": rdkit_cols_selected,
    "RDKit_membrana_seleccionada": rdkit_cols_selected + mem_cols_selected,
    "Membrana_solo": mem_cols_selected,
}

# Versión filtrada automática: todas las numéricas menos targets,
# quitando columnas vacías/constantes después.
all_candidate_cols = [
    c for c in numeric_cols
    if c not in TARGETS
]

feature_sets["RDKit_membrana_filtrada_auto"] = all_candidate_cols


# =========================
# ENTRENAMIENTO
# =========================

resumen_rows = []
metricas_rows = []
predicciones_rows = []

best_global_model = None
best_global_name = None
best_global_set = None
best_global_rmse = np.inf

for set_name, cols in feature_sets.items():

    cols = [c for c in cols if c in df_common.columns]

    X = df_common[cols]
    X = limpiar_columnas(X)

    Y = df_common[TARGETS]

    print("\n" + "=" * 70)
    print("Feature set:", set_name)
    print("N ligandos:", len(df_common))
    print("N variables:", X.shape[1])

    if X.shape[1] == 0:
        print("Saltando porque no hay variables.")
        continue

    X_train, X_test, Y_train, Y_test, lig_train, lig_test = train_test_split(
        X,
        Y,
        df_common["ligando"],
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
    )

    Y_train_np = Y_train.values
    Y_test_np = Y_test.values

    cv = KFold(
        n_splits=min(5, len(X_train)),
        shuffle=True,
        random_state=RANDOM_STATE
    )

    max_pls = max(1, min(8, len(X_train) - 1, X_train.shape[1]))

    pca_options = sorted(set([
        x for x in [3, 5, 8, 10, 15, 20]
        if x <= min(X_train.shape[1], len(X_train) - 1)
    ]))

    if not pca_options:
        pca_options = [1]

    modelos = {
        "Ridge": (
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("scaler", StandardScaler()),
                ("model", Ridge())
            ]),
            {
                "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]
            }
        ),

        "PLSRegression": (
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("scaler", StandardScaler()),
                ("model", PLSRegression())
            ]),
            {
                "model__n_components": list(range(1, max_pls + 1))
            }
        ),

        "PCA_Ridge": (
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("scaler", StandardScaler()),
                ("pca", PCA(random_state=RANDOM_STATE)),
                ("model", Ridge())
            ]),
            {
                "pca__n_components": pca_options,
                "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]
            }
        ),

        "RandomForest": (
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("model", RandomForestRegressor(
                    random_state=RANDOM_STATE,
                    n_jobs=-1
                ))
            ]),
            {
                "model__n_estimators": [300, 600],
                "model__max_features": ["sqrt", 0.5, 1.0],
                "model__min_samples_leaf": [1, 2, 4],
                "model__max_depth": [None, 5, 10]
            }
        ),

        "ExtraTrees": (
            Pipeline([
                ("imputer", SimpleImputer(strategy="median")),
                ("variance", VarianceThreshold()),
                ("model", ExtraTreesRegressor(
                    random_state=RANDOM_STATE,
                    n_jobs=-1
                ))
            ]),
            {
                "model__n_estimators": [300, 600],
                "model__max_features": ["sqrt", 0.5, 1.0],
                "model__min_samples_leaf": [1, 2, 4],
                "model__max_depth": [None, 5, 10]
            }
        ),
    }

    for model_name, (pipeline, param_grid) in modelos.items():

        print(f"\nEntrenando: {set_name} + {model_name}")

        grid = GridSearchCV(
            estimator=pipeline,
            param_grid=param_grid,
            scoring=rmse_scorer,
            cv=cv,
            refit=True,
            n_jobs=1,
            verbose=0
        )

        try:
            grid.fit(X_train, Y_train_np)
        except Exception as e:
            print("ERROR entrenando:", set_name, model_name)
            print(e)
            continue

        pred = grid.predict(X_test)
        pred = np.asarray(pred)

        test_mae = mean_absolute_error(Y_test_np, pred)
        test_rmse = np.sqrt(mean_squared_error(Y_test_np, pred))
        test_r2 = r2_score(Y_test_np, pred, multioutput="uniform_average")

        resumen_rows.append({
            "feature_set": set_name,
            "modelo": model_name,
            "n_ligandos": len(df_common),
            "n_features": X.shape[1],
            "CV_RMSE": -grid.best_score_,
            "TEST_MAE": test_mae,
            "TEST_RMSE": test_rmse,
            "TEST_R2": test_r2,
            "best_params": str(grid.best_params_)
        })

        metricas_rows.extend(
            calcular_metricas(set_name, model_name, Y_test_np, pred)
        )

        pred_df = pd.DataFrame({
            "feature_set": set_name,
            "modelo": model_name,
            "ligando": list(lig_test),
        })

        for i, target in enumerate(TARGETS):
            pred_df[f"real_{target}"] = Y_test_np[:, i]
            pred_df[f"pred_{target}"] = pred[:, i]
            pred_df[f"error_{target}"] = pred[:, i] - Y_test_np[:, i]

        predicciones_rows.append(pred_df)

        # Guardar modelo individual
        safe_name = f"{set_name}_{model_name}".replace("/", "_").replace(" ", "_")
        joblib.dump(grid.best_estimator_, OUTDIR / f"best_{safe_name}.joblib")

        # Guardar resultados CV
        pd.DataFrame(grid.cv_results_).to_csv(
            OUTDIR / f"cv_results_{safe_name}.csv",
            index=False
        )

        print("Best params:", grid.best_params_)
        print(f"CV RMSE: {-grid.best_score_:.4f}")
        print(f"TEST RMSE: {test_rmse:.4f}")
        print(f"TEST MAE: {test_mae:.4f}")
        print(f"TEST R2: {test_r2:.4f}")

        if test_rmse < best_global_rmse:
            best_global_rmse = test_rmse
            best_global_model = grid.best_estimator_
            best_global_name = model_name
            best_global_set = set_name


# =========================
# GUARDAR RESULTADOS
# =========================

resumen = pd.DataFrame(resumen_rows)
metricas = pd.DataFrame(metricas_rows)

if predicciones_rows:
    predicciones = pd.concat(predicciones_rows, ignore_index=True)
else:
    predicciones = pd.DataFrame()

resumen = resumen.sort_values("TEST_RMSE")

resumen.to_csv(OUTDIR / "resumen_comparacion_modelos.csv", index=False)
metricas.to_csv(OUTDIR / "metricas_por_proteina.csv", index=False)
predicciones.to_csv(OUTDIR / "predicciones_test.csv", index=False)

if best_global_model is not None:
    joblib.dump(best_global_model, OUTDIR / "mejor_modelo_global.joblib")

    with open(OUTDIR / "mejor_modelo_global.txt", "w", encoding="utf-8") as f:
        f.write(f"Mejor feature_set: {best_global_set}\n")
        f.write(f"Mejor modelo: {best_global_name}\n")
        f.write(f"TEST_RMSE: {best_global_rmse:.6f}\n")

print("\n" + "=" * 70)
print("FIN")
print("Mejor combinación:")
print("Feature set:", best_global_set)
print("Modelo:", best_global_name)
print("TEST RMSE:", best_global_rmse)

print("\nTop resultados:")
print(resumen.head(20).to_string(index=False))

print("\nArchivos guardados en:")
print(OUTDIR)