from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import colormaps

try:
    import joblib
except ImportError:
    joblib = None

# ============================================================
# CONFIGURACIÓN
# ============================================================

st.set_page_config(
    page_title="TFM · Docking y nanopartículas",
    layout="wide"
)

st.title("TFM · Docking, membrana, nanopartículas y aprendizaje automático")
st.caption(
    "Dashboard auxiliar para revisar resultados, modelos, variables de membrana, "
    "nanopartículas y ranking contextual."
)

# ============================================================
# RUTAS ROBUSTAS
# La app puede ejecutarse desde:
#   streamlit run app/app_tfm.py
# o incluso desde dentro de app/
# ============================================================

APP_DIR = Path(__file__).resolve().parent
CWD = Path.cwd().resolve()

CANDIDATE_ROOTS = []
for p in [
    CWD,
    CWD.parent,
    APP_DIR,
    APP_DIR.parent,
]:
    p = p.resolve()
    if p not in CANDIDATE_ROOTS:
        CANDIDATE_ROOTS.append(p)


def find_project_root():
    """
    Intenta localizar la raíz del proyecto buscando carpetas esperadas.
    """
    for root in CANDIDATE_ROOTS:
        if (root / "datasets").exists() or (root / "resultados").exists():
            return root

    return CWD


ROOT = find_project_root()


def first_existing(relative_candidates):
    """
    Devuelve el primer archivo existente a partir de una lista de rutas relativas.
    Si ninguno existe, devuelve la primera ruta candidata bajo ROOT para mostrar errores claros.
    """
    for rel in relative_candidates:
        path = ROOT / rel
        if path.exists():
            return path

    return ROOT / relative_candidates[0]


PATHS = {
    # Dataset
    "dataset": first_existing([
        Path("datasets/dataset_descriptores_membrana.csv"),
        Path("dataset_descriptores_membrana.csv"),
    ]),

    # Resultados de modelos
    "resumen": first_existing([
        Path("resultados/modelos/resultados_comparacion_membrana/resumen_comparacion_modelos.csv"),
        Path("resultados/modelos/resumen_comparacion_modelos.csv"),
        Path("resultados/tablas_generadas/tabla_01_resultados_modelos_ordenada.csv"),
        Path("resultados_comparacion_membrana/resumen_comparacion_modelos.csv"),
        Path("resumen_comparacion_modelos.csv"),
    ]),
    "metricas": first_existing([
        Path("resultados/modelos/resultados_comparacion_membrana/metricas_por_proteina.csv"),
        Path("resultados/modelos/metricas_por_proteina.csv"),
        Path("resultados/tablas_generadas/tabla_04_metricas_mejor_modelo_por_proteina.csv"),
        Path("resultados_comparacion_membrana/metricas_por_proteina.csv"),
        Path("metricas_por_proteina.csv"),
    ]),
    "predicciones": first_existing([
        Path("resultados/modelos/resultados_comparacion_membrana/predicciones_test.csv"),
        Path("resultados/modelos/predicciones_test.csv"),
        Path("resultados/tablas_generadas/tabla_06_predicciones_test_mejor_modelo.csv"),
        Path("resultados_comparacion_membrana/predicciones_test.csv"),
        Path("predicciones_test.csv"),
    ]),

    # Diagnósticos
    "correlaciones": first_existing([
        Path("resultados/diagnosticos/correlaciones_features_targets.csv"),
        Path("resultados/tablas_generadas/tabla_10_top30_variables_correlacion_targets.csv"),
        Path("correlaciones_features_targets.csv"),
    ]),

    # Nanopartículas
    "nano": first_existing([
        Path("resultados/nanoparticulas/resultados_docking_nanoparticulas_corregido.csv"),
        Path("resultados_docking_nanoparticulas_corregido.csv"),
    ]),
    "nano_wide": first_existing([
        Path("resultados/nanoparticulas/tabla_nano_proteina_contexto_wide.csv"),
        Path("resultados/nanoparticulas/tabla_nano_proteina_wide.csv"),
        Path("resultados/tablas_generadas/tabla_12_nano_proteina_wide.csv"),
        Path("tabla_nano_proteina_contexto_wide.csv"),
        Path("tabla_nano_proteina_wide.csv"),
    ]),
    "ranking": first_existing([
        Path("resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina_ordenado.csv"),
        Path("resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina.csv"),
        Path("resultados/tablas_generadas/tabla_14_ranking_contextual_completo.csv"),
        Path("ranking_contextual_ligando_nano_proteina.csv"),
    ]),

    # Carpetas de figuras/tablas
    "figuras": ROOT / "resultados" / "figuras_generadas",
    "tablas": ROOT / "resultados" / "tablas_generadas",

    # Modelo final para predicción desde SMILES
    "modelo_final": first_existing([
        Path("resultados/modelos/modelo_final_rdkit/modelo_final_extratrees_rdkit.joblib"),
        Path("modelo_final_rdkit/modelo_final_extratrees_rdkit.joblib"),
    ]),
    "features_modelo_final": first_existing([
        Path("resultados/modelos/modelo_final_rdkit/features_modelo_final.csv"),
        Path("modelo_final_rdkit/features_modelo_final.csv"),
    ]),
}


# ============================================================
# CARGA
# ============================================================

@st.cache_data(show_spinner=False)
def load_csv(path_str):
    path = Path(path_str)
    if path.exists():
        return pd.read_csv(path)
    return None


def load_all():
    return {
        "dataset": load_csv(str(PATHS["dataset"])),
        "resumen": load_csv(str(PATHS["resumen"])),
        "metricas": load_csv(str(PATHS["metricas"])),
        "predicciones": load_csv(str(PATHS["predicciones"])),
        "correlaciones": load_csv(str(PATHS["correlaciones"])),
        "nano": load_csv(str(PATHS["nano"])),
        "nano_wide": load_csv(str(PATHS["nano_wide"])),
        "ranking": load_csv(str(PATHS["ranking"])),
    }


data = load_all()
dataset = data["dataset"]
resumen = data["resumen"]
metricas = data["metricas"]
predicciones = data["predicciones"]
correlaciones = data["correlaciones"]
nano = data["nano"]
nano_wide = data["nano_wide"]
ranking = data["ranking"]


# ============================================================
# FUNCIONES DE PLOT
# ============================================================

PALETTE = {
    "blue": "#1f77b4",
    "orange": "#ff9f1c",
    "green": "#2ca02c",
    "red": "#d62728",
    "purple": "#7b2cbf",
    "gray": "#6c757d",
}


def bar_chart(df, x, y, title, rotation=45, horizontal=False, color=None):
    fig, ax = plt.subplots(figsize=(10, 5))

    if color is None:
        color = PALETTE["blue"]

    if horizontal:
        ax.barh(df[x].astype(str), pd.to_numeric(df[y], errors="coerce"), color=color)
        ax.set_xlabel(y)
        ax.invert_yaxis()
    else:
        ax.bar(df[x].astype(str), pd.to_numeric(df[y], errors="coerce"), color=color)
        ax.set_ylabel(y)
        ax.tick_params(axis="x", labelrotation=rotation)

    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    return fig


def scatter_chart(real, pred, title):
    real = pd.to_numeric(pd.Series(real), errors="coerce")
    pred = pd.to_numeric(pd.Series(pred), errors="coerce")
    ok = real.notna() & pred.notna()

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(real[ok], pred[ok], color=PALETTE["blue"], alpha=0.85)

    if ok.sum() > 0:
        mn = min(real[ok].min(), pred[ok].min())
        mx = max(real[ok].max(), pred[ok].max())
        ax.plot([mn, mx], [mn, mx], linestyle="--", color=PALETTE["red"], linewidth=1.5)

    ax.set_xlabel("Real")
    ax.set_ylabel("Predicho")
    ax.set_title(title)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def heatmap_from_wide(df, title):
    if df is None or df.empty or "nanoparticula" not in df.columns:
        return None

    value_cols = [c for c in df.columns if c != "nanoparticula"]
    if not value_cols:
        return None

    mat = df[value_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    labels_y = df["nanoparticula"].astype(str).tolist()
    labels_x = value_cols

    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(mat, cmap="viridis")

    ax.set_xticks(np.arange(len(labels_x)))
    ax.set_xticklabels(labels_x)
    ax.set_yticks(np.arange(len(labels_y)))
    ax.set_yticklabels(labels_y)

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", color="white", fontsize=9)

    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    return fig


def show_missing_box():
    st.warning("Algunos archivos no se han encontrado. Revisa la pestaña 'Archivos y rutas'.")


def dataframe_download(df, filename):
    if df is not None:
        st.download_button(
            f"Descargar {filename}",
            df.to_csv(index=False).encode("utf-8"),
            filename,
            "text/csv"
        )


def choose_existing_columns(df, preferred):
    return [c for c in preferred if c in df.columns]


# ============================================================
# PREDICCIÓN DESDE SMILES
# ============================================================

@st.cache_resource(show_spinner=False)
def load_final_model(model_path_str):
    if joblib is None:
        return None, "joblib no está instalado."

    path = Path(model_path_str)
    if not path.exists():
        return None, f"No existe el modelo: {path}"

    try:
        return joblib.load(path), ""
    except Exception as e:
        return None, str(e)


@st.cache_data(show_spinner=False)
def load_final_features(features_path_str):
    path = Path(features_path_str)
    if not path.exists():
        return None

    df = pd.read_csv(path)
    if "feature" in df.columns:
        return df["feature"].astype(str).tolist()

    return df.iloc[:, 0].astype(str).tolist()


def calculate_rdkit_features_from_smiles(smiles, feature_cols):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    except ImportError as e:
        raise ImportError("RDKit no está instalado en este entorno.") from e

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES no válido.")

    mol_h = Chem.AddHs(mol)

    # Generar conformación 3D para descriptores geométricos.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol_h, params)

    if status == 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)
            except Exception:
                pass

    features = {}

    # Descriptores 2D/topológicos
    features["MolWt"] = Descriptors.MolWt(mol)
    features["ExactMolWt"] = Descriptors.ExactMolWt(mol)
    features["MolLogP"] = Descriptors.MolLogP(mol)
    features["TPSA"] = rdMolDescriptors.CalcTPSA(mol)
    features["NumHDonors"] = rdMolDescriptors.CalcNumHBD(mol)
    features["NumHAcceptors"] = rdMolDescriptors.CalcNumHBA(mol)
    features["NumRotatableBonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
    features["NumHeavyAtoms"] = mol.GetNumHeavyAtoms()
    features["NumAtoms"] = mol_h.GetNumAtoms()
    features["FractionCSP3"] = rdMolDescriptors.CalcFractionCSP3(mol)
    features["NumRings"] = rdMolDescriptors.CalcNumRings(mol)
    features["FormalCharge"] = Chem.GetFormalCharge(mol)
    features["Chi0n"] = rdMolDescriptors.CalcChi0n(mol)
    features["Chi1n"] = rdMolDescriptors.CalcChi1n(mol)
    features["Kappa1"] = rdMolDescriptors.CalcKappa1(mol)
    features["Kappa2"] = rdMolDescriptors.CalcKappa2(mol)
    features["LabuteASA"] = rdMolDescriptors.CalcLabuteASA(mol_h)

    # Descriptores 3D. Si falló la conformación, quedan NaN y el pipeline imputará.
    three_d_funcs = {
        "Asphericity": rdMolDescriptors.CalcAsphericity,
        "Eccentricity": rdMolDescriptors.CalcEccentricity,
        "InertialShapeFactor": rdMolDescriptors.CalcInertialShapeFactor,
        "NPR1": rdMolDescriptors.CalcNPR1,
        "NPR2": rdMolDescriptors.CalcNPR2,
        "PMI1": rdMolDescriptors.CalcPMI1,
        "PMI2": rdMolDescriptors.CalcPMI2,
        "PMI3": rdMolDescriptors.CalcPMI3,
        "RadiusOfGyration": rdMolDescriptors.CalcRadiusOfGyration,
        "SpherocityIndex": rdMolDescriptors.CalcSpherocityIndex,
    }

    for name, func in three_d_funcs.items():
        try:
            features[name] = func(mol_h)
        except Exception:
            features[name] = np.nan

    X = pd.DataFrame([{c: features.get(c, np.nan) for c in feature_cols}])
    return X, features


# ============================================================
# SIDEBAR
# ============================================================

with st.sidebar:
    st.header("Navegación")
    section = st.radio(
        "Sección",
        [
            "Resumen",
            "Modelos",
            "Métricas por proteína",
            "Predicciones",
            "Membrana",
            "Nanopartículas",
            "Ranking contextual",
            "Predicción nuevo ligando",
            "Figuras generadas",
            "Dataset",
            "Archivos y rutas",
        ]
    )

    st.divider()
    st.caption("Raíz detectada")
    st.code(str(ROOT), language="text")


# ============================================================
# SECCIONES
# ============================================================

if section == "Resumen":
    st.header("Resumen general")

    c1, c2, c3, c4 = st.columns(4)

    if dataset is not None:
        y_cols = [c for c in dataset.columns if c.startswith("y_")]
        mem_cols = [c for c in dataset.columns if c.startswith("mem_")]

        with c1:
            st.metric("Ligandos", dataset.shape[0])
        with c2:
            st.metric("Variables", dataset.shape[1])
        with c3:
            st.metric("Proteínas objetivo", len(y_cols))
        with c4:
            st.metric("Variables membrana", len(mem_cols))
    else:
        show_missing_box()

    if resumen is not None and "TEST_RMSE" in resumen.columns:
        st.subheader("Mejor modelo")
        df_res = resumen.copy()
        df_res["TEST_RMSE"] = pd.to_numeric(df_res["TEST_RMSE"], errors="coerce")
        best = df_res.sort_values("TEST_RMSE").iloc[0]

        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Modelo", str(best.get("modelo", "")))
        with c2:
            st.metric("Bloque", str(best.get("feature_set", "")))
        with c3:
            st.metric("TEST_RMSE", f"{best.get('TEST_RMSE', np.nan):.4f}")
        with c4:
            if "TEST_R2" in best.index:
                st.metric("TEST_R²", f"{best.get('TEST_R2', np.nan):.4f}")

        st.dataframe(df_res.sort_values("TEST_RMSE").head(10), use_container_width=True)
    else:
        show_missing_box()

    st.info(
        "Esta app solo visualiza resultados ya calculados. "
        "Las tablas base siguen siendo los CSV guardados en `datasets/` y `resultados/`."
    )


elif section == "Modelos":
    st.header("Comparación de modelos")

    if resumen is None:
        st.error(f"No encuentro el resumen de modelos: {PATHS['resumen']}")
    else:
        df = resumen.copy()

        if "TEST_RMSE" in df.columns:
            df["TEST_RMSE"] = pd.to_numeric(df["TEST_RMSE"], errors="coerce")
            df = df.sort_values("TEST_RMSE")

        if all(c in df.columns for c in ["feature_set", "modelo"]):
            df["modelo_completo"] = df["feature_set"].astype(str) + " · " + df["modelo"].astype(str)

        st.dataframe(df, use_container_width=True)
        dataframe_download(df, "resumen_modelos.csv")

        if "TEST_RMSE" in df.columns and "modelo_completo" in df.columns:
            topn = st.slider("Número de modelos a mostrar", 5, min(30, len(df)), min(12, len(df)))
            plot_df = df.head(topn).copy()
            st.pyplot(
                bar_chart(
                    plot_df,
                    "modelo_completo",
                    "TEST_RMSE",
                    "Top modelos por RMSE en test",
                    horizontal=True,
                    color=PALETTE["blue"]
                )
            )

        if all(c in df.columns for c in ["feature_set", "TEST_RMSE"]):
            st.subheader("Mejor modelo por bloque de variables")
            best_by = df.sort_values("TEST_RMSE").groupby("feature_set", as_index=False).first()
            st.dataframe(best_by, use_container_width=True)

            st.pyplot(
                bar_chart(
                    best_by,
                    "feature_set",
                    "TEST_RMSE",
                    "Mejor RMSE por bloque de variables",
                    rotation=35,
                    color=PALETTE["green"]
                )
            )


elif section == "Métricas por proteína":
    st.header("Métricas por proteína")

    if metricas is None:
        st.error(f"No encuentro métricas por proteína: {PATHS['metricas']}")
    else:
        df = metricas.copy()

        if resumen is not None and all(c in resumen.columns for c in ["feature_set", "modelo", "TEST_RMSE"]):
            res = resumen.copy()
            res["TEST_RMSE"] = pd.to_numeric(res["TEST_RMSE"], errors="coerce")
            best = res.sort_values("TEST_RMSE").iloc[0]
            best_fset = best["feature_set"]
            best_model = best["modelo"]

            if all(c in df.columns for c in ["feature_set", "modelo"]):
                only_best = st.checkbox(
                    f"Mostrar solo mejor modelo: {best_fset} + {best_model}",
                    value=True
                )
                if only_best:
                    df = df[(df["feature_set"] == best_fset) & (df["modelo"] == best_model)].copy()

        st.dataframe(df, use_container_width=True)
        dataframe_download(df, "metricas_por_proteina.csv")

        if all(c in df.columns for c in ["salida", "RMSE"]):
            sub = df[df["salida"].astype(str).str.lower() != "media"].copy()
            st.pyplot(
                bar_chart(
                    sub,
                    "salida",
                    "RMSE",
                    "RMSE por proteína",
                    color=PALETTE["orange"]
                )
            )

        if all(c in df.columns for c in ["salida", "R2"]):
            sub = df[df["salida"].astype(str).str.lower() != "media"].copy()
            st.pyplot(
                bar_chart(
                    sub,
                    "salida",
                    "R2",
                    "R² por proteína",
                    color=PALETTE["green"]
                )
            )


elif section == "Predicciones":
    st.header("Predicciones real vs predicho")

    if predicciones is None:
        st.error(f"No encuentro predicciones: {PATHS['predicciones']}")
    else:
        df = predicciones.copy()

        if "feature_set" in df.columns and "modelo" in df.columns:
            c1, c2, c3 = st.columns(3)
            with c1:
                fset = st.selectbox("Bloque", sorted(df["feature_set"].dropna().unique()))
            with c2:
                modelos_disponibles = sorted(df[df["feature_set"] == fset]["modelo"].dropna().unique())
                mod = st.selectbox("Modelo", modelos_disponibles)
            with c3:
                targets = [c.replace("real_", "") for c in df.columns if c.startswith("real_y_")]
                if not targets:
                    targets = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]
                target = st.selectbox("Proteína", targets)

            sub = df[(df["feature_set"] == fset) & (df["modelo"] == mod)].copy()
        else:
            c1, c2 = st.columns(2)
            with c1:
                targets = [c.replace("real_", "") for c in df.columns if c.startswith("real_y_")]
                if not targets:
                    targets = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]
                target = st.selectbox("Proteína", targets)
            sub = df.copy()
            fset = ""
            mod = ""

        real_col = f"real_{target}"
        pred_col = f"pred_{target}"
        err_col = f"error_{target}"

        cols = choose_existing_columns(sub, ["ligando", real_col, pred_col, err_col])
        if cols:
            st.dataframe(sub[cols], use_container_width=True)
            dataframe_download(sub[cols], f"predicciones_{target}.csv")

        if real_col in sub.columns and pred_col in sub.columns:
            st.pyplot(
                scatter_chart(
                    sub[real_col].values,
                    sub[pred_col].values,
                    f"{target} · {fset} · {mod}"
                )
            )
        else:
            st.warning("No encuentro columnas real_/pred_ para esa proteína.")


elif section == "Membrana":
    st.header("Variables de membrana")

    if dataset is None:
        st.error(f"No encuentro dataset: {PATHS['dataset']}")
    else:
        mem_cols = [c for c in dataset.columns if c.startswith("mem_")]
        st.write(f"Variables de membrana disponibles: **{len(mem_cols)}**")

        if mem_cols:
            show_cols = ["ligando"] + mem_cols if "ligando" in dataset.columns else mem_cols
            st.dataframe(dataset[show_cols], use_container_width=True)
            dataframe_download(dataset[show_cols], "variables_membrana.csv")

            numeric_mem = dataset[mem_cols].apply(pd.to_numeric, errors="coerce")
            numeric_cols = [c for c in numeric_mem.columns if numeric_mem[c].notna().sum() > 0]

            if numeric_cols:
                col = st.selectbox("Variable de membrana", numeric_cols)
                plot_df = dataset[["ligando", col]].copy() if "ligando" in dataset.columns else dataset[[col]].copy()
                plot_df[col] = pd.to_numeric(plot_df[col], errors="coerce")
                plot_df = plot_df.dropna().sort_values(col)

                if "ligando" not in plot_df.columns:
                    plot_df["ligando"] = plot_df.index.astype(str)

                st.pyplot(
                    bar_chart(
                        plot_df.tail(20),
                        "ligando",
                        col,
                        f"Top 20 valores de {col}",
                        rotation=70,
                        color=PALETTE["purple"]
                    )
                )

            if all(c in dataset.columns for c in ["mem_logP", "mem_logPerm"]):
                fig, ax = plt.subplots(figsize=(7, 5))
                x = pd.to_numeric(dataset["mem_logP"], errors="coerce")
                y = pd.to_numeric(dataset["mem_logPerm"], errors="coerce")
                ok = x.notna() & y.notna()
                ax.scatter(x[ok], y[ok], color=PALETTE["blue"], alpha=0.85)
                ax.set_xlabel("mem_logP")
                ax.set_ylabel("mem_logPerm")
                ax.set_title("mem_logP vs mem_logPerm")
                ax.grid(alpha=0.25)
                fig.tight_layout()
                st.pyplot(fig)
        else:
            st.warning("El dataset no contiene columnas que empiecen por `mem_`.")


elif section == "Nanopartículas":
    st.header("Docking nanopartícula-proteína")

    if nano_wide is not None:
        st.subheader("Tabla nanopartícula-proteína")
        st.dataframe(nano_wide, use_container_width=True)
        dataframe_download(nano_wide, "tabla_nano_proteina_wide.csv")

        fig = heatmap_from_wide(nano_wide, "Energía docking nanopartícula-proteína")
        if fig is not None:
            st.pyplot(fig)
    else:
        st.warning(f"No encuentro tabla wide: {PATHS['nano_wide']}")

    if nano is not None:
        st.subheader("Resultados long")
        st.dataframe(nano, use_container_width=True)
        dataframe_download(nano, "resultados_docking_nanoparticulas.csv")

        if all(c in nano.columns for c in ["proteina", "nanoparticula"]):
            energy_col = None
            for candidate in ["binding_energy_min", "E_nano_proteina", "energia", "energy"]:
                if candidate in nano.columns:
                    energy_col = candidate
                    break

            if energy_col is not None:
                prot = st.selectbox("Proteína", sorted(nano["proteina"].dropna().unique()))
                sub = nano[nano["proteina"] == prot].copy()
                sub[energy_col] = pd.to_numeric(sub[energy_col], errors="coerce")
                sub = sub.sort_values(energy_col)

                st.pyplot(
                    bar_chart(
                        sub,
                        "nanoparticula",
                        energy_col,
                        f"Docking nanopartícula-proteína · {prot}",
                        rotation=0,
                        color=PALETTE["orange"]
                    )
                )
    else:
        st.warning(f"No encuentro tabla long: {PATHS['nano']}")


elif section == "Ranking contextual":
    st.header("Ranking contextual ligando-nanopartícula-proteína")

    if ranking is None:
        st.error(f"No encuentro ranking contextual: {PATHS['ranking']}")
    else:
        df = ranking.copy()

        c1, c2, c3 = st.columns(3)

        with c1:
            if "proteina" in df.columns:
                prot_opt = ["Todas"] + sorted(df["proteina"].dropna().unique())
                prot = st.selectbox("Proteína", prot_opt)
            else:
                prot = "Todas"

        with c2:
            if "nanoparticula" in df.columns:
                nano_opt = ["Todas"] + sorted(df["nanoparticula"].dropna().unique())
                nano_sel = st.selectbox("Nanopartícula", nano_opt)
            else:
                nano_sel = "Todas"

        with c3:
            topn = st.slider("Top N", 5, 100, 30)

        if prot != "Todas":
            df = df[df["proteina"] == prot]
        if nano_sel != "Todas":
            df = df[df["nanoparticula"] == nano_sel]

        if "score_contextual" in df.columns:
            df["score_contextual"] = pd.to_numeric(df["score_contextual"], errors="coerce")
            df = df.sort_values("score_contextual", ascending=False)

        st.dataframe(df.head(topn), use_container_width=True)
        dataframe_download(df, "ranking_contextual_filtrado.csv")

        if all(c in df.columns for c in ["ligando", "nanoparticula", "proteina", "score_contextual"]):
            plot_df = df.head(min(15, len(df))).copy()
            plot_df["label"] = (
                plot_df["ligando"].astype(str).str.slice(0, 18)
                + " · "
                + plot_df["nanoparticula"].astype(str)
                + " · "
                + plot_df["proteina"].astype(str)
            )

            st.pyplot(
                bar_chart(
                    plot_df,
                    "label",
                    "score_contextual",
                    "Top score contextual",
                    horizontal=True,
                    color=PALETTE["green"]
                )
            )



elif section == "Predicción nuevo ligando":
    st.header("Predicción de afinidad desde SMILES")

    st.write(
        "Esta pestaña convierte un SMILES en descriptores RDKit, carga el modelo final "
        "ExtraTrees y predice las afinidades frente a las cuatro proteínas. "
        "Es una extensión práctica de la app, útil como cribado preliminar."
    )

    if joblib is None:
        st.error("Falta `joblib`. Instala con: `pip install joblib`.")
    else:
        model, model_error = load_final_model(str(PATHS["modelo_final"]))
        feature_cols = load_final_features(str(PATHS["features_modelo_final"]))

        c1, c2 = st.columns(2)
        with c1:
            st.caption("Modelo detectado")
            st.code(str(PATHS["modelo_final"]), language="text")
        with c2:
            st.caption("Features detectadas")
            st.code(str(PATHS["features_modelo_final"]), language="text")

        if model is None:
            st.error(model_error)
            st.info("Primero ejecuta `python scripts\\modelos\\entrenar_modelo_final_rdkit.py`.")
        elif feature_cols is None:
            st.error("No encuentro `features_modelo_final.csv`.")
        else:
            smiles = st.text_input(
                "Introduce un SMILES",
                value="O=C(O)c1ccc(N)cc1",
                help="Ejemplo: ácido p-aminobenzoico."
            )

            nombre = st.text_input("Nombre opcional", value="nuevo_ligando")

            if st.button("Predecir afinidades"):
                try:
                    X, raw_features = calculate_rdkit_features_from_smiles(smiles, feature_cols)
                    pred = model.predict(X)[0]

                    targets = ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"]

                    pred_df = pd.DataFrame({
                        "proteina": targets,
                        "afinidad_predicha": pred,
                    })

                    pred_df["afinidad_predicha"] = pd.to_numeric(
                        pred_df["afinidad_predicha"],
                        errors="coerce"
                    ).round(4)

                    st.success("Predicción completada.")
                    st.dataframe(pred_df, use_container_width=True)

                    fig, ax = plt.subplots(figsize=(7, 4))
                    ax.bar(pred_df["proteina"], pred_df["afinidad_predicha"], color=PALETTE["blue"])
                    ax.set_ylabel("Afinidad predicha")
                    ax.set_title(f"Perfil predicho · {nombre}")
                    ax.grid(axis="y", alpha=0.25)
                    fig.tight_layout()
                    st.pyplot(fig)

                    st.subheader("Descriptores calculados")
                    X_show = X.T.reset_index()
                    X_show.columns = ["descriptor", "valor"]
                    st.dataframe(X_show, use_container_width=True)

                    out = pred_df.copy()
                    out.insert(0, "ligando", nombre)
                    out.insert(1, "SMILES", smiles)

                    st.download_button(
                        "Descargar predicción CSV",
                        out.to_csv(index=False).encode("utf-8"),
                        f"prediccion_{nombre}.csv",
                        "text/csv"
                    )

                    st.warning(
                        "Estas predicciones son orientativas. El modelo se entrenó con el espacio químico "
                        "del dataset del TFM; moléculas muy distintas pueden quedar fuera de su dominio de aplicabilidad."
                    )

                except Exception as e:
                    st.error(str(e))

elif section == "Figuras generadas":
    st.header("Figuras generadas para la memoria")

    figs_dir = PATHS["figuras"]

    if not figs_dir.exists():
        st.warning(
            f"No existe la carpeta {figs_dir}. "
            "Ejecuta antes `Rscript scripts\\utilidades\\generar_figuras_color_tfm.R`."
        )
    else:
        images = sorted(list(figs_dir.glob("*.png")))

        if not images:
            st.warning(f"No hay PNG en {figs_dir}.")
        else:
            st.write(f"Figuras encontradas: **{len(images)}**")

            selected = st.selectbox("Figura", images, format_func=lambda p: p.name)

            st.image(str(selected), use_container_width=True)

            with open(selected, "rb") as f:
                st.download_button(
                    "Descargar figura",
                    f.read(),
                    file_name=selected.name,
                    mime="image/png"
                )

            st.subheader("Galería")
            cols = st.columns(2)
            for i, img in enumerate(images):
                with cols[i % 2]:
                    st.image(str(img), caption=img.name, use_container_width=True)


elif section == "Dataset":
    st.header("Dataset completo")

    if dataset is None:
        st.error(f"No encuentro dataset: {PATHS['dataset']}")
    else:
        st.dataframe(dataset, use_container_width=True)
        dataframe_download(dataset, "dataset_descriptores_membrana.csv")

        numeric_cols = dataset.select_dtypes(include=[np.number]).columns.tolist()
        if numeric_cols:
            st.subheader("Resumen numérico")
            st.dataframe(dataset[numeric_cols].describe().T, use_container_width=True)


elif section == "Archivos y rutas":
    st.header("Archivos y rutas detectadas")

    st.write("Raíz detectada:")
    st.code(str(ROOT), language="text")

    rows = []
    for name, path in PATHS.items():
        exists = path.exists() if isinstance(path, Path) else False
        rows.append({
            "clave": name,
            "existe": exists,
            "ruta": str(path),
        })

    status = pd.DataFrame(rows)
    st.dataframe(status, use_container_width=True)

    missing = status[~status["existe"]]
    if len(missing) > 0:
        st.warning("Hay rutas que no existen. Si todavía no has regenerado CSV/figuras, es normal.")
    else:
        st.success("Todas las rutas principales existen.")
