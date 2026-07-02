
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="TFM Docking Nanopartículas", layout="wide")
st.title("TFM · Docking, membrana y nanopartículas")
st.caption("Dashboard sencillo para explorar modelos, membrana, docking nanopartícula-proteína y ranking contextual")

ROOT = Path(".")
PATHS = {
    "resumen": ROOT / "resultados_comparacion_membrana" / "resumen_comparacion_modelos.csv",
    "metricas": ROOT / "resultados_comparacion_membrana" / "metricas_por_proteina.csv",
    "predicciones": ROOT / "resultados_comparacion_membrana" / "predicciones_test.csv",
    "ranking": ROOT / "ranking_contextual_ligando_nano_proteina.csv",
    "nano_wide": ROOT / "tabla_nano_proteina_contexto_wide.csv",
    "nano": ROOT / "resultados_docking_nanoparticulas_corregido.csv",
    "dataset": ROOT / "dataset_descriptores_membrana.csv",
}

def load_csv(path):
    if path.exists():
        return pd.read_csv(path)
    return None

def bar_chart(df, x, y, title, rotation=45):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df[x].astype(str), df[y])
    ax.set_title(title)
    ax.set_ylabel(y)
    ax.tick_params(axis="x", labelrotation=rotation)
    fig.tight_layout()
    return fig

def scatter_chart(real, pred, title):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(real, pred)
    mn = min(np.nanmin(real), np.nanmin(pred))
    mx = max(np.nanmax(real), np.nanmax(pred))
    ax.plot([mn, mx], [mn, mx], linestyle="--")
    ax.set_xlabel("Real")
    ax.set_ylabel("Predicho")
    ax.set_title(title)
    fig.tight_layout()
    return fig

resumen = load_csv(PATHS["resumen"])
metricas = load_csv(PATHS["metricas"])
predicciones = load_csv(PATHS["predicciones"])
ranking = load_csv(PATHS["ranking"])
nano_wide = load_csv(PATHS["nano_wide"])
nano = load_csv(PATHS["nano"])
dataset = load_csv(PATHS["dataset"])

section = st.sidebar.radio("Sección", [
    "Resumen",
    "Modelos",
    "Predicciones",
    "Membrana",
    "Nanopartículas",
    "Ranking contextual",
    "Dataset"
])

if section == "Resumen":
    st.header("Resumen del TFM")
    c1, c2, c3 = st.columns(3)
    if dataset is not None:
        y_cols = [c for c in dataset.columns if c.startswith("y_")]
        mem_cols = [c for c in dataset.columns if c.startswith("mem_")]
        with c1:
            st.metric("Ligandos", dataset.shape[0])
            st.metric("Variables", dataset.shape[1])
        with c2:
            st.metric("Proteínas objetivo", len(y_cols))
            st.metric("Variables membrana", len(mem_cols))
    if resumen is not None:
        best = resumen.sort_values("TEST_RMSE").iloc[0]
        with c3:
            st.metric("Mejor modelo", str(best["modelo"]))
            st.metric("Mejor bloque", str(best["feature_set"]))
            st.metric("TEST_RMSE", f"{best['TEST_RMSE']:.4f}")
        st.subheader("Top 10 modelos")
        st.dataframe(resumen.sort_values("TEST_RMSE").head(10), use_container_width=True)
    st.info("La app es una herramienta auxiliar de visualización; el modelo final y las tablas originales siguen estando en los CSV.")

elif section == "Modelos":
    st.header("Comparación de modelos")
    if resumen is None:
        st.error("No encuentro resultados_comparacion_membrana/resumen_comparacion_modelos.csv")
    else:
        df = resumen.sort_values("TEST_RMSE").copy()
        df["modelo_completo"] = df["feature_set"] + " · " + df["modelo"]
        st.dataframe(df, use_container_width=True)
        topn = st.slider("Número de modelos", 5, min(30, len(df)), 12)
        st.pyplot(bar_chart(df.head(topn), "modelo_completo", "TEST_RMSE", "RMSE en test", 70))
        st.subheader("Mejor modelo por bloque")
        best_by = df.groupby("feature_set", as_index=False).first()
        st.dataframe(best_by, use_container_width=True)
        st.pyplot(bar_chart(best_by, "feature_set", "TEST_RMSE", "Mejor RMSE por bloque", 45))

elif section == "Predicciones":
    st.header("Predicciones real vs predicho")
    if predicciones is None:
        st.error("No encuentro resultados_comparacion_membrana/predicciones_test.csv")
    else:
        df = predicciones.copy()
        c1, c2, c3 = st.columns(3)
        with c1:
            fset = st.selectbox("Bloque", sorted(df["feature_set"].unique()))
        with c2:
            mod = st.selectbox("Modelo", sorted(df["modelo"].unique()))
        with c3:
            target = st.selectbox("Proteína", ["y_1ao6", "y_1hzh", "y_2hav", "y_3ghg"])
        sub = df[(df["feature_set"] == fset) & (df["modelo"] == mod)].copy()
        real_col, pred_col, err_col = f"real_{target}", f"pred_{target}", f"error_{target}"
        st.dataframe(sub[["ligando", real_col, pred_col, err_col]], use_container_width=True)
        st.pyplot(scatter_chart(sub[real_col].values, sub[pred_col].values, f"{target} · {fset} · {mod}"))

elif section == "Membrana":
    st.header("Variables de membrana")
    if dataset is None:
        st.error("No encuentro dataset_descriptores_membrana.csv")
    else:
        mem_cols = [c for c in dataset.columns if c.startswith("mem_")]
        st.write(f"Variables de membrana disponibles: {len(mem_cols)}")
        show_cols = ["ligando"] + mem_cols[:]
        st.dataframe(dataset[show_cols], use_container_width=True)
        numeric_mem = dataset[mem_cols].select_dtypes(include=[np.number])
        if len(numeric_mem.columns) > 0:
            col = st.selectbox("Variable", numeric_mem.columns)
            plot_df = dataset[["ligando", col]].dropna().sort_values(col).tail(20)
            st.pyplot(bar_chart(plot_df, "ligando", col, f"Top 20 {col}", 75))

elif section == "Nanopartículas":
    st.header("Docking nanopartícula-proteína")
    if nano_wide is not None:
        st.subheader("Tabla wide")
        st.dataframe(nano_wide, use_container_width=True)
    if nano is not None:
        st.subheader("Tabla long")
        st.dataframe(nano, use_container_width=True)
        prot = st.selectbox("Proteína", sorted(nano["proteina"].unique()))
        sub = nano[nano["proteina"] == prot].sort_values("binding_energy_min")
        st.pyplot(bar_chart(sub, "nanoparticula", "binding_energy_min", f"Nanopartícula-proteína · {prot}", 0))

elif section == "Ranking contextual":
    st.header("Ranking contextual ligando-nanopartícula-proteína")
    if ranking is None:
        st.error("No encuentro ranking_contextual_ligando_nano_proteina.csv")
    else:
        df = ranking.copy()
        c1, c2 = st.columns(2)
        with c1:
            prot_opt = ["Todas"] + sorted(df["proteina"].dropna().unique())
            prot = st.selectbox("Proteína", prot_opt)
        with c2:
            nano_opt = ["Todas"] + sorted(df["nanoparticula"].dropna().unique())
            nano_sel = st.selectbox("Nanopartícula", nano_opt)
        if prot != "Todas":
            df = df[df["proteina"] == prot]
        if nano_sel != "Todas":
            df = df[df["nanoparticula"] == nano_sel]
        df = df.sort_values("score_contextual", ascending=False)
        st.dataframe(df, use_container_width=True)
        plot_df = df.head(15).copy()
        plot_df["label"] = plot_df["ligando"].astype(str).str.slice(0, 20) + " · " + plot_df["nanoparticula"].astype(str) + " · " + plot_df["proteina"].astype(str)
        st.pyplot(bar_chart(plot_df, "label", "score_contextual", "Top 15 score contextual", 75))

elif section == "Dataset":
    st.header("Dataset completo")
    if dataset is None:
        st.error("No encuentro dataset_descriptores_membrana.csv")
    else:
        st.dataframe(dataset, use_container_width=True)
        st.download_button("Descargar CSV", dataset.to_csv(index=False).encode("utf-8"), "dataset_descriptores_membrana.csv", "text/csv")
