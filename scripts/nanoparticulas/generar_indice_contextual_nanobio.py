# generar_indice_contextual_nanobio.py
# ============================================================
# Renombra y refuerza el ranking contextual como índice propio:
#   ICNB = Índice Contextual Nano-Bio
#
# Entrada:
#   resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina.csv
#
# Salidas:
#   resultados/nanoparticulas/indice_contextual_nanobio.csv
#   resultados/nanoparticulas/top30_indice_contextual_nanobio.csv
#   resultados/nanoparticulas/top10_icnb_por_proteina.csv
#   resultados/nanoparticulas/mejor_combinacion_icnb_por_proteina.csv
#   resultados/figuras_generadas/fig_13_top20_icnb.png
#
# Ejecutar desde raíz:
#   python scripts\nanoparticulas\generar_indice_contextual_nanobio.py
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd()

INPUTS = [
    ROOT / "resultados" / "nanoparticulas" / "ranking_contextual_ligando_nano_proteina.csv",
    ROOT / "resultados" / "nanoparticulas" / "ranking_contextual_ligando_nano_proteina_ordenado.csv",
    ROOT / "ranking_contextual_ligando_nano_proteina.csv",
]

OUT_NANO = ROOT / "resultados" / "nanoparticulas"
OUT_FIG = ROOT / "resultados" / "figuras_generadas"

OUT_NANO.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)


def first_existing(paths):
    for p in paths:
        if p.exists():
            return p
    return None


def minmax_0_100(x):
    x = pd.Series(x, dtype=float)
    mn = x.min(skipna=True)
    mx = x.max(skipna=True)

    if pd.isna(mn) or pd.isna(mx) or mx == mn:
        return pd.Series(np.full(len(x), 50.0), index=x.index)

    return 100.0 * (x - mn) / (mx - mn)


def main():
    path = first_existing(INPUTS)
    if path is None:
        raise FileNotFoundError("No encuentro ranking_contextual_ligando_nano_proteina.csv")

    print("Leyendo:", path)
    df = pd.read_csv(path)

    required = [
        "ligando",
        "nanoparticula",
        "proteina",
        "score_ligando_proteina",
        "score_nano_proteina",
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError("Faltan columnas: " + ", ".join(missing))

    df = df.copy()

    df["score_ligando_proteina"] = pd.to_numeric(df["score_ligando_proteina"], errors="coerce")
    df["score_nano_proteina"] = pd.to_numeric(df["score_nano_proteina"], errors="coerce")

    # Índice contextual bruto: media de las dos puntuaciones normalizadas.
    df["ICNB_bruto"] = (
        df["score_ligando_proteina"] +
        df["score_nano_proteina"]
    ) / 2

    # Escala 0-100 global.
    df["ICNB_global_0_100"] = minmax_0_100(df["ICNB_bruto"])

    # Escala 0-100 dentro de cada proteína. Esta es la más clara para memoria.
    df["ICNB_0_100"] = df.groupby("proteina")["ICNB_bruto"].transform(minmax_0_100)

    df["ranking_ICNB_global"] = df["ICNB_global_0_100"].rank(ascending=False, method="dense").astype(int)
    df["ranking_ICNB_por_proteina"] = df.groupby("proteina")["ICNB_0_100"].rank(ascending=False, method="dense").astype(int)

    # Orden principal: por proteína y luego ICNB.
    df = df.sort_values(
        ["proteina", "ICNB_0_100", "ICNB_global_0_100"],
        ascending=[True, False, False]
    )

    cols_first = [
        "ligando",
        "nanoparticula",
        "proteina",
        "ICNB_0_100",
        "ICNB_global_0_100",
        "ICNB_bruto",
        "ranking_ICNB_por_proteina",
        "ranking_ICNB_global",
        "y_ligando_proteina",
        "E_nano_proteina",
        "score_ligando_proteina",
        "score_nano_proteina",
    ]

    cols = [c for c in cols_first if c in df.columns] + [c for c in df.columns if c not in cols_first]
    df = df[cols]

    out_full = OUT_NANO / "indice_contextual_nanobio.csv"
    df.to_csv(out_full, index=False)

    # Top global
    top30 = df.sort_values("ICNB_global_0_100", ascending=False).head(30)
    top30.to_csv(OUT_NANO / "top30_indice_contextual_nanobio.csv", index=False)

    # Top por proteína
    top10_prot = (
        df.sort_values(["proteina", "ICNB_0_100"], ascending=[True, False])
        .groupby("proteina", as_index=False)
        .head(10)
    )
    top10_prot.to_csv(OUT_NANO / "top10_icnb_por_proteina.csv", index=False)

    best_prot = (
        df.sort_values(["proteina", "ICNB_0_100"], ascending=[True, False])
        .groupby("proteina", as_index=False)
        .head(1)
    )
    best_prot.to_csv(OUT_NANO / "mejor_combinacion_icnb_por_proteina.csv", index=False)

    # Figura top 20 global
    plot_df = top30.head(20).copy()
    plot_df = plot_df.sort_values("ICNB_global_0_100", ascending=True)
    plot_df["label"] = (
        plot_df["ligando"].astype(str).str.slice(0, 24)
        + " | "
        + plot_df["nanoparticula"].astype(str)
        + " | "
        + plot_df["proteina"].astype(str)
    )

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.barh(plot_df["label"], plot_df["ICNB_global_0_100"])
    ax.set_xlabel("ICNB global (0-100)")
    ax.set_title("Top 20 combinaciones por Índice Contextual Nano-Bio")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()

    fig_path = OUT_FIG / "fig_13_top20_icnb.png"
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)

    print("\nCreados:")
    print(out_full)
    print(OUT_NANO / "top30_indice_contextual_nanobio.csv")
    print(OUT_NANO / "top10_icnb_por_proteina.csv")
    print(OUT_NANO / "mejor_combinacion_icnb_por_proteina.csv")
    print(fig_path)

    print("\nTop 10 global:")
    show_cols = ["ligando", "nanoparticula", "proteina", "ICNB_global_0_100", "ICNB_0_100"]
    print(top30[show_cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
