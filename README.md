# TFM_Docking

Repositorio del Trabajo de Fin de Máster:

**Modelización computacional de la interacción de agentes de recubrimiento de nanopartículas metálicas (Au, Ag, ZnO) con proteínas plasmáticas mediante acoplamiento molecular y aprendizaje automático**

**Autor:** Pablo Fuentes de Mateo  
**Máster:** Ingeniería Biomédica, Universidad de Burgos  
**Tutores:** Santiago Aparicio Martínez y Pedro A. Marcos Villa

## Descripción

Este repositorio contiene el flujo computacional desarrollado para estudiar la interacción entre agentes de recubrimiento de nanopartículas metálicas y proteínas plasmáticas mediante docking molecular, descriptores moleculares, variables de membrana y aprendizaje automático.

El trabajo integra:

- docking ligando-proteína frente a proteínas plasmáticas;
- descriptores moleculares RDKit 2D/3D;
- variables de membrana COSMOtherm/COSMOmic;
- modelos de regresión multisalida;
- docking nanopartícula-proteína;
- Índice Contextual Nano-Bio (ICNB);
- análisis estructural preliminar mediante PLIP;
- aplicación Streamlit para visualización y predicción desde SMILES.

## Estructura del repositorio

```text
TFM_Docking/
├── app/                         # Aplicación Streamlit
├── datasets/                    # Datasets de entrada e integrados
├── scripts/
│   ├── modelos/                 # Entrenamiento y comparación de modelos
│   ├── membrana/                # Procesamiento de variables de membrana
│   ├── nanoparticulas/          # Docking NP-proteína e ICNB
│   ├── interacciones/           # Preparación/análisis PLIP
│   └── tablas/                  # Generación de tablas y figuras
├── resultados/
│   ├── modelos/
│   ├── nanoparticulas/
│   ├── diagnosticos/
│   ├── tablas_generadas/
│   ├── figuras_generadas/
│   └── plip/
└── PlantillaTFM_Quarto_TFM/     # Memoria y anexos en Quarto
```

## Sistemas estudiados

Proteínas plasmáticas utilizadas como salidas del modelo:

| Código PDB | Proteína | Variable |
|---|---|---|
| 1AO6 | Albúmina sérica humana | `y_1ao6` |
| 1HZH | Anticuerpo monoclonal humano IgG1 b12 | `y_1hzh` |
| 2HAV | Transferrina sérica humana en forma apo | `y_2hav` |
| 3GHG | Fibrinógeno humano | `y_3ghg` |

Modelos simplificados de nanopartículas:

| Modelo interno | Notación |
|---|---|
| `Ag_13` | (Ag)<sub>13</sub> |
| `Au_13` | (Au)<sub>13</sub> |
| `ZnO_12_0` | (ZnO)<sub>12</sub> |

## Software principal

Versiones registradas durante el trabajo:

| Software | Versión |
|---|---|
| Python | 3.10.20 |
| scikit-learn | 1.7.2 |
| Open Babel | 3.1.0 |
| PLIP | 3.0.0 |
| TURBOMOLE | v7.8.1 |
| ORCA | 5.0 |
| AutoDock Vina | 1.2.0 |
| AutoDock | 4.2.6 |
| MGLTools | 2011 |
| RDKit | 2026.03.3 |

También se utilizaron R, Quarto, LaTeX, Streamlit y COSMOtherm/COSMOmic.

## Instalación básica

El trabajo se desarrolló principalmente en un entorno conda llamado `docking`.

```bash
conda create -n docking python=3.10
conda activate docking
pip install pandas numpy matplotlib scikit-learn joblib streamlit
conda install -c conda-forge rdkit
```

Para la aplicación Streamlit:

```bash
pip install -r app/requirements_app.txt
```

## Ejecución de la aplicación

Desde la raíz del repositorio:

```bash
conda activate docking
streamlit run app/app.py
```

La aplicación permite consultar resultados, comparar modelos, revisar métricas, explorar el ICNB y predecir afinidades de nuevos ligandos a partir de SMILES.

## Modelo predictivo

El modelo final seleccionado fue:

```text
RDKit_solo + ExtraTreesRegressor
```

Archivos esperados por la aplicación:

```text
resultados/modelos/modelo_final_rdkit/modelo_final_extratrees_rdkit.joblib
resultados/modelos/modelo_final_rdkit/features_modelo_final.csv
```

## Generación de resultados

Los scripts principales se encuentran en `scripts/`. Las tablas y figuras finales se generan desde `scripts/tablas/` y se guardan en:

```text
resultados/tablas_generadas/
resultados/figuras_generadas/
```

## Resultados principales

- El mejor rendimiento predictivo se obtuvo con `RDKit_solo + ExtraTreesRegressor`.
- Las variables de membrana aportaron información fisicoquímica complementaria, pero no mejoraron claramente el error frente al bloque RDKit.
- El modelo (ZnO)<sub>12</sub> mostró energías nanopartícula-proteína más favorables que (Ag)<sub>13</sub> y (Au)<sub>13</sub> dentro del protocolo utilizado.
- El ácido fólico destacó como ligando mejor posicionado en el análisis contextual.
- El ICNB se utilizó como métrica comparativa de priorización, no como energía física absoluta.

## Cita

```text
Fuentes de Mateo, P. Modelización computacional de la interacción de agentes de recubrimiento de nanopartículas metálicas (Au, Ag, ZnO) con proteínas plasmáticas mediante acoplamiento molecular y aprendizaje automático. Trabajo de Fin de Máster, Universidad de Burgos, 2026.
```
