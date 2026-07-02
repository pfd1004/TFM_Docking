# Borrador para memoria: resultados de modelado, membrana y nanopartículas

## Metodología resumida del modelado predictivo

A partir de las estructuras optimizadas de los ligandos se calcularon descriptores moleculares 2D/3D mediante RDKit. Las variables respuesta fueron las afinidades de docking obtenidas frente a cuatro proteínas plasmáticas (`y_1ao6`, `y_1hzh`, `y_2hav` y `y_3ghg`). Además, se incorporó un bloque de variables de membrana procedentes del archivo de permeabilidad COSMOtherm, incluyendo `mem_logP`, `mem_logP_kg_l`, `mem_logPerm` y estadísticos derivados de los perfiles de distribución, energía libre, entropía media y difusión a lo largo de la bicapa de DMPC.

Se compararon varios modelos de regresión multisalida: Ridge, PLSRegression, PCA+Ridge, RandomForest y ExtraTrees. La evaluación se realizó mediante partición entrenamiento/test y validación cruzada interna para la selección de hiperparámetros. Las métricas empleadas fueron MAE, RMSE y R² medio sobre las cuatro proteínas.

## Resultados principales del modelo

El mejor rendimiento en test se obtuvo con **ExtraTreesRegressor usando únicamente descriptores RDKit**, con un RMSE medio de **0.412**, MAE de **0.319** y R² medio de **0.884**. El modelo ExtraTrees con RDKit y variables de membrana seleccionadas obtuvo un rendimiento prácticamente equivalente, con RMSE de **0.412**, MAE de **0.339** y R² de **0.886**.

Por tanto, la incorporación de variables de membrana no produjo una mejora clara y consistente del rendimiento predictivo frente al uso exclusivo de descriptores RDKit. Esto sugiere que las variables de membrana aportan información fisicoquímica complementaria sobre el comportamiento del soluto en bicapa lipídica, pero no explican directamente la afinidad proteína-ligando con mayor eficacia que los descriptores moleculares clásicos en este conjunto de datos.

PLSRegression también mostró un comportamiento competitivo con RDKit, con RMSE de **0.426** y R² de **0.886**, lo que es coherente con la utilidad de métodos latentes en conjuntos de datos pequeños y con variables correlacionadas. La reducción dimensional mediante PCA no mejoró el rendimiento, ya que PCA+Ridge obtuvo un RMSE de **0.467** con RDKit.

### Tabla resumen de mejores bloques

| Bloque de variables | Mejor modelo | RMSE test | MAE test | R² test |
|---|---:|---:|---:|---:|
| RDKit solo | ExtraTrees | 0.412 | 0.319 | 0.884 |
| RDKit + membrana | ExtraTrees | 0.412 | 0.339 | 0.886 |
| Membrana solo | RandomForest | 0.731 | 0.621 | 0.650 |

## Interpretación de las variables de membrana

Las variables de membrana por sí solas alcanzaron un rendimiento inferior al de los descriptores RDKit, con RMSE de **0.731** para el mejor modelo de este bloque. Esto indica que contienen señal parcial, pero no suficiente para sustituir a los descriptores estructurales. Químicamente, este resultado es razonable: propiedades como logP, permeabilidad, distribución en membrana o barrera de energía libre describen la partición y movilidad en una bicapa lipídica, mientras que las afinidades objetivo proceden de acoplamiento molecular con proteínas plasmáticas.

## Análisis nanopartícula-proteína

Los resultados de docking nanopartícula-proteína se analizaron como una capa independiente respecto al modelo ligando-proteína. Las energías mínimas obtenidas fueron:

| nanoparticula   |   1ao6 |   1hzh |   2hav |
|:----------------|-------:|-------:|-------:|
| Ag_13           |  -2.24 |  -2.13 |  -2.17 |
| Au_13           |  -2.34 |  -2.12 |  -2.34 |
| ZnO_12_0        |  -7.76 |  -7.38 |  -8.44 |

La nanopartícula **ZnO_12_0** presentó energías más negativas frente a las tres proteínas evaluadas, especialmente en `2hav`, lo que sugiere una interacción directa más favorable en las condiciones de docking empleadas. Estos valores no se introdujeron como variables directas del modelo principal porque la unidad experimental del modelo era el ligando y no se dispone de complejos ligando-nanopartícula-proteína para todos los ligandos.

## Ranking contextual ligando-nanopartícula-proteína

Para integrar de forma exploratoria la información de docking ligando-proteína y nanopartícula-proteína se construyó un ranking contextual. En lugar de sumar energías brutas procedentes de cálculos diferentes, se normalizaron las puntuaciones por proteína y se combinó la afinidad ligando-proteína con la afinidad nanopartícula-proteína. El objetivo fue generar una priorización exploratoria de combinaciones, no un modelo predictivo de complejos funcionalizados.

Top 10 del ranking contextual:

| ligando                                 | nanoparticula   | proteina   |   y_ligando_proteina |   E_nano_proteina |   score_contextual |
|:----------------------------------------|:----------------|:-----------|---------------------:|------------------:|-------------------:|
| Folic_acid                              | ZnO_12_0        | 1ao6       |               -9.21  |             -7.76 |            4.36194 |
| Doxorubicin_proxy_daunorubicin_aglycone | ZnO_12_0        | 1ao6       |               -8.36  |             -7.76 |            3.7936  |
| Curcumin                                | ZnO_12_0        | 1ao6       |               -7.742 |             -7.76 |            3.38038 |
| Methotrexate                            | ZnO_12_0        | 1ao6       |               -7.615 |             -7.76 |            3.29546 |
| Tannic_acid_approx                      | ZnO_12_0        | 1ao6       |               -7.531 |             -7.76 |            3.23929 |
| TAT_peptide_fragment_YGRKKK_6_aa        | ZnO_12_0        | 1ao6       |               -7.4   |             -7.76 |            3.1517  |
| Gemcitabine                             | ZnO_12_0        | 1ao6       |               -6.74  |             -7.76 |            2.7104  |
| Biotin                                  | ZnO_12_0        | 1ao6       |               -6.604 |             -7.76 |            2.61947 |
| RGD_peptide_Arg-Gly-Asp                 | ZnO_12_0        | 1ao6       |               -6.553 |             -7.76 |            2.58537 |
| Folic_acid                              | Au_13           | 1ao6       |               -9.21  |             -2.34 |            2.26027 |

Este análisis debe interpretarse como exploratorio. Para entrenar un modelo específico de complejos ligando-nanopartícula-proteína sería necesario disponer de cálculos completos de docking para múltiples combinaciones de ligandos, nanopartículas y proteínas.

## Frecuencias y archivos Molden

Las frecuencias vibracionales y los archivos Molden no se incluyeron como variables directas en el CSV de descriptores. Su papel en el flujo de trabajo fue validar las geometrías optimizadas y conservar información estructural/orbital. En la memoria puede indicarse que las geometrías se validaron mediante análisis vibracional, comprobando la ausencia de frecuencias imaginarias cuando correspondía. Si se quisieran incorporar en el futuro, deberían extraerse descriptores escalares como número de frecuencias imaginarias, frecuencia mínima, HOMO, LUMO, gap HOMO-LUMO o energía total, pero no introducir el archivo Molden completo como variable.

## Figuras generadas

- `figuras/fig_modelos_top12_rmse.png`: comparación global de modelos.
- `figuras/fig_mejor_por_bloque_rmse.png`: mejor modelo por bloque de variables.
- `figuras/fig_rmse_por_proteina_modelo_final.png`: RMSE por proteína para el modelo final.
- `figuras/fig_real_vs_pred_y_*.png`: real frente a predicho para cada proteína.
- `figuras/fig_nano_proteina_*.png`: docking nanopartícula-proteína.
- `figuras/fig_top15_ranking_contextual.png`: ranking contextual.
