# Paquete de cierre TFM

Este paquete contiene los archivos para cerrar los cuatro bloques:

1. Modelo final y comparación RDKit vs RDKit + membrana.
2. Análisis nanopartícula-proteína y ranking contextual.
3. App Streamlit para visualización.
4. Borrador de texto para memoria.

## Ejecutar la app

En la carpeta principal del proyecto, copia `app_tfm.py` y ejecuta:

```cmd
pip install -r requirements_app.txt
streamlit run app_tfm.py
```

La app espera encontrar estos archivos/rutas:

- `dataset_descriptores_membrana.csv`
- `ranking_contextual_ligando_nano_proteina.csv`
- `tabla_nano_proteina_contexto_wide.csv`
- `resultados_docking_nanoparticulas_corregido.csv`
- `resultados_comparacion_membrana/resumen_comparacion_modelos.csv`
- `resultados_comparacion_membrana/metricas_por_proteina.csv`
- `resultados_comparacion_membrana/predicciones_test.csv`

## Entrenar modelo final

```cmd
python entrenar_modelo_final_rdkit.py
```

Genera `modelo_final_rdkit/modelo_final_extratrees_rdkit.joblib`.

## Resultado seleccionado

Modelo final recomendado: **RDKit_solo + ExtraTreesRegressor**.

- TEST_RMSE: 0.4117
- TEST_MAE: 0.3188
- TEST_R2: 0.8840

Membrana no mejora de forma clara:

- RDKit + membrana TEST_RMSE: 0.4117
- RDKit + membrana TEST_R2: 0.8862

## Archivos importantes

- `borrador_memoria_modelos_nanoparticulas.md`: texto aprovechable para memoria.
- `figuras/`: figuras ya exportadas en PNG.
- `tabla_resultados_modelos_ordenada.csv`: modelos ordenados por RMSE.
- `tabla_top50_ranking_contextual.csv`: top 50 combinaciones ligando-nano-proteína.
