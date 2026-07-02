Estructura del proyecto TFM
===========================

Ejecuta siempre los scripts desde la raiz TFM_DOCKING.

Carpetas:
- datasets: CSV de entrada.
- scripts/modelos: entrenamiento y comparacion de modelos.
- scripts/membrana: scripts de membrana y COSMOtherm.
- scripts/nanoparticulas: docking de nanoparticulas y ranking contextual.
- scripts/diagnostico: diagnostico de columnas y correlaciones.
- resultados/modelos: resultados de ML.
- resultados/nanoparticulas: tablas nano-proteina y ranking.
- resultados/diagnosticos: diagnosticos.
- app: app Streamlit.

Comandos:
python scripts\modelos\entrenar_comparativo_membrana.py
python scripts\modelos\entrenar_modelo_final.py
streamlit run app\app_tfm.py

Los .xlsx antiguos se han movido a resultados/xlsx_antiguos.
Los scripts han sido modificados para no generar nuevos .xlsx.
