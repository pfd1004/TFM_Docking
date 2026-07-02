
# -------------------------------------------------------------------------
# Utilidades generales
# -------------------------------------------------------------------------

crear_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

buscar_archivo <- function(candidatos, obligatorio = TRUE) {
  for (x in candidatos) {
    if (file.exists(x)) {
      return(x)
    }
  }

  if (obligatorio) {
    stop(
      paste0(
        "No encuentro ninguno de estos archivos:\n",
        paste0(" - ", candidatos, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  return(NA_character_)
}

leer_csv <- function(candidatos, obligatorio = TRUE) {
  path <- buscar_archivo(candidatos, obligatorio = obligatorio)

  if (is.na(path)) {
    return(NULL)
  }

  cat("Leyendo:", path, "\n")
  read.csv(path, check.names = FALSE)
}

guardar_csv <- function(df, path) {
  crear_dir(dirname(path))
  write.csv(df, path, row.names = FALSE, na = "")
  cat("Creado:", path, "\n")
}

redondear_columnas <- function(df, cols, digits = 4) {
  for (c in cols) {
    if (c %in% names(df)) {
      df[[c]] <- round(as.numeric(df[[c]]), digits)
    }
  }
  df
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - m) / s
}

limpiar_nombre_archivo <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

abrir_png <- function(path, width = 2200, height = 1600, res = 300) {
  crear_dir(dirname(path))
  png(filename = path, width = width, height = height, res = res)
  par(mar = c(9, 5, 4, 2) + 0.1)
}

cerrar_png <- function(path) {
  dev.off()
  cat("Creada figura:", path, "\n")
}

# -------------------------------------------------------------------------
# Rutas de salida
# -------------------------------------------------------------------------

OUT_TABLAS <- "resultados/tablas_generadas"
OUT_FIGS <- "resultados/figuras_generadas"
OUT_MODELOS <- "resultados/modelos"
OUT_NANO <- "resultados/nanoparticulas"
OUT_DIAG <- "resultados/diagnosticos"

crear_dir(OUT_TABLAS)
crear_dir(OUT_FIGS)
crear_dir(OUT_MODELOS)
crear_dir(OUT_NANO)
crear_dir(OUT_DIAG)

# -------------------------------------------------------------------------
# 1. Cargar archivos principales
# -------------------------------------------------------------------------

resumen <- leer_csv(c(
  "resultados/modelos/resultados_comparacion_membrana/resumen_comparacion_modelos.csv",
  "resultados/modelos/resumen_comparacion_modelos.csv",
  "resultados_comparacion_membrana/resumen_comparacion_modelos.csv",
  "resumen_comparacion_modelos.csv"
), obligatorio = FALSE)

metricas <- leer_csv(c(
  "resultados/modelos/resultados_comparacion_membrana/metricas_por_proteina.csv",
  "resultados/modelos/metricas_por_proteina.csv",
  "resultados_comparacion_membrana/metricas_por_proteina.csv",
  "metricas_por_proteina.csv"
), obligatorio = FALSE)

predicciones <- leer_csv(c(
  "resultados/modelos/resultados_comparacion_membrana/predicciones_test.csv",
  "resultados/modelos/predicciones_test.csv",
  "resultados_comparacion_membrana/predicciones_test.csv",
  "predicciones_test.csv"
), obligatorio = FALSE)

dataset <- leer_csv(c(
  "datasets/dataset_descriptores_membrana.csv",
  "dataset_descriptores_membrana.csv"
), obligatorio = FALSE)

nano <- leer_csv(c(
  "resultados/nanoparticulas/resultados_docking_nanoparticulas_corregido.csv",
  "resultados_docking_nanoparticulas_corregido.csv"
), obligatorio = FALSE)

ranking <- leer_csv(c(
  "resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina.csv",
  "ranking_contextual_ligando_nano_proteina.csv"
), obligatorio = FALSE)

correlaciones <- leer_csv(c(
  "resultados/diagnosticos/correlaciones_features_targets.csv",
  "correlaciones_features_targets.csv"
), obligatorio = FALSE)

# Objetos que se van creando
resumen_ordenado <- NULL
mejor_global <- NULL
mejor_feature_set <- NULL
mejor_modelo <- NULL
metricas_mejor <- NULL
pred_mejor <- NULL
tabla_nano_long <- NULL
nano_wide <- NULL
ranking_limpio <- NULL

# -------------------------------------------------------------------------
# 2. Tablas de modelos
# -------------------------------------------------------------------------

if (!is.null(resumen)) {
  required <- c("feature_set", "modelo", "TEST_RMSE")
  if (!all(required %in% names(resumen))) {
    stop("El resumen de modelos no tiene las columnas mínimas: feature_set, modelo, TEST_RMSE", call. = FALSE)
  }

  resumen$TEST_RMSE <- as.numeric(resumen$TEST_RMSE)
  if ("TEST_MAE" %in% names(resumen)) resumen$TEST_MAE <- as.numeric(resumen$TEST_MAE)
  if ("TEST_R2" %in% names(resumen)) resumen$TEST_R2 <- as.numeric(resumen$TEST_R2)
  if ("CV_RMSE" %in% names(resumen)) resumen$CV_RMSE <- as.numeric(resumen$CV_RMSE)

  resumen_ordenado <- resumen[order(resumen$TEST_RMSE), ]
  resumen_ordenado <- redondear_columnas(
    resumen_ordenado,
    c("CV_RMSE", "TEST_MAE", "TEST_RMSE", "TEST_R2"),
    digits = 4
  )

  guardar_csv(
    resumen_ordenado,
    file.path(OUT_MODELOS, "tabla_resultados_modelos_ordenada.csv")
  )

  guardar_csv(
    resumen_ordenado,
    file.path(OUT_TABLAS, "tabla_01_resultados_modelos_ordenada.csv")
  )

  # Mejor modelo por bloque de variables
  bloques <- unique(resumen_ordenado$feature_set)
  mejor_por_bloque <- do.call(
    rbind,
    lapply(bloques, function(b) {
      sub <- resumen_ordenado[resumen_ordenado$feature_set == b, ]
      sub[which.min(sub$TEST_RMSE), ]
    })
  )

  mejor_por_bloque <- mejor_por_bloque[order(mejor_por_bloque$TEST_RMSE), ]

  guardar_csv(
    mejor_por_bloque,
    file.path(OUT_MODELOS, "tabla_mejor_modelo_por_bloque.csv")
  )

  guardar_csv(
    mejor_por_bloque,
    file.path(OUT_TABLAS, "tabla_02_mejor_modelo_por_bloque.csv")
  )

  # Mejor modelo global
  mejor_global <- resumen_ordenado[1, , drop = FALSE]
  mejor_feature_set <- mejor_global$feature_set[1]
  mejor_modelo <- mejor_global$modelo[1]

  guardar_csv(
    mejor_global,
    file.path(OUT_TABLAS, "tabla_03_mejor_modelo_global.csv")
  )

  cat("\nMejor modelo global:\n")
  print(mejor_global[, intersect(c(
    "feature_set", "modelo", "CV_RMSE", "TEST_MAE", "TEST_RMSE", "TEST_R2", "best_params"
  ), names(mejor_global)), drop = FALSE])

} else {
  cat("AVISO: No se ha encontrado resumen_comparacion_modelos.csv. Salto tablas de modelos.\n")
}

# -------------------------------------------------------------------------
# 3. Métricas por proteína para el mejor modelo global
# -------------------------------------------------------------------------

if (!is.null(metricas) && !is.null(resumen_ordenado)) {
  metricas_mejor <- metricas[
    metricas$feature_set == mejor_feature_set &
      metricas$modelo == mejor_modelo,
    ,
    drop = FALSE
  ]

  metricas_mejor <- redondear_columnas(
    metricas_mejor,
    c("MAE", "RMSE", "R2"),
    digits = 4
  )

  guardar_csv(
    metricas_mejor,
    file.path(OUT_MODELOS, "tabla_metricas_mejor_modelo_por_proteina.csv")
  )

  guardar_csv(
    metricas_mejor,
    file.path(OUT_TABLAS, "tabla_04_metricas_mejor_modelo_por_proteina.csv")
  )

  if ("salida" %in% names(metricas_mejor)) {
    metricas_proteinas <- metricas_mejor[metricas_mejor$salida != "media", , drop = FALSE]
    guardar_csv(
      metricas_proteinas,
      file.path(OUT_TABLAS, "tabla_05_metricas_proteinas_sin_media.csv")
    )
  }

} else {
  cat("AVISO: No se ha encontrado metricas_por_proteina.csv o resumen. Salto métricas por proteína.\n")
}

# -------------------------------------------------------------------------
# 4. Tabla de errores/predicciones del mejor modelo
# -------------------------------------------------------------------------

if (!is.null(predicciones) && !is.null(resumen_ordenado)) {
  pred_mejor <- predicciones[
    predicciones$feature_set == mejor_feature_set &
      predicciones$modelo == mejor_modelo,
    ,
    drop = FALSE
  ]

  error_cols <- grep("^error_y_", names(pred_mejor), value = TRUE)

  if (length(error_cols) > 0) {
    pred_mejor$error_abs_medio <- rowMeans(abs(pred_mejor[, error_cols, drop = FALSE]), na.rm = TRUE)
    pred_mejor <- pred_mejor[order(-pred_mejor$error_abs_medio), ]
  }

  pred_mejor <- redondear_columnas(
    pred_mejor,
    grep("^(real_|pred_|error_|error_abs_medio)", names(pred_mejor), value = TRUE),
    digits = 4
  )

  guardar_csv(
    pred_mejor,
    file.path(OUT_TABLAS, "tabla_06_predicciones_test_mejor_modelo.csv")
  )

  if ("error_abs_medio" %in% names(pred_mejor)) {
    top_errores <- head(pred_mejor, 20)
    guardar_csv(
      top_errores,
      file.path(OUT_TABLAS, "tabla_07_top20_ligandos_mayor_error.csv")
    )
  }

} else {
  cat("AVISO: No se ha encontrado predicciones_test.csv o resumen. Salto predicciones.\n")
}

# -------------------------------------------------------------------------
# 5. Resumen del dataset
# -------------------------------------------------------------------------

if (!is.null(dataset)) {
  y_cols <- grep("^y_", names(dataset), value = TRUE)
  mem_cols <- grep("^mem_", names(dataset), value = TRUE)
  numeric_cols <- names(dataset)[sapply(dataset, is.numeric)]

  tabla_dataset <- data.frame(
    metrica = c(
      "n_filas",
      "n_columnas",
      "n_ligandos_unicos",
      "n_variables_numericas",
      "n_variables_membrana",
      "n_targets"
    ),
    valor = c(
      nrow(dataset),
      ncol(dataset),
      if ("ligando" %in% names(dataset)) length(unique(dataset$ligando)) else NA,
      length(numeric_cols),
      length(mem_cols),
      length(y_cols)
    )
  )

  guardar_csv(
    tabla_dataset,
    file.path(OUT_TABLAS, "tabla_08_resumen_dataset.csv")
  )

  missing <- data.frame(
    columna = names(dataset),
    n_missing = sapply(dataset, function(x) sum(is.na(x))),
    pct_missing = round(sapply(dataset, function(x) mean(is.na(x)) * 100), 2),
    stringsAsFactors = FALSE
  )

  missing <- missing[order(-missing$n_missing), ]

  guardar_csv(
    missing,
    file.path(OUT_TABLAS, "tabla_09_missing_por_columna.csv")
  )

} else {
  cat("AVISO: No se ha encontrado dataset_descriptores_membrana.csv. Salto resumen dataset.\n")
}

# -------------------------------------------------------------------------
# 6. Variables más correlacionadas con las salidas
# -------------------------------------------------------------------------

if (!is.null(correlaciones)) {
  if ("abs_spearman_media_y" %in% names(correlaciones)) {
    correlaciones$abs_spearman_media_y <- as.numeric(correlaciones$abs_spearman_media_y)
    correlaciones <- correlaciones[order(-correlaciones$abs_spearman_media_y), ]
  }

  top_correlaciones <- head(correlaciones, 30)
  top_correlaciones <- redondear_columnas(
    top_correlaciones,
    c("abs_spearman_media_y", "abs_spearman_max_y"),
    digits = 4
  )

  guardar_csv(
    top_correlaciones,
    file.path(OUT_TABLAS, "tabla_10_top30_variables_correlacion_targets.csv")
  )

} else {
  cat("AVISO: No se ha encontrado correlaciones_features_targets.csv. Salto correlaciones.\n")
}

# -------------------------------------------------------------------------
# 7. Nanopartícula-proteína
# -------------------------------------------------------------------------

if (!is.null(nano)) {
  energy_col <- NA_character_

  for (candidate in c("binding_energy_min", "E_nano_proteina", "energia", "energy")) {
    if (candidate %in% names(nano)) {
      energy_col <- candidate
      break
    }
  }

  if (is.na(energy_col)) {
    cat("AVISO: No encuentro columna de energía en nanopartículas. Salto tabla wide.\n")
  } else {
    nano_ok <- nano
    if ("estado" %in% names(nano_ok)) {
      nano_ok <- nano_ok[tolower(as.character(nano_ok$estado)) == "ok", , drop = FALSE]
    }

    if (!all(c("nanoparticula", "proteina") %in% names(nano_ok))) {
      cat("AVISO: Faltan columnas nanoparticula/proteina. Salto tabla nano.\n")
    } else {
      nano_ok[[energy_col]] <- as.numeric(nano_ok[[energy_col]])

      tabla_nano_long <- nano_ok[, intersect(c(
        "nanoparticula", "proteina", energy_col, "n_energias_detectadas", "archivo", "estado"
      ), names(nano_ok)), drop = FALSE]

      names(tabla_nano_long)[names(tabla_nano_long) == energy_col] <- "E_nano_proteina"

      tabla_nano_long <- tabla_nano_long[order(tabla_nano_long$proteina, tabla_nano_long$E_nano_proteina), ]
      tabla_nano_long <- redondear_columnas(tabla_nano_long, c("E_nano_proteina"), digits = 4)

      guardar_csv(
        tabla_nano_long,
        file.path(OUT_NANO, "tabla_nano_proteina_long.csv")
      )

      guardar_csv(
        tabla_nano_long,
        file.path(OUT_TABLAS, "tabla_11_nano_proteina_long.csv")
      )

      tmp <- tabla_nano_long[, c("nanoparticula", "proteina", "E_nano_proteina")]
      nano_wide <- reshape(
        tmp,
        idvar = "nanoparticula",
        timevar = "proteina",
        direction = "wide"
      )

      names(nano_wide) <- sub("^E_nano_proteina\\.", "", names(nano_wide))
      nano_wide <- redondear_columnas(nano_wide, setdiff(names(nano_wide), "nanoparticula"), digits = 4)

      guardar_csv(
        nano_wide,
        file.path(OUT_NANO, "tabla_nano_proteina_wide.csv")
      )

      guardar_csv(
        nano_wide,
        file.path(OUT_TABLAS, "tabla_12_nano_proteina_wide.csv")
      )

      proteinas <- unique(tabla_nano_long$proteina)
      mejor_nano <- do.call(
        rbind,
        lapply(proteinas, function(p) {
          sub <- tabla_nano_long[tabla_nano_long$proteina == p, , drop = FALSE]
          sub[which.min(sub$E_nano_proteina), ]
        })
      )

      guardar_csv(
        mejor_nano,
        file.path(OUT_TABLAS, "tabla_13_mejor_nanoparticula_por_proteina.csv")
      )
    }
  }

} else {
  cat("AVISO: No se ha encontrado resultados_docking_nanoparticulas_corregido.csv. Salto nanopartículas.\n")
}

# -------------------------------------------------------------------------
# 8. Ranking contextual
# -------------------------------------------------------------------------

if (!is.null(ranking)) {
  if ("score_contextual" %in% names(ranking)) {
    ranking$score_contextual <- as.numeric(ranking$score_contextual)
    ranking <- ranking[order(-ranking$score_contextual), ]
  }

  cols_ranking <- intersect(c(
    "ligando",
    "nanoparticula",
    "proteina",
    "y_ligando_proteina",
    "E_nano_proteina",
    "score_ligando_proteina",
    "score_nano_proteina",
    "score_contextual",
    "score_contextual_medio"
  ), names(ranking))

  ranking_limpio <- ranking[, cols_ranking, drop = FALSE]

  ranking_limpio <- redondear_columnas(
    ranking_limpio,
    setdiff(names(ranking_limpio), c("ligando", "nanoparticula", "proteina")),
    digits = 4
  )

  guardar_csv(
    ranking_limpio,
    file.path(OUT_NANO, "ranking_contextual_ligando_nano_proteina_ordenado.csv")
  )

  guardar_csv(
    ranking_limpio,
    file.path(OUT_TABLAS, "tabla_14_ranking_contextual_completo.csv")
  )

  top50 <- head(ranking_limpio, 50)

  guardar_csv(
    top50,
    file.path(OUT_NANO, "tabla_top50_ranking_contextual.csv")
  )

  guardar_csv(
    top50,
    file.path(OUT_TABLAS, "tabla_15_top50_ranking_contextual.csv")
  )

  if ("proteina" %in% names(ranking_limpio) && "score_contextual" %in% names(ranking_limpio)) {
    proteinas <- unique(ranking_limpio$proteina)
    top_por_proteina <- do.call(
      rbind,
      lapply(proteinas, function(p) {
        sub <- ranking_limpio[ranking_limpio$proteina == p, , drop = FALSE]
        head(sub[order(-as.numeric(sub$score_contextual)), ], 10)
      })
    )

    guardar_csv(
      top_por_proteina,
      file.path(OUT_TABLAS, "tabla_16_top10_ranking_por_proteina.csv")
    )
  }

  if ("nanoparticula" %in% names(ranking_limpio) && "score_contextual" %in% names(ranking_limpio)) {
    nanos <- unique(ranking_limpio$nanoparticula)
    top_por_nano <- do.call(
      rbind,
      lapply(nanos, function(n) {
        sub <- ranking_limpio[ranking_limpio$nanoparticula == n, , drop = FALSE]
        head(sub[order(-as.numeric(sub$score_contextual)), ], 10)
      })
    )

    guardar_csv(
      top_por_nano,
      file.path(OUT_TABLAS, "tabla_17_top10_ranking_por_nanoparticula.csv")
    )
  }

} else {
  cat("AVISO: No se ha encontrado ranking_contextual_ligando_nano_proteina.csv. Salto ranking contextual.\n")
}

# -------------------------------------------------------------------------
# 9. FIGURAS
# -------------------------------------------------------------------------

cat("\n=== Generando figuras ===\n")

# -------------------------------------------------------------------------
# Figura 1: Top modelos por RMSE test
# -------------------------------------------------------------------------

if (!is.null(resumen_ordenado)) {
  fig_path <- file.path(OUT_FIGS, "fig_01_top_modelos_rmse_test.png")
  abrir_png(fig_path, width = 2600, height = 1700, res = 300)

  top <- head(resumen_ordenado, 15)
  etiquetas <- paste(top$feature_set, top$modelo, sep = " + ")

  barplot(
    height = rev(top$TEST_RMSE),
    names.arg = rev(etiquetas),
    horiz = TRUE,
    las = 1,
    xlab = "RMSE en test",
    main = "Comparación de modelos por RMSE en test",
    cex.names = 0.75
  )

  cerrar_png(fig_path)
}

# -------------------------------------------------------------------------
# Figura 2: Mejor modelo por bloque de variables
# -------------------------------------------------------------------------

if (!is.null(resumen_ordenado)) {
  bloques <- unique(resumen_ordenado$feature_set)
  mejor_por_bloque_fig <- do.call(
    rbind,
    lapply(bloques, function(b) {
      sub <- resumen_ordenado[resumen_ordenado$feature_set == b, ]
      sub[which.min(sub$TEST_RMSE), ]
    })
  )

  mejor_por_bloque_fig <- mejor_por_bloque_fig[order(mejor_por_bloque_fig$TEST_RMSE), ]

  fig_path <- file.path(OUT_FIGS, "fig_02_mejor_rmse_por_bloque_variables.png")
  abrir_png(fig_path, width = 2300, height = 1600, res = 300)

  barplot(
    height = mejor_por_bloque_fig$TEST_RMSE,
    names.arg = mejor_por_bloque_fig$feature_set,
    las = 2,
    ylab = "RMSE en test",
    main = "Mejor RMSE por bloque de variables",
    cex.names = 0.8
  )

  cerrar_png(fig_path)
}

# -------------------------------------------------------------------------
# Figura 3: RMSE por proteína del mejor modelo
# -------------------------------------------------------------------------

if (!is.null(metricas_mejor) && "salida" %in% names(metricas_mejor)) {
  sub <- metricas_mejor[metricas_mejor$salida != "media", , drop = FALSE]

  if (nrow(sub) > 0 && "RMSE" %in% names(sub)) {
    fig_path <- file.path(OUT_FIGS, "fig_03_rmse_por_proteina_mejor_modelo.png")
    abrir_png(fig_path, width = 2000, height = 1500, res = 300)

    barplot(
      height = as.numeric(sub$RMSE),
      names.arg = sub$salida,
      las = 1,
      ylab = "RMSE en test",
      main = paste("RMSE por proteína ·", mejor_feature_set, "+", mejor_modelo)
    )

    cerrar_png(fig_path)
  }
}

# -------------------------------------------------------------------------
# Figura 4: R2 por proteína del mejor modelo
# -------------------------------------------------------------------------

if (!is.null(metricas_mejor) && "salida" %in% names(metricas_mejor)) {
  sub <- metricas_mejor[metricas_mejor$salida != "media", , drop = FALSE]

  if (nrow(sub) > 0 && "R2" %in% names(sub)) {
    fig_path <- file.path(OUT_FIGS, "fig_04_r2_por_proteina_mejor_modelo.png")
    abrir_png(fig_path, width = 2000, height = 1500, res = 300)

    barplot(
      height = as.numeric(sub$R2),
      names.arg = sub$salida,
      las = 1,
      ylab = expression(R^2),
      main = paste("R² por proteína ·", mejor_feature_set, "+", mejor_modelo)
    )
    abline(h = 0, lty = 2)

    cerrar_png(fig_path)
  }
}

# -------------------------------------------------------------------------
# Figuras 5.x: Real vs predicho para cada proteína
# -------------------------------------------------------------------------

if (!is.null(pred_mejor)) {
  targets <- c("y_1ao6", "y_1hzh", "y_2hav", "y_3ghg")

  for (target in targets) {
    real_col <- paste0("real_", target)
    pred_col <- paste0("pred_", target)

    if (all(c(real_col, pred_col) %in% names(pred_mejor))) {
      real <- as.numeric(pred_mejor[[real_col]])
      pred <- as.numeric(pred_mejor[[pred_col]])

      ok <- is.finite(real) & is.finite(pred)

      if (sum(ok) >= 2) {
        fig_path <- file.path(
          OUT_FIGS,
          paste0("fig_05_real_vs_predicho_", target, ".png")
        )

        abrir_png(fig_path, width = 1800, height = 1800, res = 300)
        par(mar = c(5, 5, 4, 2) + 0.1)

        plot(
          real[ok],
          pred[ok],
          pch = 19,
          xlab = "Afinidad real por docking",
          ylab = "Afinidad predicha",
          main = paste("Real vs predicho ·", target)
        )

        lim <- range(c(real[ok], pred[ok]), na.rm = TRUE)
        abline(a = 0, b = 1, lty = 2)
        grid()

        cerrar_png(fig_path)
      }
    }
  }
}

# -------------------------------------------------------------------------
# Figura 6: Top variables correlacionadas con targets
# -------------------------------------------------------------------------

if (!is.null(correlaciones) && "abs_spearman_media_y" %in% names(correlaciones)) {
  corr <- correlaciones
  corr$abs_spearman_media_y <- as.numeric(corr$abs_spearman_media_y)
  corr <- corr[order(-corr$abs_spearman_media_y), ]
  top <- head(corr, 20)

  name_col <- if ("columna" %in% names(top)) "columna" else names(top)[1]

  fig_path <- file.path(OUT_FIGS, "fig_06_top_variables_correlacion_targets.png")
  abrir_png(fig_path, width = 2600, height = 1700, res = 300)

  barplot(
    height = rev(top$abs_spearman_media_y),
    names.arg = rev(top[[name_col]]),
    horiz = TRUE,
    las = 1,
    xlab = "|Spearman| medio con las afinidades",
    main = "Variables más correlacionadas con las afinidades",
    cex.names = 0.75
  )

  cerrar_png(fig_path)
}

# -------------------------------------------------------------------------
# Figura 7: Heatmap nanopartícula-proteína
# -------------------------------------------------------------------------

if (!is.null(nano_wide)) {
  if ("nanoparticula" %in% names(nano_wide) && ncol(nano_wide) > 1) {
    mat <- as.matrix(nano_wide[, setdiff(names(nano_wide), "nanoparticula"), drop = FALSE])
    mode(mat) <- "numeric"
    rownames(mat) <- nano_wide$nanoparticula

    if (all(dim(mat) > 0)) {
      fig_path <- file.path(OUT_FIGS, "fig_07_heatmap_nano_proteina.png")
      abrir_png(fig_path, width = 1800, height = 1500, res = 300)

      par(mar = c(5, 8, 4, 5) + 0.1)
      image(
        x = seq_len(ncol(mat)),
        y = seq_len(nrow(mat)),
        z = t(mat[nrow(mat):1, , drop = FALSE]),
        axes = FALSE,
        xlab = "Proteína",
        ylab = "",
        main = "Energía de docking nanopartícula-proteína"
      )

      axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
      axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 1)

      for (i in seq_len(nrow(mat))) {
        for (j in seq_len(ncol(mat))) {
          text(j, nrow(mat) - i + 1, labels = round(mat[i, j], 2), cex = 0.9)
        }
      }

      cerrar_png(fig_path)
    }
  }
}

# -------------------------------------------------------------------------
# Figura 8: Barras de docking nano-proteína por proteína
# -------------------------------------------------------------------------

if (!is.null(tabla_nano_long) && all(c("nanoparticula", "proteina", "E_nano_proteina") %in% names(tabla_nano_long))) {
  proteinas <- unique(tabla_nano_long$proteina)

  for (p in proteinas) {
    sub <- tabla_nano_long[tabla_nano_long$proteina == p, , drop = FALSE]
    sub <- sub[order(sub$E_nano_proteina), ]

    fig_path <- file.path(
      OUT_FIGS,
      paste0("fig_08_docking_nano_proteina_", limpiar_nombre_archivo(p), ".png")
    )

    abrir_png(fig_path, width = 1800, height = 1400, res = 300)

    barplot(
      height = sub$E_nano_proteina,
      names.arg = sub$nanoparticula,
      las = 2,
      ylab = "Energía de docking",
      main = paste("Docking nanopartícula-proteína ·", p)
    )
    abline(h = 0, lty = 2)

    cerrar_png(fig_path)
  }
}

# -------------------------------------------------------------------------
# Figura 9: Top ranking contextual
# -------------------------------------------------------------------------

if (!is.null(ranking_limpio) && "score_contextual" %in% names(ranking_limpio)) {
  top <- head(ranking_limpio[order(-as.numeric(ranking_limpio$score_contextual)), ], 20)

  if (nrow(top) > 0) {
    etiquetas <- paste(top$ligando, top$nanoparticula, top$proteina, sep = " | ")

    fig_path <- file.path(OUT_FIGS, "fig_09_top20_ranking_contextual.png")
    abrir_png(fig_path, width = 3000, height = 1900, res = 300)

    barplot(
      height = rev(as.numeric(top$score_contextual)),
      names.arg = rev(etiquetas),
      horiz = TRUE,
      las = 1,
      xlab = "Score contextual",
      main = "Top 20 combinaciones ligando-nanopartícula-proteína",
      cex.names = 0.65
    )

    cerrar_png(fig_path)
  }
}

# -------------------------------------------------------------------------
# Figura 10: Distribución de afinidades y por proteína
# -------------------------------------------------------------------------

if (!is.null(dataset)) {
  y_cols <- grep("^y_", names(dataset), value = TRUE)

  if (length(y_cols) > 0) {
    fig_path <- file.path(OUT_FIGS, "fig_10_distribucion_afinidades_docking.png")
    abrir_png(fig_path, width = 2200, height = 1600, res = 300)

    vals <- unlist(dataset[, y_cols, drop = FALSE])
    vals <- as.numeric(vals)
    vals <- vals[is.finite(vals)]

    hist(
      vals,
      breaks = 20,
      xlab = "Afinidad de docking",
      main = "Distribución global de afinidades ligando-proteína"
    )

    cerrar_png(fig_path)

    # Boxplot por proteína
    fig_path <- file.path(OUT_FIGS, "fig_11_boxplot_afinidades_por_proteina.png")
    abrir_png(fig_path, width = 2100, height = 1600, res = 300)

    boxplot(
      dataset[, y_cols, drop = FALSE],
      las = 2,
      ylab = "Afinidad de docking",
      main = "Distribución de afinidades por proteína"
    )

    cerrar_png(fig_path)
  }
}

# -------------------------------------------------------------------------
# Figura 12: Membrana logP vs logPerm
# -------------------------------------------------------------------------

if (!is.null(dataset)) {
  if (all(c("mem_logP", "mem_logPerm") %in% names(dataset))) {
    x <- as.numeric(dataset$mem_logP)
    y <- as.numeric(dataset$mem_logPerm)
    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) >= 3) {
      fig_path <- file.path(OUT_FIGS, "fig_12_membrana_logp_vs_logperm.png")
      abrir_png(fig_path, width = 1800, height = 1700, res = 300)
      par(mar = c(5, 5, 4, 2) + 0.1)

      plot(
        x[ok],
        y[ok],
        pch = 19,
        xlab = "mem_logP",
        ylab = "mem_logPerm",
        main = "Relación entre logP de membrana y permeabilidad"
      )
      grid()

      if ("ligando" %in% names(dataset)) {
        # Etiquetar solo algunos extremos para que no quede ilegible
        ord <- order(abs(scale(x[ok])) + abs(scale(y[ok])), decreasing = TRUE)
        idx_ok <- which(ok)
        idx_label <- idx_ok[head(ord, min(8, length(ord)))]
        text(
          x[idx_label],
          y[idx_label],
          labels = dataset$ligando[idx_label],
          pos = 4,
          cex = 0.55
        )
      }

      cerrar_png(fig_path)
    }
  }
}

# -------------------------------------------------------------------------
# 10. Índices de salidas
# -------------------------------------------------------------------------

tablas_generadas <- list.files(OUT_TABLAS, pattern = "\\.csv$", full.names = FALSE)
figuras_generadas <- list.files(OUT_FIGS, pattern = "\\.png$", full.names = FALSE)

indice_tablas <- data.frame(
  archivo = tablas_generadas,
  ruta = file.path(OUT_TABLAS, tablas_generadas),
  stringsAsFactors = FALSE
)

indice_figuras <- data.frame(
  archivo = figuras_generadas,
  ruta = file.path(OUT_FIGS, figuras_generadas),
  stringsAsFactors = FALSE
)

guardar_csv(
  indice_tablas,
  file.path(OUT_TABLAS, "indice_tablas_generadas.csv")
)

guardar_csv(
  indice_figuras,
  file.path(OUT_FIGS, "indice_figuras_generadas.csv")
)

cat("\n=== FIN ===\n")
cat("Tablas generadas en:", OUT_TABLAS, "\n")
cat("Figuras generadas en:", OUT_FIGS, "\n")
cat("Total tablas:", nrow(indice_tablas), "\n")
cat("Total figuras:", nrow(indice_figuras), "\n\n")
