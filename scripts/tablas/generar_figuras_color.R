options(stringsAsFactors = FALSE)

cat("\n=== Generador de figuras con color TFM ===\n")
cat("Directorio de trabajo:", getwd(), "\n\n")

OUT_FIGS <- "resultados/figuras_generadas"
if (!dir.exists(OUT_FIGS)) dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)

# Paleta sobria para memoria
COL_BLUE <- "#12355B"
COL_GOLD <- "#D8AD6E"
COL_GREEN <- "#2E8B57"
COL_RED <- "#B23A48"
COL_PURPLE <- "#6A4C93"
COL_CYAN <- "#2A9D8F"
COL_GREY <- "#52616B"
COL_LIGHT <- "#F6F8FA"
PALETTE <- c(COL_BLUE, COL_GOLD, COL_GREEN, COL_RED, COL_PURPLE, COL_CYAN, COL_GREY)

buscar <- function(candidatos) {
  for (x in candidatos) if (file.exists(x)) return(x)
  return(NA_character_)
}

leer <- function(candidatos) {
  p <- buscar(candidatos)
  if (is.na(p)) {
    cat("AVISO: no encuentro:", paste(candidatos, collapse = " | "), "\n")
    return(NULL)
  }
  cat("Leyendo:", p, "\n")
  read.csv(p, check.names = FALSE)
}

png_open <- function(path, width = 2400, height = 1600, res = 300) {
  png(path, width = width, height = height, res = res)
  par(bg = "white", mar = c(9, 5, 4, 2) + 0.1)
}

png_close <- function(path) {
  dev.off()
  cat("Figura creada:", path, "\n")
}

short <- function(x, n = 36) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

# Etiquetas más legibles para figuras de la memoria
label_bloque <- function(x) {
  x <- as.character(x)
  d <- c(
    "RDKit_solo" = "RDKit solo",
    "RDKit_membrana_seleccionada" = "RDKit + membrana",
    "RDKit_membrana_filtrada_auto" = "RDKit + membrana filtrada",
    "Membrana_solo" = "Membrana solo"
  )
  out <- x
  m <- x %in% names(d)
  out[m] <- d[x[m]]
  out
}

label_bloque_eje <- function(x) {
  x <- as.character(x)
  d <- c(
    "RDKit_solo" = "RDKit\nsolo",
    "RDKit_membrana_seleccionada" = "RDKit +\nmembrana",
    "RDKit_membrana_filtrada_auto" = "RDKit +\nmembrana filtrada",
    "Membrana_solo" = "Membrana\nsolo"
  )
  out <- x
  m <- x %in% names(d)
  out[m] <- d[x[m]]
  out
}

label_bloque_corto <- function(x) {
  x <- as.character(x)
  d <- c(
    "RDKit_solo" = "RDKit",
    "RDKit_membrana_seleccionada" = "RDKit+mem",
    "RDKit_membrana_filtrada_auto" = "RDKit+mem filt.",
    "Membrana_solo" = "Membrana"
  )
  out <- x
  m <- x %in% names(d)
  out[m] <- d[x[m]]
  out
}

label_variable_corto <- function(x) {
  x <- as.character(x)
  x <- gsub("^mem_", "mem_", x)
  x <- gsub("diffusion_E9", "diff_E9", x)
  x <- gsub("meanentropy", "entropy", x)
  x <- gsub("free_energy", "free_E", x)
  x <- gsub("distribution", "distrib", x)
  x <- gsub("InertialShapeFactor", "InertialShape", x)
  x
}

label_np_plotmath <- function(x) {
  labs <- lapply(as.character(x), function(v) {
    if (v == "Ag_13") {
      quote((plain(Ag))[13])
    } else if (v == "Au_13") {
      quote((plain(Au))[13])
    } else if (v == "ZnO_12_0") {
      quote((plain(ZnO))[12])
    } else {
      as.character(v)
    }
  })
  as.expression(labs)
}
num <- function(x) as.numeric(as.character(x))

# -------------------------------------------------------------------------
# Cargar datos
# -------------------------------------------------------------------------

resumen <- leer(c(
  "resultados/modelos/resultados_comparacion_membrana/resumen_comparacion_modelos.csv",
  "resultados/modelos/resumen_comparacion_modelos.csv",
  "resultados_comparacion_membrana/resumen_comparacion_modelos.csv"
))

metricas <- leer(c(
  "resultados/modelos/resultados_comparacion_membrana/metricas_por_proteina.csv",
  "resultados/modelos/metricas_por_proteina.csv",
  "resultados_comparacion_membrana/metricas_por_proteina.csv"
))

predicciones <- leer(c(
  "resultados/modelos/resultados_comparacion_membrana/predicciones_test.csv",
  "resultados/modelos/predicciones_test.csv",
  "resultados_comparacion_membrana/predicciones_test.csv"
))

dataset <- leer(c(
  "datasets/dataset_descriptores_membrana.csv",
  "dataset_descriptores_membrana.csv"
))

nano <- leer(c(
  "resultados/nanoparticulas/resultados_docking_nanoparticulas_corregido.csv",
  "resultados_docking_nanoparticulas_corregido.csv"
))

ranking <- leer(c(
  "resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina.csv",
  "ranking_contextual_ligando_nano_proteina.csv"
))

correlaciones <- leer(c(
  "resultados/diagnosticos/correlaciones_features_targets.csv",
  "correlaciones_features_targets.csv"
))

# Determinar mejor modelo
mejor_feature_set <- NA_character_
mejor_modelo <- NA_character_
resumen_ord <- NULL
if (!is.null(resumen) && all(c("feature_set", "modelo", "TEST_RMSE") %in% names(resumen))) {
  resumen$TEST_RMSE <- num(resumen$TEST_RMSE)
  if ("TEST_R2" %in% names(resumen)) resumen$TEST_R2 <- num(resumen$TEST_R2)
  if ("TEST_MAE" %in% names(resumen)) resumen$TEST_MAE <- num(resumen$TEST_MAE)
  if ("CV_RMSE" %in% names(resumen)) resumen$CV_RMSE <- num(resumen$CV_RMSE)
  resumen_ord <- resumen[order(resumen$TEST_RMSE), ]
  mejor_feature_set <- resumen_ord$feature_set[1]
  mejor_modelo <- resumen_ord$modelo[1]
}

# -------------------------------------------------------------------------
# 1. Top modelos por RMSE, color por bloque de variables
# -------------------------------------------------------------------------

if (!is.null(resumen_ord)) {
  # Se muestran menos modelos y con etiquetas más grandes para que el eje Y
  # sea legible al insertar la figura en el PDF de la memoria.
  top <- head(resumen_ord, 12)
  labels <- paste(label_bloque_corto(top$feature_set), top$modelo, sep = " + ")
  labels <- short(labels, 34)

  grupos <- as.factor(top$feature_set)
  cols <- PALETTE[as.integer(grupos)]

  path <- file.path(OUT_FIGS, "fig_01_top_modelos_rmse_test.png")
  png_open(path, width = 3600, height = 2300)
  par(
    mar = c(5.8, 20.0, 4.8, 2.2) + 0.1,
    xaxs = "i",
    yaxs = "i",
    cex.axis = 1.15,
    cex.lab = 1.25,
    cex.main = 1.30
  )

  barplot(
    rev(top$TEST_RMSE),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    col = rev(cols),
    border = NA,
    xlab = "RMSE en test",
    main = "Comparación de modelos por RMSE en test",
    cex.names = 1.15,
    xlim = c(0, max(top$TEST_RMSE, na.rm = TRUE) * 1.26)
  )

  grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
  legend(
    "topright",
    legend = label_bloque(levels(grupos)),
    fill = PALETTE[seq_along(levels(grupos))],
    bty = "n",
    bg = "white",
    cex = 1.05,
    inset = 0.02
  )

  png_close(path)
}

# -------------------------------------------------------------------------
# 2. Mejor modelo por bloque de variables
# -------------------------------------------------------------------------

if (!is.null(resumen_ord)) {
  bloques <- unique(resumen_ord$feature_set)
  best <- do.call(rbind, lapply(bloques, function(b) {
    s <- resumen_ord[resumen_ord$feature_set == b, ]
    s[which.min(s$TEST_RMSE), ]
  }))
  best <- best[order(best$TEST_RMSE), ]

  path <- file.path(OUT_FIGS, "fig_02_mejor_rmse_por_bloque_variables.png")
  png_open(path, width = 2000, height = 1450)
  par(mar = c(6.2, 5.2, 4.0, 1.5) + 0.1)

  barplot(
    best$TEST_RMSE,
    names.arg = label_bloque_eje(best$feature_set),
    las = 1,
    col = PALETTE[seq_len(nrow(best))],
    border = NA,
    ylab = "RMSE en test",
    main = "Mejor RMSE por bloque de variables",
    cex.names = 0.78,
    ylim = c(0, max(best$TEST_RMSE, na.rm = TRUE) * 1.15)
  )
  grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)

  png_close(path)
}

# -------------------------------------------------------------------------
# 3 y 4. Métricas por proteína
# -------------------------------------------------------------------------

if (!is.null(metricas) && !is.na(mejor_feature_set)) {
  sub <- metricas[metricas$feature_set == mejor_feature_set & metricas$modelo == mejor_modelo, ]
  sub <- sub[sub$salida != "media", ]
  if (nrow(sub) > 0) {
    sub$RMSE <- num(sub$RMSE)
    sub$R2 <- num(sub$R2)
    sub$MAE <- num(sub$MAE)

    path <- file.path(OUT_FIGS, "fig_03_rmse_por_proteina_mejor_modelo.png")
    png_open(path, width = 2000, height = 1500)
    barplot(sub$RMSE, names.arg = sub$salida, col = COL_BLUE, border = NA,
            ylab = "RMSE", main = paste("RMSE por proteína ·", mejor_modelo))
    grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
    png_close(path)

    path <- file.path(OUT_FIGS, "fig_04_r2_por_proteina_mejor_modelo.png")
    png_open(path, width = 2000, height = 1500)
    cols <- ifelse(sub$R2 >= 0, COL_GREEN, COL_RED)
    barplot(sub$R2, names.arg = sub$salida, col = cols, border = NA,
            ylab = "R²", main = paste("R² por proteína ·", mejor_modelo))
    abline(h = 0, lty = 2, col = COL_GREY)
    grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
    png_close(path)
  }
}

# -------------------------------------------------------------------------
# 5. Real vs predicho por proteína
# -------------------------------------------------------------------------

if (!is.null(predicciones) && !is.na(mejor_feature_set)) {
  pred <- predicciones[predicciones$feature_set == mejor_feature_set & predicciones$modelo == mejor_modelo, ]
  targets <- c("y_1ao6", "y_1hzh", "y_2hav", "y_3ghg")

  for (i in seq_along(targets)) {
    target <- targets[i]
    rc <- paste0("real_", target)
    pc <- paste0("pred_", target)

    if (all(c(rc, pc) %in% names(pred))) {
      real <- num(pred[[rc]])
      pr <- num(pred[[pc]])
      ok <- is.finite(real) & is.finite(pr)
      if (sum(ok) >= 2) {
        path <- file.path(OUT_FIGS, paste0("fig_05_real_vs_predicho_", target, ".png"))
        png_open(path, width = 1800, height = 1800)
        par(mar = c(5,5,4,2)+0.1)
        plot(real[ok], pr[ok], pch = 19, col = PALETTE[i],
             xlab = "Afinidad real", ylab = "Afinidad predicha",
             main = paste("Real vs predicho ·", target))
        lim <- range(c(real[ok], pr[ok]), na.rm = TRUE)
        abline(a = 0, b = 1, lty = 2, col = COL_GREY, lwd = 2)
        grid(col = "#DDDDDD", lty = 3)
        png_close(path)
      }
    }
  }
}

# -------------------------------------------------------------------------
# 6. Variables más correlacionadas
# -------------------------------------------------------------------------
# 6. Top variables correlacionadas con las afinidades
# -------------------------------------------------------------------------

if (!is.null(correlaciones) && "abs_spearman_media_y" %in% names(correlaciones)) {
  correlaciones$abs_spearman_media_y <- num(correlaciones$abs_spearman_media_y)
  corr <- correlaciones[order(-correlaciones$abs_spearman_media_y), ]
  # Menos variables y texto más grande para que las etiquetas del eje Y
  # sean legibles en la memoria.
  top <- head(corr, 14)
  name_col <- if ("columna" %in% names(top)) "columna" else names(top)[1]

  is_mem <- if ("es_membrana" %in% names(top)) as.integer(top$es_membrana) == 1 else grepl("^mem_", top[[name_col]])
  cols <- ifelse(is_mem, COL_GOLD, COL_BLUE)

  path <- file.path(OUT_FIGS, "fig_06_top_variables_correlacion_targets.png")
  png_open(path, width = 3600, height = 2300)
  par(
    mar = c(5.8, 18.5, 4.8, 2.5) + 0.1,
    xaxs = "i",
    yaxs = "i",
    cex.axis = 1.15,
    cex.lab = 1.25,
    cex.main = 1.30
  )

  barplot(
    rev(top$abs_spearman_media_y),
    names.arg = rev(short(label_variable_corto(top[[name_col]]), 28)),
    horiz = TRUE,
    las = 1,
    col = rev(cols),
    border = NA,
    xlab = "|Spearman| medio",
    main = "Variables más correlacionadas con las afinidades",
    cex.names = 1.15,
    xlim = c(0, max(top$abs_spearman_media_y, na.rm = TRUE) * 1.25)
  )

  legend(
    "bottomright",
    legend = c("RDKit/geométricas", "Membrana"),
    fill = c(COL_BLUE, COL_GOLD),
    bty = "n",
    bg = "white",
    cex = 1.05
  )
  grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
  png_close(path)
}

# -------------------------------------------------------------------------
# 7 y 8. Nanopartícula-proteína
# -------------------------------------------------------------------------

nano_long <- NULL
if (!is.null(nano)) {
  e_col <- NA_character_
  for (cc in c("binding_energy_min", "E_nano_proteina", "energia", "energy")) if (cc %in% names(nano)) e_col <- cc
  if (!is.na(e_col) && all(c("nanoparticula", "proteina") %in% names(nano))) {
    nano_long <- nano
    if ("estado" %in% names(nano_long)) nano_long <- nano_long[tolower(as.character(nano_long$estado)) == "ok", ]
    nano_long$E <- num(nano_long[[e_col]])

    # heatmap con colores custom
    prot <- unique(nano_long$proteina)
    naps <- unique(nano_long$nanoparticula)
    mat <- matrix(NA_real_, nrow = length(naps), ncol = length(prot), dimnames = list(naps, prot))
    for (i in seq_len(nrow(nano_long))) mat[nano_long$nanoparticula[i], nano_long$proteina[i]] <- nano_long$E[i]

    path <- file.path(OUT_FIGS, "fig_07_heatmap_nano_proteina.png")
    png_open(path, width = 1700, height = 1300)
    par(mar = c(4.8, 7.2, 3.8, 4.2) + 0.1)
    pal <- colorRampPalette(c(COL_BLUE, "white", COL_GOLD))(80)
    image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1,,drop=FALSE]),
          axes = FALSE, col = pal, xlab = "Proteína", ylab = "", main = "Docking nanopartícula-proteína")
    axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), cex.axis = 0.85)
    axis(
      2,
      at = seq_len(nrow(mat)),
      labels = label_np_plotmath(rev(rownames(mat))),
      las = 1,
      cex.axis = 1.05
    )    
    for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) text(j, nrow(mat)-i+1, round(mat[i,j],2), cex=.85)
    png_close(path)

    for (p in prot) {
      s <- nano_long[nano_long$proteina == p, ]
      s <- s[order(s$E), ]
      path <- file.path(OUT_FIGS, paste0("fig_08_docking_nano_proteina_", p, ".png"))
      png_open(path, width = 1800, height = 1400)
      barplot(s$E, names.arg = s$nanoparticula, las = 2, col = PALETTE[seq_len(nrow(s))], border = NA,
              ylab = "Energía de docking", main = paste("Docking nanopartícula-proteína ·", p))
      abline(h = 0, lty = 2, col = COL_GREY)
      grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
      png_close(path)
    }
  }
}

# -------------------------------------------------------------------------
# 9. Ranking contextual
# -------------------------------------------------------------------------

if (!is.null(ranking) && "score_contextual" %in% names(ranking)) {
  ranking$score_contextual <- num(ranking$score_contextual)
  top <- head(ranking[order(-ranking$score_contextual), ], 20)
  labs <- short(paste(top$ligando, top$nanoparticula, top$proteina, sep = " | "), 62)
  cols <- PALETTE[as.integer(as.factor(top$nanoparticula))]

  path <- file.path(OUT_FIGS, "fig_09_top20_ranking_contextual.png")
  png_open(path, width = 3100, height = 1950)
  barplot(rev(top$score_contextual), names.arg = rev(labs), horiz = TRUE, las = 1,
          col = rev(cols), border = NA, xlab = "Score contextual",
          main = "Top 20 combinaciones ligando-nanopartícula-proteína",
          cex.names = 0.65)
  legend("bottomright", legend = levels(as.factor(top$nanoparticula)), fill = PALETTE[seq_along(levels(as.factor(top$nanoparticula)))], bty = "n", cex = .8)
  grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
  png_close(path)
}

# -------------------------------------------------------------------------
# 10 y 11. Distribución de afinidades
# -------------------------------------------------------------------------

if (!is.null(dataset)) {
  y_cols <- grep("^y_", names(dataset), value = TRUE)
  if (length(y_cols) > 0) {
    vals <- as.numeric(unlist(dataset[, y_cols, drop = FALSE])); vals <- vals[is.finite(vals)]
    path <- file.path(OUT_FIGS, "fig_10_distribucion_afinidades_docking.png")
    png_open(path, width = 2200, height = 1600)
    hist(vals, breaks = 20, col = COL_BLUE, border = "white", xlab = "Afinidad de docking", main = "Distribución global de afinidades")
    grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
    png_close(path)

    path <- file.path(OUT_FIGS, "fig_11_boxplot_afinidades_por_proteina.png")
    png_open(path, width = 2100, height = 1600)
    boxplot(dataset[, y_cols, drop=FALSE], las = 2, col = PALETTE[seq_along(y_cols)], border = COL_GREY,
            ylab = "Afinidad de docking", main = "Distribución de afinidades por proteína")
    grid(nx = NA, ny = NULL, col = "#DDDDDD", lty = 3)
    png_close(path)
  }
}

# -------------------------------------------------------------------------
# 12. Membrana logP vs logPerm
# -------------------------------------------------------------------------

if (!is.null(dataset) && all(c("mem_logP", "mem_logPerm") %in% names(dataset))) {
  x <- num(dataset$mem_logP); y <- num(dataset$mem_logPerm)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 3) {
    path <- file.path(OUT_FIGS, "fig_12_membrana_logp_vs_logperm.png")
    png_open(path, width = 1650, height = 1500)
    par(mar = c(5,5,4,2)+0.1)
    plot(x[ok], y[ok], pch = 19, col = COL_GOLD, xlab = "mem_logP", ylab = "mem_logPerm", main = "Relación entre logP de membrana y permeabilidad")
    grid(col = "#DDDDDD", lty = 3)
    fit <- lm(y[ok] ~ x[ok]); abline(fit, col = COL_BLUE, lwd = 2)
    png_close(path)
  }
}

# índice
figs <- list.files(OUT_FIGS, pattern = "\\.png$", full.names = FALSE)
write.csv(data.frame(archivo = figs, ruta = file.path(OUT_FIGS, figs)), file.path(OUT_FIGS, "indice_figuras_generadas.csv"), row.names = FALSE)
cat("\nFIN. Figuras creadas en:", OUT_FIGS, "\n")
