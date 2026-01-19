#!/usr/bin/env Rscript
# ============================================================
# run_topkat_job.R
# Ejecuta TopKAT + scale_importance a partir de un subset ya creado
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(TopKAT)
})

# ============================================================
# === PARSEO DE ARGUMENTOS ===
# ============================================================
option_list <- list(
  make_option(c("-s", "--subset_dir"),
              type = "character",
              help = "Directorio del subset_TOPKAT_* generado por la app",
              metavar = "PATH"),
  make_option(c("-o", "--omega"),
              type = "character",
              default = "0,0.5,1",
              help = "Valores de omega separados por coma [default %default]"),
  make_option(c("-t", "--threshold"),
              type = "numeric",
              default = 500,
              help = "Threshold para scale_importance [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$subset_dir)) {
  stop("ERROR: Debes especificar --subset_dir")
}

subset_dir <- normalizePath(opt$subset_dir)

# ============================================================
# === ARCHIVOS ESPERADOS ===
# ============================================================
rips_file <- file.path(subset_dir, "rips_list_subset.rds")
K0_file   <- file.path(subset_dir, "K_dim0_subset.rds")
K1_file   <- file.path(subset_dir, "K_dim1_subset.rds")
PIDs_file <- file.path(subset_dir, "PIDs_subset.rds")
pid_csv   <- file.path(subset_dir, "PID_seleccionados_final.csv")

required_files <- c(rips_file, K0_file, K1_file, PIDs_file, pid_csv)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop(
    "ERROR: Faltan archivos en el subset:\n",
    paste(missing, collapse = "\n")
  )
}

cat("✔ Subset válido encontrado en:\n", subset_dir, "\n")

# ============================================================
# === CARGA DE DATOS ===
# ============================================================
cat("→ Cargando datos...\n")

rips_sub <- readRDS(rips_file)
K0_sub   <- readRDS(K0_file)
K1_sub   <- readRDS(K1_file)
PIDs_sub <- readRDS(PIDs_file)

pid_df <- read.csv(pid_csv, stringsAsFactors = FALSE)

if (!"grupo" %in% colnames(pid_df)) {
  stop("ERROR: El CSV no tiene columna 'grupo'")
}

pid_df <- pid_df %>%
  arrange(match(as.character(PID), as.character(PIDs_sub)))

y <- ifelse(pid_df$grupo == "A", 1, 0)

# ============================================================
# === TOPKAT ===
# ============================================================
omega.list <- as.numeric(strsplit(opt$omega, ",")[[1]])
K.list <- list(K0 = K0_sub, K1 = K1_sub)

cat("→ Ejecutando TopKAT...\n")

res_topkat <- TopKAT(
  y = y,
  X = NULL,
  K.list = K.list,
  omega.list = omega.list,
  outcome.type = "binary"
)

saveRDS(res_topkat, file.path(subset_dir, "TopKAT_result.rds"))

cat("✔ TopKAT terminado\n")
cat("  p-value global:", signif(res_topkat$overall.pval, 4), "\n")

# ============================================================
# === SCALE IMPORTANCE ===
# ============================================================
cat("→ Ejecutando scale_importance...\n")

res_scale <- scale_importance(
  pd.list    = rips_sub,
  y          = y,
  omega.list = omega.list,
  threshold  = opt$threshold,
  PIDs       = seq_along(PIDs_sub),
  outcome.type = "binary"
)

saveRDS(res_scale, file.path(subset_dir, "scale_importance_result.rds"))

cat("✔ scale_importance terminado\n")
cat("  min.thresh:", res_scale$min.thresh, "\n")

# ============================================================
# === GRAFICAR RESULTADOS ===
# ============================================================
df_plot <- data.frame(
  epsilon = res_scale$threshold.seq,
  pval    = res_scale$pvals
)

# --- PNG normal ---
png_normal <- file.path(subset_dir, "scale_importance_normal.png")
png(png_normal, width = 1800, height = 1200, res = 200)
print(
  ggplot(df_plot, aes(epsilon, pval)) +
    geom_point() +
    theme_bw() +
    geom_vline(xintercept = res_scale$min.thresh, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
    labs(x = expression(epsilon), y = "p-valor")
)
dev.off()

# --- PNG log10 ---
png_log <- file.path(subset_dir, "scale_importance_log10.png")
png(png_log, width = 1800, height = 1200, res = 200)
print(
  ggplot(df_plot, aes(epsilon, pval)) +
    geom_point() +
    scale_y_log10() +
    theme_bw() +
    geom_vline(xintercept = res_scale$min.thresh, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
    labs(x = expression(epsilon), y = expression(log[10](p)))
)
dev.off()

cat("✔ PNGs generados\n")
cat("✔ Job finalizado correctamente\n")
