# R/helpers_logic.R
library(dplyr)
library(stringr)
library(purrr)

# Clasificación extendida similar a tu R
clasificar_categoria_v2 <- function(nombre) {
  n <- tolower(nombre %||% "")
  n <- str_replace_all(n, "[\\-\\s]+", "_")
  if (str_detect(n, "stroma_ad_carcinoma") || str_detect(n, "stroma_ad_carc")) return("stroma_ad_carcinoma")
  if (str_detect(n, "stroma_ad_dysplasia") || str_detect(n, "stroma_ad_dyspl") || str_detect(n, "stroma_ad_dys")) return("stroma_ad_dysplasia")
  if (str_detect(n, "stroma") && str_detect(n, "carcinoma")) return("stroma_ad_carcinoma")
  if (str_detect(n, "stroma") && str_detect(n, "dysplasia")) return("stroma_ad_dysplasia")
  if (str_detect(n, "carcinoma")) return("carcinoma")
  if (str_detect(n, "dysplasia")) return("dysplasia")
  return("otro")
}

get_FA_status <- function(nombre) {
  n <- toupper(nombre %||% "")
  if (str_detect(n, "_F") || str_detect(n, "_FA") || str_detect(n, " FA") || str_detect(n, "\\(F\\)") ) return("FA")
  return("NOFA")
}

# Construir pid_info (una fila por PID)
build_pid_info <- function(data1.df, archivo_col = "archivo_base") {
  pid_info <- data1.df %>%
    group_by(PID) %>%
    summarise(archivo_base = unique(.data[[archivo_col]])[1], .groups = "drop") %>%
    mutate(
      category_raw = archivo_base,
      category = sapply(archivo_base, clasificar_categoria_v2),
      FA_status = sapply(archivo_base, get_FA_status)
    )
  pid_info
}

# Dadas las definiciones de comparaciones (lista), filtra disponibilidad
filter_by_cmp <- function(pid_info, cmpA, cmpB) {
  A <- pid_info %>% filter(category == cmpA$type & (is.null(cmpA$FA) || FA_status == cmpA$FA))
  B <- pid_info %>% filter(category == cmpB$type & (is.null(cmpB$FA) || FA_status == cmpB$FA))
  list(A = A, B = B)
}

# Selección reproducible
select_reproducible <- function(df, n, seed = 123) {
  if (nrow(df) < n) stop("No hay suficientes muestras")
  set.seed(seed)
  df %>% slice_sample(n = n)
}

# Crear subset rips + kernels + PIDs
create_subset_files <- function(selected_pids_chr, rips_full, K0_full, K1_full, PIDs_full, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  # Check existence
  missing <- setdiff(selected_pids_chr, PIDs_full)
  if (length(missing) > 0) stop("Some selected PIDs not in PIDs_full: ", paste(missing, collapse = ", "))
  # rips_full is a named list with names = PIDs_full
  rips_subset <- rips_full[selected_pids_chr]
  K0_subset <- K0_full[selected_pids_chr, selected_pids_chr, drop = FALSE]
  K1_subset <- K1_full[selected_pids_chr, selected_pids_chr, drop = FALSE]
  saveRDS(rips_subset, file.path(outdir, "rips_list_subset.rds"))
  saveRDS(K0_subset, file.path(outdir, "K_dim0_subset.rds"))
  saveRDS(K1_subset, file.path(outdir, "K_dim1_subset.rds"))
  saveRDS(selected_pids_chr, file.path(outdir, "PIDs_subset.rds"))
  return(list(rips = rips_subset, K0 = K0_subset, K1 = K1_subset, PIDs = selected_pids_chr))
}
