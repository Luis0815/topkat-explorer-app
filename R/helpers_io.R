# R/helpers_io.R
library(readr)
library(tools)

safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("Archivo CSV no encontrado: ", path)
  ext <- tolower(file_ext(path))
  if (ext == "gz") {
    read_csv(path)
  } else {
    read_csv(path)
  }
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) stop("Archivo RDS no encontrado: ", path)
  readRDS(path)
}

safe_save_rds <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  message("✔ Guardado: ", path)
}

safe_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_csv(df, path)
  message("✔ Guardado CSV: ", path)
}
