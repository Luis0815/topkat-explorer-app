# TopKAT + Connectivity Explorer

Aplicación **Shiny** que combina dos análisis para estudios de **organización espacial celular**.

---

##  Nota importante sobre los datos

Por razones de tamaño, **los datos no se versionan directamente en este repositorio**.

* Las carpetas `data_selection/` y `data_circulos/` **  contienen archivos comprimidos (`.gz`)**..

Hay que descomprimirlos con gunzip *.gz

### 📥 Descarga obligatoria de `data_selection/`

La carpeta `data_selection/` **debe descargarse externamente** y luego descomprimirse de forma local.

Descárgala desde el siguiente enlace de Google Drive:

👉 [https://drive.google.com/drive/folders/1g9hDSDqCRh4XwGtCR1__IxoikGjuidRk?usp=drive_link](https://drive.google.com/drive/folders/1g9hDSDqCRh4XwGtCR1__IxoikGjuidRk?usp=drive_link)

Una vez descargada:

1. Copia la carpeta `data_selection/` dentro del repositorio.
2. Descomprime los archivos `.gz` manteniendo los nombres originales esperados por la aplicación.

>  **No subas estos archivos al repositorio**, ya que exceden los límites de GitHub.

---

## 🧬 Módulo **TopKAT Explorer**

Permite:

* Seleccionar la carpeta de datos de entrada que contiene los archivos requeridos para el análisis.
* Ejecutar **TopKAT** y **scale_importance** sobre diagramas de persistencia (PDs).
* Seleccionar automáticamente **subconjuntos de muestras** según comparaciones definidas.
* Generar carpetas de trabajo reproducibles con todos los resultados (**RDS + PNG + CSV**).
* Visualizar resultados dentro de la aplicación.

---

## 🔗 Módulo de **Conectividad Celular**

* Calcula matrices de **conectividad promedio** entre tipos celulares usando distancias espaciales *(x,y)*.
* Produce **heatmaps comparativos** entre Grupo A y Grupo B con **misma escala global**.
* Guarda las imágenes dentro del mismo *subset* generado por TopKAT.

---

# ✨ Características principales

## ✔ TopKAT Explorer

* Carga metadata y archivos RDS (rips, kernels, PIDs).
* Permite elegir una comparación predefinida (FA vs NOFA, carcinoma vs dysplasia, etc.).
* Selecciona **nA y nB** muestras con **semillas reproducibles**.
* Genera un directorio:

```
subset_TOPKAT_YYYYMMDD_HHMMSS/
```

Con:

* `rips_list_subset.rds`
* `K_dim0_subset.rds`
* `K_dim1_subset.rds`
* `PIDs_subset.rds`
* `PID_seleccionados_final.csv`
* `scale_importance_normal.png`
* `scale_importance_log10.png`

---

## ✔ Conectividad Celular

Para cada PID:

* Calcula distancias euclidianas.
* Determina conexiones (`dist ≤ ε`).
* Construye matrices **A** y **B** promedio.

Genera:

* `connectivity_groupA.svg`
* `connectivity_groupB.svg`
* `connectivity_A_vs_B.svg`

Los heatmaps muestran:

* **Eje X:** tipo celular 1
* **Eje Y:** tipo celular 2
* **Fill:** número de conexiones

Todo se guarda directamente en el subset generado por TopKAT.

---

# 📦 Requisitos

## 🔹 Versión recomendada de R

**R ≥ 4.2**

## 🔹 Librerías necesarias

```r
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(readr)
library(purrr)
library(viridis)
library(igraph)
library(here)
library(base64enc)
library(TopKAT)
```

### Instalación rápida

```r
install.packages(c(
  "shiny", "shinyWidgets", "shinycssloaders", "dplyr", "tidyr",
  "ggplot2", "reshape2", "stringr", "readr", "purrr", "viridis",
  "igraph", "here", "base64enc","svglite","DT"
))
```

---

### Instalar TopKAT (si no está instalada)

```r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("sarahsamorodnitsky/TopKAT")
```

O bien:

```r
library(remotes)
remotes::install_github("sarahsamorodnitsky/TopKAT")
```

---

## 📂 Datasets esperados

El archivo de **metadata** debe incluir al menos las columnas:

* `PID`
* `type`
* `FA_status`
* `sample_id`
* `x`
* `y`

El módulo **TopKAT** requiere además los archivos:

* `rips_list.rds`
* `K_dim0.rds`
* `K_dim1.rds`
* `PIDs.rds`

*(todos provenientes de la carpeta `data_selection/` descargada externamente)*

---

##  Cómo clonar y ejecutar el repositorio

### Clonar (SSH)

```bash
git clone git@github.com:USUARIO/NOMBRE_REPO.git
cd NOMBRE_REPO
```

### Ejecutar desde R (servidor)

```r
setwd("/ruta/a/tu/NOMBRE_REPO")

shiny::runApp(
  appDir = "app.R",
  host = "132.132.132.132",   # ejemplo
  port = 1111                  # ejemplo
)
```

### Ejecutar desde R (local)

```r
setwd("/ruta/a/tu/NOMBRE_REPO")

shiny::runApp("app.R")
```
