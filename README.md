# TopKAT + Connectivity Explorer

Aplicación **Shiny** que combina dos análisis para estudios de **organización espacial celular**:

---

##  Módulo **TopKAT Explorer**

Permite:

- Ejecutar **TopKAT** y **scale_importance** sobre diagramas de persistencia (PDs).
- Seleccionar automáticamente **subconjuntos de muestras** según comparaciones definidas.
- Generar carpetas de trabajo reproducibles con todos los resultados (**RDS + PNG + CSV**).
- Visualizar resultados dentro de la aplicación.

---

##  Módulo de **Conectividad Celular**

- Calcula matrices de **conectividad promedio** entre tipos celulares usando distancias espaciales *(x,y)*.
- Produce **heatmaps comparativos** entre Grupo A y Grupo B con **misma escala global**.
- Guarda las imágenes dentro del mismo *subset* generado por TopKAT.

---

#  Características principales

## ✔ TopKAT Explorer
- Carga metadata y archivos RDS (rips, kernels, PIDs).
- Permite elegir una comparación predefinida (FA vs NOFA, carcinoma vs dysplasia, etc.).
- Selecciona **nA y nB** muestras con **semillas reproducibles**.
- Genera un directorio: subset_TOPKAT_YYYYMMDD_HHMMSS/


Con:

- `rips_list_subset.rds`
- `K_dim0_subset.rds`
- `K_dim1_subset.rds`
- `PIDs_subset.rds`
- `PID_seleccionados_final.csv`
- `scale_importance_normal.png`
- `scale_importance_log10.png`

---

## ✔ Conectividad Celular
Para cada PID:

- Calcula distancias euclidianas.
- Determina conexiones (`dist ≤ ε`).
- Construye matrices **A** y **B** promedio.

Genera:

- `connectivity_groupA.svg`
- `connectivity_groupB.svg`
- `connectivity_A_vs_B.svg`

Los heatmaps muestran:

- **Eje X:** tipo celular 1  
- **Eje Y:** tipo celular 2  
- **Fill:** número de conexiones  

Todo se guarda directamente en el subset generado por TopKAT.

---


---

# 📦 Requisitos

## 🔹 Versión recomendada de R
**R ≥ 4.2**

## 🔹 Librerías necesarias
La app utiliza:

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
---

Instalación rápida
```r
install.packages(c(
  "shiny", "shinyWidgets", "shinycssloaders", "dplyr", "tidyr",
  "ggplot2", "reshape2", "stringr", "readr", "purrr", "viridis",
  "igraph", "here", "base64enc"
))
```

---
Instalar TopKAT (si no está instalada)

```r
# First, install devtools
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install from Github
devtools::install_github("sarahsamorodnitsky/TopKAT")
```

o

```r
library(remotes)
remotes::install_github("sarahsamorodnitsky/TopKAT")
```

##  Datasets esperados

El archivo de **metadata** debe incluir al menos las columnas:

- `PID`
- `type`
- `FA_status`
- `sample_id`
- `x`
- `y`

El módulo **TopKAT** requiere además los archivos:

- `rips_list.rds`
- `K_dim0.rds`
- `K_dim1.rds`
- `PIDs.rds`

---

##  Cómo clonar el repositorio

###  SSH

```bash
git clone git@github.com:USUARIO/NOMBRE_REPO.git

cd NOMBRE_REPO
```

desde R

```R

setwd("/ruta/a/tu/NOMBRE_REPO")


shiny::runApp(
  appDir = "app_v3.R",
  host = "132.132.132.132",   #example
  port = 1111   #example 
)

```
