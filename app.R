# app.R — Super-App: TopKAT Explorer + Conectividad + Conectividad Marcadores
options(shiny.useRagg = TRUE)

library(shiny)
library(shinyWidgets)
library(DT)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(stringr)
library(shinycssloaders)
library(base64enc)
library(tidyr)
library(reshape2)
library(viridis)
library(igraph)

# ============================================================
# === CONFIGURACIÓN GENERAL ===
# ============================================================
BASE_DIR <- getwd()

# Defaults solo como fallback (opcional)
DEFAULT_DATA_FILE <- NULL
DEFAULT_OUT_DIR   <- NULL
DEFAULT_PID_FILE  <- NULL
DEFAULT_OUTPUT    <- BASE_DIR


# Si usas helpers externos (recomendado), descomenta:
# source("R/helpers_io.R")
# source("R/helpers_logic.R")

# ============================================================
# === comparaciones_list ===
# ============================================================
comparaciones_list <- list(
  "carcinoma FA vs carcinoma NOFA" = list(A=list(type="carcinoma",FA="FA"), B=list(type="carcinoma",FA="NOFA")),
    
  "dysplasia FA vs dysplasia NOFA" = list(A=list(type="dysplasia",FA="FA"),B=list(type="dysplasia",FA="NOFA")),
    
  "stroma_ad_carcinoma FA vs NOFA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),B=list(type="stroma_ad_carcinoma",FA="NOFA")),
    
  "stroma_ad_dysplasia FA vs NOFA" = list(A=list(type="stroma_ad_dysplasia",FA="FA"),B=list(type="stroma_ad_dysplasia",FA="NOFA")),
    
  "carcinoma vs dysplasia (ignorar FA)" = list(A=list(type="carcinoma",FA=NULL),B=list(type="dysplasia",FA=NULL)),
    
  "carcinoma FA vs dysplasia FA" = list(A=list(type="carcinoma",FA="FA"),B=list(type="dysplasia",FA="FA")),
    
  "carcinoma NOFA vs dysplasia NOFA" = list(A=list(type="carcinoma",FA="NOFA"),B=list(type="dysplasia",FA="NOFA")),
    
  "stroma_ad_carcinoma FA vs carcinoma FA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),B=list(type="carcinoma",FA="FA")),
    
  "stroma_ad_carcinoma FA vs stroma_ad_dysplasia FA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),B=list(type="stroma_ad_dysplasia",FA="FA")),
    
  "carcinoma invasive vs carcinoma invasive (split)" = list(
    A = list(type = "carcinoma_invasive", FA = NULL),
    B = list(type = "carcinoma_invasive", FA = NULL),
    split = TRUE
  ),

  "carcinoma_in_situ vs carcinoma_in_situ (split)" =
    list(
      A = list(type = "carcinoma_in_situ", FA = NULL),
      B = list(type = "carcinoma_in_situ", FA = NULL),
      split = TRUE
    ),
    
  "carcinoma in situ vs carcinoma invasive" = list(
    A = list(type = "carcinoma_in_situ", FA = NULL),
    B = list(type = "carcinoma_invasive", FA = NULL)
  ),

  "HG dysplasia vs LG dysplasia" = list(
    A = list(type = "HG_dysplasia", FA = NULL),
    B = list(type = "LG_dysplasia", FA = NULL)
  ),

  "HG dysplasia vs HG dysplasia (split)" = list(
    A = list(type = "HG_dysplasia", FA = NULL),
    B = list(type = "HG_dysplasia", FA = NULL),
    split = TRUE
  ),

    "LG_dysplasia vs LG_dysplasia (split)" = list(
  A = list(type = "LG_dysplasia", FA = NULL),
  B = list(type = "LG_dysplasia", FA = NULL),
  split = TRUE
),

    "stroma_ad_carcinoma_invasive vs stroma_ad_carcinoma_invasive (split)" = list(
  A = list(type = "stroma_ad_carcinoma_invasive", FA = NULL),
  B = list(type = "stroma_ad_carcinoma_invasive", FA = NULL),
  split = TRUE
),

"stroma_ad_carcinoma_in_situ vs stroma_ad_carcinoma_invasive" = list(
  A = list(type = "stroma_ad_carcinoma_in_situ", FA = NULL),
  B = list(type = "stroma_ad_carcinoma_invasive", FA = NULL)
),

"stroma_ad_carcinoma_in_situ vs stroma_ad_carcinoma_in_situ (split)" = list(
  A = list(type = "stroma_ad_carcinoma_in_situ", FA = NULL),
  B = list(type = "stroma_ad_carcinoma_in_situ", FA = NULL),
  split = TRUE
),

    "stroma_ad_HG_dysplasia vs stroma_ad_HG_dysplasia (split)" = list(
  A = list(type = "stroma_ad_HG_dysplasia", FA = NULL),
  B = list(type = "stroma_ad_HG_dysplasia", FA = NULL),
  split = TRUE
),

"stroma_ad_LG_dysplasia vs stroma_ad_LG_dysplasia (split)" = list(
  A = list(type = "stroma_ad_LG_dysplasia", FA = NULL),
  B = list(type = "stroma_ad_LG_dysplasia", FA = NULL),
  split = TRUE
),


 "stroma_ad_HG_dysplasia vs stroma_ad_LG_dysplasia" =
    list(
      A = list(type = "stroma_ad_HG_dysplasia", FA = NULL),
      B = list(type = "stroma_ad_LG_dysplasia", FA = NULL)
    )

)

# ============================================================
# === FUNCIÓN: generar etiqueta legible desde cmp element ===
# ============================================================
make_label_from_cmp <- function(cmp_part) {
  if (is.null(cmp_part$FA)) {
    return(as.character(cmp_part$type))
  } else {
    return(paste0(cmp_part$type, " ", cmp_part$FA))
  }
}

# ============================================================
# === FUNCIÓN: inferir subcategoría desde nombre de archivo ===
# ============================================================
infer_subcategory_from_filename <- function(x) {
  x <- tolower(x)

  case_when(
    str_detect(x, "stroma") & str_detect(x, "carcinoma") & str_detect(x, "invasive") ~
      "stroma_ad_carcinoma_invasive",

    str_detect(x, "stroma") & str_detect(x, "carcinoma") & str_detect(x, "in_situ|insitu") ~
      "stroma_ad_carcinoma_in_situ",

    str_detect(x, "stroma") & str_detect(x, "hg") ~
      "stroma_ad_HG_dysplasia",

    str_detect(x, "stroma") & str_detect(x, "lg") ~
      "stroma_ad_LG_dysplasia",

    str_detect(x, "carcinoma") & str_detect(x, "invasive") ~
      "carcinoma_invasive",

    str_detect(x, "carcinoma") & str_detect(x, "in_situ|insitu") ~
      "carcinoma_in_situ",

    str_detect(x, "dysplasia") & str_detect(x, "hg") ~
      "HG_dysplasia",

    str_detect(x, "dysplasia") & str_detect(x, "lg") ~
      "LG_dysplasia",

    TRUE ~ "other"
  )
}

# ============================================================ 
# === UI del módulo TopKAT === 
# ============================================================
topkat_ui <- function(id) {
  ns <- NS(id)

  # ---- lógica previa (FUERA del UI) ----
  allowed_dirs <- c("data_selection", "data_circulos")

  base_dirs <- list.dirs(getwd(), recursive = FALSE, full.names = TRUE)
  base_dirs <- base_dirs[basename(base_dirs) %in% allowed_dirs]

  fluidPage(
    titlePanel("TopKAT Explorer — módulo"),
    sidebarLayout(
      sidebarPanel(
        width = 3,

        pickerInput(
          ns("base_dir"),
          "Carpeta base del análisis",
          choices = base_dirs,
          selected = if (length(base_dirs) > 0) base_dirs[1] else NULL
        ),

        uiOutput(ns("csv_selector")),
        uiOutput(ns("rds_dir_selector")),

        actionButton(ns("load_all"), "Cargar metadata y RDS"),
        hr(),
        uiOutput(ns("available_categories_ui")),
        hr(),
        selectInput(ns("cmp_choice"), "Elegir comparación", choices = names(comparaciones_list)),
        numericInput(ns("nA"), "N muestras grupo A", value = 3, min = 1),
        numericInput(ns("nB"), "N muestras grupo B", value = 3, min = 1),
        numericInput(ns("seedA"), "Seed Grupo A", value = 123),
        numericInput(ns("seedB"), "Seed Grupo B", value = 456),
        actionButton(ns("run_select"), "Seleccionar muestras y generar subset"),
        hr(),
        actionButton(ns("run_topkat"), "Ejecutar TopKAT y scale_importance"),
        hr(),
        verbatimTextOutput(ns("status"))
      ),
      mainPanel(
        width = 9,
        tabsetPanel(
          tabPanel("Resumen",
            h4("Tabla pid_info"),
            DTOutput(ns("pid_info_table")) %>% withSpinner(),
            hr(),
            h4("Conteos por categoría y FA"),
            DTOutput(ns("counts_table")) %>% withSpinner(),
            hr(),
            h4("Conteos por subcategoría y FA"),
            DTOutput(ns("subcounts_table")) %>% withSpinner()
          ),
          tabPanel("Selección",
            h4("Disponibilidad (Grupo A / Grupo B)"),
            DTOutput(ns("avail_A")),
            DTOutput(ns("avail_B")),
            hr(),
            h4("PIDs seleccionados"),
            DTOutput(ns("selected_pids_table"))
          ),
          tabPanel("Resultados",
            h4("TopKAT summary"),
            verbatimTextOutput(ns("topkat_summary")),
            hr(),
            h4("scale_importance (PNG generado)"),
            uiOutput(ns("scale_png_ui"))
          )
        )
      )
    )
  )
}

# ============================================================
# === Server del módulo TopKAT ===
# ============================================================
topkat_server <- function(id, shared_pid_df, shared_metadata, shared_outdir) {
  moduleServer(id, function(input, output, session) {
    output$sample_warning <- renderUI(NULL)
    ns <- session$ns

    rv <- reactiveValues(
      data_loaded = FALSE,
      data1 = NULL,
      pid_info = NULL,
      rips_full = NULL,
      K0_full = NULL,
      K1_full = NULL,
      PIDs_full = NULL,
      availA = NULL,
      availB = NULL,
      pid_df_final = NULL,
      subset_dir = NULL,
      topkat_res = NULL,
      res_scale_import = NULL,
      scale_png_path_normal = NULL,
      scale_png_path_log = NULL,
      A_label = NULL,
      B_label = NULL
    )

    output$csv_selector <- renderUI({
      req(input$base_dir)

      csvs <- list.files(
        input$base_dir,
        pattern = "\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )

      pickerInput(
        ns("data_file"),
        "Archivo CSV de metadata",
        choices = csvs
      )
    })

    output$rds_dir_selector <- renderUI({
      req(input$base_dir)

      dirs <- list.dirs(
        input$base_dir,
        recursive = TRUE,
        full.names = TRUE
      )

      pickerInput(
        ns("out_dir"),
        "Carpeta que contiene los RDS",
        choices = dirs
      )
    })

    output$status <- renderText("Listo (TopKAT).")

    observeEvent(input$load_all, {
      output$status <- renderText("Cargando archivos ...")
      tryCatch({
        data_file <- input$data_file
        out_dir   <- input$out_dir

        data1.df  <- safe_read_csv(data_file)
        pid_info  <- build_pid_info(data1.df, archivo_col="archivo_base")
        
        # ---- NUEVO: subcategoría derivada del nombre del archivo ----
        pid_info <- pid_info %>%
          mutate(
            subcategory = infer_subcategory_from_filename(archivo_base)
          )

        rv$data1 <- data1.df
        shared_metadata(data1.df)

        rv$pid_info <- pid_info
        rv$rips_full <- safe_read_rds(file.path(out_dir,"rips_list.rds"))
        rv$K0_full   <- safe_read_rds(file.path(out_dir,"K_dim0.rds"))
        rv$K1_full   <- safe_read_rds(file.path(out_dir,"K_dim1.rds"))
        rv$PIDs_full <- safe_read_rds(file.path(out_dir,"PIDs.rds"))
        rv$data_loaded <- TRUE

        output$status <- renderText("Carga completa")
      },
      error=function(e){ output$status <- renderText(paste("Error:", e$message)) })
    })

    output$available_categories_ui <- renderUI({
      req(rv$data_loaded)
      cats <- sort(unique(rv$pid_info$category))
      pickerInput(ns("cats"),"Filtrar categorías:",
                  choices=cats, selected=cats, multiple=TRUE,
                  options=list(`actions-box`=TRUE))
    })

    output$pid_info_table <- renderDT({
      req(rv$data_loaded)
      df <- rv$pid_info
      if (!is.null(input$cats)) df <- df %>% filter(category %in% input$cats)
      datatable(df)
    })

    output$counts_table <- renderDT({
      req(rv$data_loaded)
      df <- rv$pid_info
      tab <- as.data.frame(table(df$category, df$FA_status))
      colnames(tab) <- c("category","FA","n")
      datatable(tab)
    })

    output$subcounts_table <- renderDT({
      req(rv$data_loaded)
      df <- rv$pid_info
        
      tab <- df %>%
        group_by(subcategory, FA_status) %>%
        summarise(n = n(), .groups = "drop")
        
      datatable(tab)
    })

    observeEvent(input$run_select,{
      req(rv$data_loaded)
      output$status <- renderText("Determinando disponibilidad...")

      cmp_sel <- comparaciones_list[[input$cmp_choice]]
      pid_info <- rv$pid_info
      is_split <- isTRUE(cmp_sel$split)
      if (!is.null(input$cats)) pid_info <- pid_info %>% filter(category %in% input$cats)

      cmpA <- cmp_sel$A
      cmpB <- cmp_sel$B

      same_type <- identical(cmpA$type, cmpB$type) &&
                   identical(cmpA$FA, cmpB$FA)  

      dispA <- pid_info %>% filter(category == cmpA$type)
      if (!is.null(cmpA$FA)) dispA <- dispA %>% filter(FA_status==cmpA$FA)

      dispB <- pid_info %>% filter(category == cmpB$type)
      if (!is.null(cmpB$FA)) dispB <- dispB %>% filter(FA_status==cmpB$FA)

        # ------------------------------------------------
        # Caso especial: comparación intra-categoría
        # ------------------------------------------------
        if (is_split && same_type) {
          dispB <- dispA  }

      rv$availA <- dispA
      rv$availB <- dispB

      if (input$nA > nrow(dispA) || input$nB > nrow(dispB)) {
        output$status <- renderText("nA o nB mayor que disponibilidad.")
        return(NULL)
      }

      if (is_split && same_type) {
        total_needed <- input$nA + input$nB
        if (total_needed > nrow(dispA)) {
          output$status <- renderText("No hay suficientes muestras para split disjunto.")
          return(NULL)
        }

        selA <- select_reproducible(dispA, input$nA, input$seedA)

        remaining <- dispA %>%
          filter(!PID %in% selA$PID)

        selB <- select_reproducible(remaining, input$nB, input$seedB)
      } else {
        selA <- select_reproducible(dispA, input$nA, input$seedA)
        selB <- select_reproducible(dispB, input$nB, input$seedB)
      }

      # --- definimos etiquetas A/B desde la comparación elegida (incluye FA si aplica) ---
      A_label <- make_label_from_cmp(cmpA)
      B_label <- make_label_from_cmp(cmpB)

      # añadimos columnas 'grupo' (A/B) y 'group_name' (etiqueta legible)
      selA <- selA %>% mutate(grupo = "A", group_name = A_label)
      selB <- selB %>% mutate(grupo = "B", group_name = B_label)

      pid_df_final <- bind_rows(selA, selB) %>% mutate(PID_chr = as.character(PID))

      subset_dir <- file.path(input$out_dir,
                              paste0("subset_TOPKAT_",format(Sys.time(),"%Y%m%d_%H%M%S")))
      dir.create(subset_dir, recursive=TRUE)
      shared_outdir(subset_dir)
        
      rv$subset_dir <- subset_dir
      rv$data_file_path <- input$data_file

      selected_pids_chr <- as.character(pid_df_final$PID)

      create_subset_files(
        selected_pids_chr,
        rv$rips_full, rv$K0_full, rv$K1_full, rv$PIDs_full,
        subset_dir
      )

      # Guardamos el CSV con la columna group_name incluida
      safe_write_csv(pid_df_final, file.path(subset_dir,"PID_seleccionados_final.csv"))

      # Guardar etiquetas y subset_dir en reactiveValues para referencia
      rv$pid_df_final <- pid_df_final
      shared_pid_df(pid_df_final)    
      rv$A_label <- A_label
      rv$B_label <- B_label

      output$status <- renderText(paste("✔ Subset creado en:", subset_dir,
                                        "\nA =", A_label, " | B =", B_label))
    })

    output$avail_A <- renderDT({ req(rv$availA); datatable(rv$availA) })
    output$avail_B <- renderDT({ req(rv$availB); datatable(rv$availB) })
    output$selected_pids_table <- renderDT({
      req(rv$pid_df_final); datatable(rv$pid_df_final)
    })

    observeEvent(input$run_topkat,{
      req(rv$subset_dir)

      withProgress(message="Ejecutando análisis TopKAT...", value=0,{
        incProgress(0.10, detail="Leyendo subset...")

        rips_sub <- readRDS(file.path(rv$subset_dir,"rips_list_subset.rds"))
        K0_sub   <- readRDS(file.path(rv$subset_dir,"K_dim0_subset.rds"))
        K1_sub   <- readRDS(file.path(rv$subset_dir,"K_dim1_subset.rds"))
        PIDs_sub <- readRDS(file.path(rv$subset_dir,"PIDs_subset.rds"))

        incProgress(0.30, detail="Preparando matrices...")
        K.list <- list(K0 = K0_sub, K1 = K1_sub)

        pid_df_final <- rv$pid_df_final
        pid_df_final <- pid_df_final %>% arrange(match(as.character(PID), PIDs_sub))
        y <- ifelse(pid_df_final$grupo=="A",1,0)

        incProgress(0.60, detail="Ejecutando TopKAT...")
        library(TopKAT)

        res_topkat <- TopKAT(
          y=y,
          X=NULL,
          K.list=K.list,
          omega.list=c(0,0.5,1),
          outcome.type="binary"
        )

        rv$topkat_res <- res_topkat
        saveRDS(res_topkat,file.path(rv$subset_dir,"TopKAT_result.rds"))

        incProgress(0.80, detail="scale_importance...")

        res_scale_import <- scale_importance(
          pd.list    = rips_sub,
          y          = y,
          omega.list = c(0,0.5,1),
          threshold  = 500,
          PIDs       = seq_along(PIDs_sub),
          outcome.type = "binary"
        )

        rv$res_scale_import <- res_scale_import
        saveRDS(res_scale_import,file.path(rv$subset_dir,"scale_importance_result.rds"))

        # -------------------------------------------------------
        #   GENERAR PNGs — NORMAL Y LOG10 (se guardan en subset_dir)
        # -------------------------------------------------------
        dfp <- data.frame(
          thresh = res_scale_import$threshold.seq,
          pval   = res_scale_import$pvals
        )

        # --- PNG 1: normal ---
        png_path_normal <- file.path(rv$subset_dir, "scale_importance_normal.png")
        
        png(png_path_normal, width = 1800, height = 1200, res = 200)
        print(
          ggplot(dfp, aes(x = thresh, y = pval)) +
            geom_point() +
            theme_bw() +
            xlab(expression(epsilon)) +
            ylab(expression(p ~ "-valor")) +
            geom_vline(xintercept = res_scale_import$min.thresh,
                       linetype = "dashed") +
            geom_hline(yintercept = 0.05, color = "red", linetype = "dotted")
        )
        dev.off()
        
        # --- PNG 2: log10 ---
        png_path_log <- file.path(rv$subset_dir, "scale_importance_log10.png")
        
        png(png_path_log, width = 1800, height = 1200, res = 200)
        print(
          ggplot(dfp, aes(x = thresh, y = pval)) +
            geom_point() +
            scale_y_log10() +
            theme_bw() +
            xlab(expression(epsilon)) +
            ylab(expression(log[10](p ~ "-valor"))) +
            geom_vline(xintercept = res_scale_import$min.thresh,
                       linetype = "dashed") +
            geom_hline(yintercept = 0.05, color = "red", linetype = "dotted")
        )
        dev.off()

        rv$scale_png_path_normal <- png_path_normal
        rv$scale_png_path_log    <- png_path_log

        incProgress(1)
      })

      output$status <- renderText(
        paste0(
          "✔ TopKAT ejecutado. p-value global = ",
          signif(rv$topkat_res$overall.pval,3),
          "\nUmbral óptimo (min.thresh) = ",
          rv$res_scale_import$min.thresh
        )
      )
    })

    output$scale_png_ui <- renderUI({
      req(rv$scale_png_path_normal, rv$scale_png_path_log)

      tagList(
        h4("Plot normal"),
        tags$img(
          src = base64enc::dataURI(
            file = rv$scale_png_path_normal,
            mime = "image/png"
          ),
          width = "800px"
        ),
        hr(),

        h4("Plot con escala log10"),
        tags$img(
          src = base64enc::dataURI(
            file = rv$scale_png_path_log,
            mime = "image/png"
          ),
          width = "800px"
        )
      )
    })

    output$topkat_summary <- renderPrint({
      req(rv$topkat_res)
      cat("p-value global:", rv$topkat_res$overall.pval, "\n\n")
      print(rv$topkat_res$p.vals)
    })
  })
}

# ============================================================
# === FUNCIONES AUXILIARES ===
# ============================================================

# Función para generar conectividad (ya existe)
generate_connectivity <- function(images.df, threshold, type.column, unique.types) {
  distances <- as.matrix(dist(images.df %>% select(x, y)))
  adj <- ifelse(distances <= threshold, 1, 0)
  diag(adj) <- 0

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  igraph::vertex_attr(g, "name") <- images.df[[type.column]]

  edges <- igraph::as_edgelist(g) %>% as.data.frame()
  colnames(edges) <- c("from","to")

  if (nrow(edges) == 0) {
    return(matrix(0, length(unique.types), length(unique.types),
                  dimnames = list(unique.types, unique.types)))
  }

  tab <- edges %>%
    group_by(from,to) %>%
    summarise(count = n(), .groups="drop")

  full <- expand.grid(unique.types, unique.types, stringsAsFactors = FALSE)
  colnames(full) <- c("from","to")

  out <- full %>%
    left_join(tab, by=c("from","to")) %>%
    mutate(count = ifelse(is.na(count),0,count))

  M <- out %>% pivot_wider(names_from = to, values_from = count)
  mat <- as.matrix(M[,-1])
  rownames(mat) <- unique.types

  mat
}

# Función para plotear matriz de conectividad (ya existe)
plot_connectivity_matrix <- function(connect, title, min_val, max_val) {
  ggplot(melt(connect), aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    viridis::scale_fill_viridis(
      option = "turbo",
      limits = c(min_val, max_val)
    ) +
    theme_minimal() +
    labs(x="Tipo celular 1", y="Tipo celular 2", fill="Conexiones") +
    theme(
      axis.text.x = element_text(size=9, angle=45, hjust=1),
      axis.text.y = element_text(size=9),
      plot.title=element_text(size=14)
    ) +
    ggtitle(title)
}

# Función para crear categorías granulares CON OPCIONES ESPECÍFICAS
create_granular_categories <- function(df, cell_types, marker_choices) {
  # Crear una nueva columna con el tipo celular base
  df$granular_type <- df[[cell_types]]
  
  # Agregar marcadores seleccionados con sus opciones específicas
  for (marker_info in marker_choices) {
    marker_name <- marker_info$marker
    selected_options <- marker_info$options
    
    if (marker_name %in% colnames(df)) {
      # Determinar el status basado en el valor de la columna
      # Suponiendo que los valores son binarios (0/1) o categóricos (PD-L1+, PD-L1-)
      if (all(df[[marker_name]] %in% c(0, 1, "0", "1", TRUE, FALSE))) {
        # Es binario
        marker_status <- ifelse(as.numeric(df[[marker_name]]) > 0, 
                               paste0(toupper(marker_name), "+"), 
                               paste0(toupper(marker_name), "-"))
      } else {
        # Ya tiene formato de etiqueta (PD-L1+, PD-L1-, etc.)
        marker_status <- as.character(df[[marker_name]])
      }
      
      # Solo incluir si está en las opciones seleccionadas
      if (length(selected_options) > 0) {
        df <- df %>% 
          mutate(
            granular_type = ifelse(
              marker_status %in% selected_options,
              paste(granular_type, marker_status),
              granular_type
            )
          )
      }
    }
  }
  
  return(df$granular_type)
}

# ============================================================
# === MÓDULO: Conectividad Marcadores ===
# ============================================================

connectivity_markers_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Configuración de conectividad con marcadores"),
        
        # Selector de tipo base
        radioButtons(
          ns("base_type"),
          "Tipo celular base:",
          choices = c("Type (original)" = "type",
                     "FinalPhenotype (refinado)" = "FinalPhenotype"),
          selected = "type"
        ),
        
        # Selector de tipos celulares específicos
        uiOutput(ns("cell_type_selector")),
        
        # Selector de marcadores con opciones específicas
        h5("Marcadores y sus opciones:"),
        uiOutput(ns("marker_selectors")),
        
        hr(),
        
        # Threshold
        numericInput(ns("threshold"), "Threshold (ε)", value = 15, min = 1),
        
        # Botón para generar
        actionButton(ns("run_markers"), "Generar conectividades con marcadores"),
        
        # Advertencia de muestra
        uiOutput(ns("sample_warning_markers")),
        
        hr(),
        
        # Información sobre categorías generadas
        h5("Resumen de categorías:"),
        verbatimTextOutput(ns("category_summary"))
      ),
      
      mainPanel(
        width = 8,
        h4("Estatus del proceso"),
        verbatimTextOutput(ns("status_markers")),
        
        # Solo dos pestañas: Matrices y Normalizadas
        tabsetPanel(
          tabPanel("Matrices",
            h4("Conectividad Grupo A"),
            imageOutput(ns("imgA_markers"), height = "auto"),
            h4("Conectividad Grupo B"),
            imageOutput(ns("imgB_markers"), height = "auto"),
            h4("Comparación A vs B"),
            imageOutput(ns("imgAB_markers"), height = "auto")
          ),
          
          tabPanel("Normalizadas",
            h4("Conectividad normalizada — Grupo A"),
            imageOutput(ns("imgA_row_markers"), height = "auto"),
            h4("Conectividad normalizada — Grupo B"),
            imageOutput(ns("imgB_row_markers"), height = "auto")
          )
        )
      )
    )
  )
}

connectivity_markers_server <- function(id, shared_pid_df, shared_metadata, shared_outdir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Estado inicial
    output$status_markers <- renderText("Listo para iniciar (Conectividad Marcadores).")
    output$sample_warning_markers <- renderUI(NULL)
    
    # Selector dinámico de tipos celulares
    output$cell_type_selector <- renderUI({
      req(shared_metadata())
      df <- shared_metadata()
      base_type <- input$base_type
      
      if (base_type %in% colnames(df)) {
        cell_types <- sort(unique(df[[base_type]]))
        pickerInput(
          ns("selected_cell_types"),
          "Tipos celulares a incluir:",
          choices = cell_types,
          selected = cell_types,
          multiple = TRUE,
          options = list(`actions-box` = TRUE, `live-search` = TRUE)
        )
      } else {
        selectInput(ns("selected_cell_types"), 
                   "Tipos celulares a incluir:",
                   choices = "Columna no disponible",
                   multiple = TRUE)
      }
    })
    
    # Selector dinámico de marcadores con sus opciones
    output$marker_selectors <- renderUI({
      req(shared_metadata())
      df <- shared_metadata()
      
      # Lista de posibles marcadores
      possible_markers <- c("tim3", "pd_l1", "pd_1", "p_s6", "ki67")
      
      # Filtrar los que existen en los datos
      existing_markers <- possible_markers[possible_markers %in% colnames(df)]
      
      if (length(existing_markers) == 0) {
        return(tags$p("No se encontraron marcadores en los datos."))
      }
      
      # Crear un selector para cada marcador
      tagList(
        lapply(existing_markers, function(marker) {
          # Obtener valores únicos para este marcador
          unique_values <- sort(unique(df[[marker]]))
          
          # Convertir a opciones legibles
          options_list <- list()
          for (val in unique_values) {
            if (is.numeric(val) || val %in% c("0", "1", "TRUE", "FALSE")) {
              # Es binario
              pos_label <- ifelse(marker == "pd_l1", "PD-L1+", paste0(toupper(marker), "+"))
              neg_label <- ifelse(marker == "pd_l1", "PD-L1-", paste0(toupper(marker), "-"))
              
              if (as.numeric(val) > 0) {
                options_list[[pos_label]] <- pos_label
              } else {
                options_list[[neg_label]] <- neg_label
              }
            } else {
              # Ya tiene formato de etiqueta
              options_list[[as.character(val)]] <- as.character(val)
            }
          }
          
          # Crear checkbox group para este marcador
          checkboxGroupInput(
            ns(paste0("marker_", marker)),
            label = ifelse(marker == "pd_l1", "PD-L1", toupper(marker)),
            choices = unique(unlist(options_list)),
            selected = character(0)
          )
        })
      )
    })
    
    # Observar cuando se ejecuta el análisis con marcadores
    observeEvent(input$run_markers, {
      req(shared_pid_df(), shared_metadata())
      
      withProgress(message = "Procesando conectividad con marcadores...", value = 0, {
        
        # 1) Obtener datos
        pid_df <- shared_pid_df()
        df <- shared_metadata()
        
        if (is.null(pid_df) || is.null(df)) {
          output$status_markers <- renderText("ERROR: No hay datos disponibles.")
          return()
        }
        
        # 2) Validar selecciones
        base_type <- input$base_type
        selected_types <- input$selected_cell_types
        
        if (is.null(selected_types) || length(selected_types) == 0) {
          output$status_markers <- renderText("ERROR: Selecciona al menos un tipo celular.")
          return()
        }
        
        if (!base_type %in% colnames(df)) {
          output$status_markers <- renderText(paste("ERROR: La columna", base_type, "no existe."))
          return()
        }
        
        # 3) Obtener selecciones de marcadores
        possible_markers <- c("tim3", "pd_l1", "pd_1", "p_s6", "ki67")
        marker_choices <- list()
        
        for (marker in possible_markers) {
          input_id <- paste0("marker_", marker)
          selected_options <- input[[input_id]]
          
          if (!is.null(selected_options) && length(selected_options) > 0) {
            marker_choices[[length(marker_choices) + 1]] <- list(
              marker = marker,
              options = selected_options
            )
          }
        }
        
        incProgress(0.2, detail = "Creando categorías con marcadores...")
        
        # 4) Filtrar datos
        df_filtered <- df %>%
          filter(PID %in% pid_df$PID) %>%
          filter(!!sym(base_type) %in% selected_types)
        
        if (nrow(df_filtered) == 0) {
          output$status_markers <- renderText("ERROR: No hay datos después del filtrado.")
          return()
        }
        
        # 5) Crear categorías con marcadores
        if (length(marker_choices) > 0) {
          df_filtered$granular_category <- create_granular_categories(
            df_filtered, 
            base_type, 
            marker_choices
          )
        } else {
          # Si no hay marcadores seleccionados, usar solo el tipo base
          df_filtered$granular_category <- df_filtered[[base_type]]
        }
        
        # Filtrar filas donde granular_category no sea NA
        df_filtered <- df_filtered %>% filter(!is.na(granular_category))
        
        # Mostrar resumen de categorías
        category_summary <- df_filtered %>%
          group_by(granular_category) %>%
          summarise(n_cells = n(), .groups = "drop") %>%
          arrange(desc(n_cells))
        
        output$category_summary <- renderPrint({
          cat("Total de categorías generadas:", nrow(category_summary), "\n")
          cat("Total de células:", sum(category_summary$n_cells), "\n\n")
          cat("Top 10 categorías:\n")
          print(head(category_summary, 10))
        })
        
        incProgress(0.4, detail = "Calculando matrices de conectividad...")
        
        # 6) Calcular matrices por grupo
        unique_categories <- sort(unique(df_filtered$granular_category))
        
        if (length(unique_categories) < 2) {
          output$status_markers <- renderText("ERROR: Se necesitan al menos 2 categorías para calcular conectividad.")
          return()
        }
        
        # Inicializar matrices
        M_A <- matrix(0, length(unique_categories), length(unique_categories),
                     dimnames = list(unique_categories, unique_categories))
        M_B <- M_A
        
        # Obtener información de grupos
        if ("grupo" %in% colnames(pid_df)) {
          # Filtrar por grupo A
          pids_A <- pid_df %>% filter(grupo == "A") %>% pull(PID)
          df_A <- df_filtered %>% filter(PID %in% pids_A)
          
          # Filtrar por grupo B
          pids_B <- pid_df %>% filter(grupo == "B") %>% pull(PID)
          df_B <- df_filtered %>% filter(PID %in% pids_B)
          
          # Calcular matrices
          if (nrow(df_A) > 0) {
            for (pid in unique(df_A$PID)) {
              sub <- df_A %>% filter(PID == pid)
              if (nrow(sub) > 0) {
                mat <- generate_connectivity(sub, input$threshold, "granular_category", unique_categories)
                M_A <- M_A + mat
              }
            }
          }
          
          if (nrow(df_B) > 0) {
            for (pid in unique(df_B$PID)) {
              sub <- df_B %>% filter(PID == pid)
              if (nrow(sub) > 0) {
                mat <- generate_connectivity(sub, input$threshold, "granular_category", unique_categories)
                M_B <- M_B + mat
              }
            }
          }
          
          # Normalizar por número de muestras
          nA <- length(unique(df_A$PID))
          nB <- length(unique(df_B$PID))
          
          M_A <- if (nA > 0) M_A / nA else M_A
          M_B <- if (nB > 0) M_B / nB else M_B
          
          incProgress(0.6, detail = "Generando visualizaciones...")
          
          # 7) Generar gráficos
          min_global <- min(c(M_A, M_B), na.rm = TRUE)
          max_global <- max(c(M_A, M_B), na.rm = TRUE)
          
          # Obtener nombres de grupo
          if ("group_name" %in% colnames(pid_df)) {
            grupo_A_name <- unique(pid_df$group_name[pid_df$grupo == "A"])[1]
            grupo_B_name <- unique(pid_df$group_name[pid_df$grupo == "B"])[1]
          } else {
            grupo_A_name <- "A"
            grupo_B_name <- "B"
          }
          
          # Crear plots
          pA <- plot_connectivity_matrix(
            M_A,
            paste0("Conectividad con Marcadores — ", grupo_A_name, " (n=", nA, ")"),
            min_global, max_global
          )
          
          pB <- plot_connectivity_matrix(
            M_B,
            paste0("Conectividad con Marcadores — ", grupo_B_name, " (n=", nB, ")"),
            min_global, max_global
          )
          
          # Versión normalizada por fila (0-100%)
          plot_row_scaled <- function(mat, title) {
            if (nrow(mat) == 0 || ncol(mat) == 0) {
              return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No hay datos")) + theme_void())
            }
            
            df_melt <- melt(mat)
            colnames(df_melt) <- c("row", "col", "value")
            df_melt$value <- as.numeric(as.character(df_melt$value))
            
            df_melt <- df_melt %>%
              group_by(row) %>%
              mutate(
                max_row = max(value, na.rm = TRUE),
                value_norm = ifelse(max_row > 0, value / max_row * 100, 0)
              ) %>%
              ungroup()
            
            ggplot(df_melt, aes(x = col, y = row, fill = value_norm)) +
              geom_tile(color = "white") +
              geom_text(aes(label = round(value_norm, 1)), size = 2.5) +
              scale_fill_viridis(option = "turbo", limits = c(0, 100), 
                               name = "% del máximo\npor fila") +
              theme_minimal() +
              labs(x = "Tipo celular destino", y = "Tipo celular origen", title = title) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                axis.text.y = element_text(size = 8),
                plot.title = element_text(size = 12)
              )
          }
          
          pA_row <- plot_row_scaled(M_A, paste0("Normalizado — ", grupo_A_name))
          pB_row <- plot_row_scaled(M_B, paste0("Normalizado — ", grupo_B_name))
          
          # Combinar A y B lado a lado
          if (requireNamespace("patchwork", quietly = TRUE) && nrow(M_A) > 0 && nrow(M_B) > 0) {
            pAB <- pA + pB
          } else {
            pAB <- pA
          }
          
          incProgress(0.8, detail = "Guardando resultados...")
          
          # 8) Guardar resultados
          pid_dir <- shared_outdir()
          if (!is.null(pid_dir)) {
            # Crear nombres de archivo
            type_tag <- make.names(base_type)
            eps_tag <- paste0("eps", input$threshold)
            
            # Crear tag para marcadores seleccionados
            markers_tag <- ""
            if (length(marker_choices) > 0) {
              marker_names <- sapply(marker_choices, function(x) x$marker)
              markers_tag <- paste0("_", paste(marker_names, collapse = "_"))
            }
            
            outA <- file.path(pid_dir, paste0("markers_A", markers_tag, "_", eps_tag, ".svg"))
            outB <- file.path(pid_dir, paste0("markers_B", markers_tag, "_", eps_tag, ".svg"))
            outAB <- file.path(pid_dir, paste0("markers_AB", markers_tag, "_", eps_tag, ".svg"))
            outA_row <- file.path(pid_dir, paste0("markers_A_row", markers_tag, "_", eps_tag, ".svg"))
            outB_row <- file.path(pid_dir, paste0("markers_B_row", markers_tag, "_", eps_tag, ".svg"))
            
            # Guardar SVG
            if (nrow(M_A) > 0) ggsave(outA, pA, width = 10, height = 8)
            if (nrow(M_B) > 0) ggsave(outB, pB, width = 10, height = 8)
            ggsave(outAB, pAB, width = 16, height = 8)
            if (nrow(M_A) > 0) ggsave(outA_row, pA_row, width = 10, height = 8)
            if (nrow(M_B) > 0) ggsave(outB_row, pB_row, width = 10, height = 8)
            
            # Guardar datos (no es neceseario)
            # if (nrow(M_A) > 0) saveRDS(M_A, file.path(pid_dir, paste0("markers_matrix_A", markers_tag, "_", eps_tag, ".rds")))
            # if (nrow(M_B) > 0) saveRDS(M_B, file.path(pid_dir, paste0("markers_matrix_B", markers_tag, "_", eps_tag, ".rds")))
            
            # Mostrar imágenes
            if (nrow(M_A) > 0) {
              output$imgA_markers <- renderImage({ 
                list(src = outA, contentType = "image/svg+xml", width = "100%") 
              }, deleteFile = FALSE)
            }
            
            if (nrow(M_B) > 0) {
              output$imgB_markers <- renderImage({ 
                list(src = outB, contentType = "image/svg+xml", width = "100%") 
              }, deleteFile = FALSE)
            }
            
            output$imgAB_markers <- renderImage({ 
              list(src = outAB, contentType = "image/svg+xml", width = "100%") 
            }, deleteFile = FALSE)
            
            if (nrow(M_A) > 0) {
              output$imgA_row_markers <- renderImage({ 
                list(src = outA_row, contentType = "image/svg+xml", width = "100%") 
              }, deleteFile = FALSE)
            }
            
            if (nrow(M_B) > 0) {
              output$imgB_row_markers <- renderImage({ 
                list(src = outB_row, contentType = "image/svg+xml", width = "100%") 
              }, deleteFile = FALSE)
            }
          }
          
          # Advertencia de muestra
          n_total <- nrow(pid_df)
          if (n_total < 50) {
            output$sample_warning_markers <- renderUI({
              tags$div(
                style = "color:#b30000; font-weight:bold; margin-top: 10px;",
                paste0(
                  "Advertencia: se seleccionaron ", n_total, 
                  " muestras. Para resultados confiables se recomiendan al menos 50."
                )
              )
            })
          } else {
            output$sample_warning_markers <- renderUI(NULL)
          }
          
          output$status_markers <- renderText(
            paste0("¡Listo! Generadas ", length(unique_categories), 
                   " categorías con marcadores. Matrices guardadas en: ", pid_dir)
          )
        } else {
          output$status_markers <- renderText("ERROR: No hay información de grupos en los datos.")
        }
        
        incProgress(1)
      })
    })
  })
}

# ============================================================
# === MÓDULO: Conectividad Promedio (ORIGINAL - SIN CAMBIOS) ===
# ============================================================

# Función para heatmap normalizado por fila (0–100%)
plot_connectivity_matrix_row_scaled <- function(mat, title) {
  df <- reshape2::melt(mat)
  colnames(df) <- c("row", "col", "value")
  
  # ---- FIX CRÍTICO: forzar valor numérico ----
  df$value <- as.numeric(as.character(df$value))
  
  # normalización por fila (0–100 %)
  df <- df %>%
    group_by(row) %>%
    mutate(
      max_row = max(value, na.rm = TRUE),
      value_norm = ifelse(max_row > 0, value / max_row * 100, 0)
    ) %>%
    ungroup()
  
  n_cols <- length(unique(df$col))
  
  df_max <- df %>%
    distinct(row, max_row) %>%
    mutate(x_pos = n_cols + 0.6)
  
  ggplot(df, aes(x = col, y = row)) +
    geom_tile(aes(fill = value_norm), color = "white") +
    geom_text(
      aes(label = round(value_norm, 1)),
      size = 3
    ) +
    geom_text(
      data = df_max,
      aes(
        x = x_pos,
        y = row,
        label = paste0("", round(max_row, 1))
      ),
      hjust = 0,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    scale_fill_viridis(
      option = "turbo",
      limits = c(0, 100),
      name = "% del máximo\npor fila"
    ) +
    scale_x_discrete(expand = expansion(add = c(0, 2))) +
    theme_minimal() +
    labs(
      x = "Tipo celular destino",
      y = "Phenotype origen",
      title = title
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14)
    )
}

connectivity_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("threshold"),"Threshold (ε)", value=15, min=1),
        selectInput(
          ns("phenotype_col"),
          "Clasificación celular a usar:",
          choices = c(
            "Type (original)" = "type",
            "FinalPhenotype (refinado)" = "FinalPhenotype"
          ),
          selected = "type"
        ),
        actionButton(ns("run"),"Generar conectividades"),
        uiOutput(ns("sample_warning"))
      ),
      mainPanel(
        h4("Estatus del proceso"),
        verbatimTextOutput(ns("status")),
        h4("Conectividad Grupo A"),
        imageOutput(ns("imgA"), height = "auto"),
        h4("Conectividad Grupo B"),
        imageOutput(ns("imgB"), height = "auto"),
        h4("Comparación A vs B"),
        imageOutput(ns("imgAB"), height = "auto"),
        h4("Conectividad normalizada por phenotype — Grupo A"),
        imageOutput(ns("imgA_row"), height = "auto"),
        h4("Conectividad normalizada por phenotype — Grupo B"),
        imageOutput(ns("imgB_row"), height = "auto")
      )
    )
  )
}

connectivity_server <- function(id, shared_pid_df, shared_metadata, shared_outdir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Estado inicial
    output$status <- renderText("Listo para iniciar (Connectivity).")
    output$sample_warning <- renderUI(NULL)

    observeEvent(input$run, {
      withProgress(message = "Procesando...", value = 0, {
        # ------------------------------------------
        # 1) Obtener selección de PIDs desde TopKAT (en memoria)
        # ------------------------------------------
        pid_df <- shared_pid_df()

        if (is.null(pid_df)) {
          output$status <- renderText("ERROR: No hay selección de PIDs desde TopKAT.")
          return()
        }

        if (!"grupo" %in% colnames(pid_df)) {
          output$status <- renderText("ERROR: La selección no tiene columna 'grupo'.")
          return()
        }

        # ------------------------------------------
        # Directorio de salida (subset_TOPKAT_*)
        # ------------------------------------------
        pid_dir <- shared_outdir()

        if (is.null(pid_dir)) {
          output$status <- renderText(
            "ERROR: No se encontró la carpeta de salida de TopKAT."
          )
          return()
        }

        incProgress(0.15)

        # ------------------------------------------
        # 2) Obtener nombres reales de los grupos desde 'group_name' si existe
        # ------------------------------------------
        if ("group_name" %in% colnames(pid_df)) {
          grupo_A_name <- unique(pid_df$group_name[pid_df$grupo == "A"])[1]
          grupo_B_name <- unique(pid_df$group_name[pid_df$grupo == "B"])[1]
          if (is.na(grupo_A_name) || is.null(grupo_A_name)) grupo_A_name <- "A"
          if (is.na(grupo_B_name) || is.null(grupo_B_name)) grupo_B_name <- "B"
        } else {
          grupos <- unique(pid_df$grupo)
          grupo_A_name <- ifelse(length(grupos) >= 1, grupos[1], "A")
          grupo_B_name <- ifelse(length(grupos) >= 2, grupos[2], "B")
        }

        output$status <- renderText(paste0("Usando nombres de grupo: A = ", grupo_A_name,
                                           " | B = ", grupo_B_name))

        incProgress(0.05)

        # ------------------------------------------
        # 3) Filtrar datos y calcular matrices promedio
        # ------------------------------------------

        df <- shared_metadata()
        
        if (is.null(df)) {
          output$status <- renderText("ERROR: No hay metadata cargada desde TopKAT.")
          return()
        }

        df <- df %>% filter(PID %in% pid_df$PID)
        type_col <- input$phenotype_col
        type_tag <- make.names(type_col)            # type / FinalPhenotype
        eps_tag  <- paste0("eps", input$threshold)  # eps15, eps30, etc.

        if (!type_col %in% colnames(df)) {
          output$status <- renderText(
            paste("ERROR: La columna", type_col, "no existe en el CSV")
          )
          return()
        }

        cell.types <- sort(unique(df[[type_col]]))

        M_A <- matrix(0, length(cell.types), length(cell.types), dimnames=list(cell.types, cell.types))
        M_B <- M_A

        incProgress(0.25)

        output$status <- renderText("Generando conectividades por PID...")

        for (i in seq_len(nrow(pid_df))) {
          pid <- pid_df$PID[i]
          grp <- pid_df$grupo[i]

          sub <- df %>% filter(PID == pid)
          if (nrow(sub) == 0) next

          mat <- generate_connectivity(sub, input$threshold, type_col, cell.types)

          if (grp == "A") M_A <- M_A + mat
          if (grp == "B") M_B <- M_B + mat
        }

        nA <- sum(pid_df$grupo == "A")
        nB <- sum(pid_df$grupo == "B")

        M_A <- M_A / max(nA,1)
        M_B <- M_B / max(nB,1)

        n_total <- nA + nB

        if (n_total < 50) {
          output$sample_warning <- renderUI({
            tags$div(
              style = "color:#b30000; font-weight:bold;",
              paste0(
                "Advertencia estadística: se seleccionaron ",
                n_total,
                " muestras (A=", nA, ", B=", nB, "). ",
                "Para resultados confiables se recomiendan al menos 50 muestras."
              )
            )
          })
        } else {
          output$sample_warning <- renderUI(NULL)
        }

        incProgress(0.25)

        # ------------------------------------------
        # 4) Graficar con la misma escala y guardar en pid_dir (subset_dir)
        # ------------------------------------------
        min_global <- min(M_A, M_B, na.rm = TRUE)
        max_global <- max(M_A, M_B, na.rm = TRUE)

        output$status <- renderText("Generando gráficas...")

        pA  <- plot_connectivity_matrix(
          M_A,
          paste0("Promedio — ", grupo_A_name, " (n=", nA, ")"),
          min_global, max_global
        )

        pB  <- plot_connectivity_matrix(
          M_B,
          paste0("Promedio — ", grupo_B_name, " (n=", nB, ")"),
          min_global, max_global
        )

        # ============================================================
        # === NUEVO: versiones normalizadas por fila (0–100%) ===
        # ============================================================
        pA_row <- plot_connectivity_matrix_row_scaled(
          M_A,
          paste0("Normalizado por phenotype — ", grupo_A_name)
        )
        
        pB_row <- plot_connectivity_matrix_row_scaled(
          M_B,
          paste0("Normalizado por phenotype — ", grupo_B_name)
        )

        # combinamos lado a lado (si patchwork está disponible)
        if (requireNamespace("patchwork", quietly = TRUE)) {
          pAB <- pA + pB
        } else {
          pAB <- pA
        }

        incProgress(0.05)
        output$status <- renderText("Guardando imágenes en carpeta...")

        # Guardar con nombres basados en etiquetas reales y en la carpeta del subset
        outA <- file.path(pid_dir,paste0("connectivity_",make.names(grupo_A_name), "_",type_tag, "_",eps_tag,".svg"))
        outB <- file.path(pid_dir,paste0("connectivity_",make.names(grupo_B_name), "_",type_tag, "_",eps_tag,".svg"))
        outAB <- file.path(pid_dir,paste0("connectivity_",make.names(grupo_A_name), "_vs_",make.names(grupo_B_name), "_",type_tag, "_",eps_tag,".svg"))

        outA_row <- file.path(pid_dir,paste0("connectivity_row_scaled_",make.names(grupo_A_name), "_",type_tag,"_",eps_tag,".svg"))
        outB_row <- file.path(pid_dir,paste0("connectivity_row_scaled_",make.names(grupo_B_name), "_",type_tag,"_",eps_tag,".svg"))

        ggsave(outA_row, pA_row, width = 9, height = 7)
        ggsave(outB_row, pB_row, width = 9, height = 7)

        ggsave(outA, pA, width = 8, height = 6)
        ggsave(outB, pB, width = 8, height = 6)
        ggsave(outAB, pAB, width = 12, height = 6)

        # ------------------------------------------
        # 5) Mostrar imágenes en la app
        # ------------------------------------------
        output$imgA  <- renderImage({ list(src = outA,  contentType = "image/svg+xml") }, deleteFile = FALSE)
        output$imgB  <- renderImage({ list(src = outB,  contentType = "image/svg+xml") }, deleteFile = FALSE)
        output$imgAB <- renderImage({ list(src = outAB, contentType = "image/svg+xml") }, deleteFile = FALSE)
        output$imgA_row <- renderImage({list(src = outA_row, contentType = "image/svg+xml")}, deleteFile = FALSE)
        output$imgB_row <- renderImage({list(src = outB_row, contentType = "image/svg+xml")}, deleteFile = FALSE)

        output$status <- renderText(paste0("¡Listo! Imágenes guardadas en: ", pid_dir))

        incProgress(0.05)
      })
    })
  })
}

# ============================================================
# === UI PRINCIPAL: NAVBAR con TRES pestañas ===
# ============================================================
ui <- navbarPage(
  "Super-App: TopKAT + Connectivity",
  tabPanel("TopKAT Explorer", topkat_ui("topkat")),
  tabPanel("Conectividad Promedio", connectivity_ui("connect")),
  tabPanel("Conectividad Marcadores", connectivity_markers_ui("markers"))
)

# ============================================================
# === SERVER PRINCIPAL: Arranca los módulos ===
# ============================================================
server <- function(input, output, session) {
  # Estado compartido entre módulos
  shared_pid_df   <- reactiveVal(NULL)
  shared_metadata <- reactiveVal(NULL)
  shared_outdir <- reactiveVal(NULL)

  # Inicializar módulos
  topkat_server("topkat", shared_pid_df, shared_metadata, shared_outdir)
  connectivity_server("connect", shared_pid_df, shared_metadata, shared_outdir)
  connectivity_markers_server("markers", shared_pid_df, shared_metadata, shared_outdir)
}

shinyApp(ui, server)