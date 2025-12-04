# app.R — Super-App: TopKAT Explorer + Conectividad (todo en un archivo)
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
BASE_DIR  <- getwd()
DATA_DIR  <- file.path(BASE_DIR, "data")

DEFAULT_DATA_FILE <- file.path(DATA_DIR, "datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv")
DEFAULT_OUT_DIR   <- DATA_DIR
DEFAULT_PID_FILE  <- file.path(DATA_DIR, "PID_seleccionados_final.csv")
DEFAULT_OUTPUT    <- DATA_DIR

# Si usas helpers externos (recomendado), descomenta:
# source("R/helpers_io.R")
# source("R/helpers_logic.R")

# ============================================================
# === comparaciones_list ===
# ============================================================
comparaciones_list <- list(
  "carcinoma FA vs carcinoma NOFA" = list(A=list(type="carcinoma",FA="FA"),
                                          B=list(type="carcinoma",FA="NOFA")),
  "dysplasia FA vs dysplasia NOFA" = list(A=list(type="dysplasia",FA="FA"),
                                          B=list(type="dysplasia",FA="NOFA")),
  "stroma_ad_carcinoma FA vs NOFA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),
                                          B=list(type="stroma_ad_carcinoma",FA="NOFA")),
  "stroma_ad_dysplasia FA vs NOFA" = list(A=list(type="stroma_ad_dysplasia",FA="FA"),
                                          B=list(type="stroma_ad_dysplasia",FA="NOFA")),
  "carcinoma vs dysplasia (ignorar FA)" = list(A=list(type="carcinoma",FA=NULL),
                                               B=list(type="dysplasia",FA=NULL)),
  "carcinoma FA vs dysplasia FA" = list(A=list(type="carcinoma",FA="FA"),
                                        B=list(type="dysplasia",FA="FA")),
  "carcinoma NOFA vs dysplasia NOFA" = list(A=list(type="carcinoma",FA="NOFA"),
                                            B=list(type="dysplasia",FA="NOFA")),
  "stroma_ad_carcinoma FA vs carcinoma FA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),
                                                  B=list(type="carcinoma",FA="FA")),
  "stroma_ad_carcinoma FA vs stroma_ad_dysplasia FA" = list(A=list(type="stroma_ad_carcinoma",FA="FA"),
                                                            B=list(type="stroma_ad_dysplasia",FA="FA"))
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
# === UI del módulo TopKAT ===
# ============================================================
topkat_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("TopKAT Explorer — módulo"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        textInput(ns("data_file"), "Ruta CSV metadata", value = DEFAULT_DATA_FILE),
        textInput(ns("out_dir"), "Directorio con RDS", value = DEFAULT_OUT_DIR),
        actionButton(ns("load_all"), "📥 Cargar metadata y RDS"),
        hr(),
        uiOutput(ns("available_categories_ui")),
        hr(),
        selectInput(ns("cmp_choice"), "Elegir comparación", choices = names(comparaciones_list)),
        numericInput(ns("nA"), "N muestras grupo A", value = 3, min = 1),
        numericInput(ns("nB"), "N muestras grupo B", value = 3, min = 1),
        numericInput(ns("seedA"), "Seed Grupo A", value = 123),
        numericInput(ns("seedB"), "Seed Grupo B", value = 456),
        actionButton(ns("run_select"), "🎯 Seleccionar muestras y generar subset"),
        hr(),
        actionButton(ns("run_topkat"), "⚙️ Ejecutar TopKAT y scale_importance"),
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
                   DTOutput(ns("counts_table")) %>% withSpinner()
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
topkat_server <- function(id) {
  moduleServer(id, function(input, output, session) {
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

    output$status <- renderText("Listo (TopKAT).")

    observeEvent(input$load_all, {
      output$status <- renderText("Cargando archivos ...")
      tryCatch({
        data_file <- input$data_file
        out_dir   <- input$out_dir

        data1.df  <- safe_read_csv(data_file)
        pid_info  <- build_pid_info(data1.df, archivo_col="archivo_base")

        rv$data1 <- data1.df
        rv$pid_info <- pid_info
        rv$rips_full <- safe_read_rds(file.path(out_dir,"rips_list.rds"))
        rv$K0_full   <- safe_read_rds(file.path(out_dir,"K_dim0.rds"))
        rv$K1_full   <- safe_read_rds(file.path(out_dir,"K_dim1.rds"))
        rv$PIDs_full <- safe_read_rds(file.path(out_dir,"PIDs.rds"))
        rv$data_loaded <- TRUE

        output$status <- renderText("✔ Carga completa")
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

    observeEvent(input$run_select,{
      req(rv$data_loaded)
      output$status <- renderText("Determinando disponibilidad...")

      cmp_sel <- comparaciones_list[[input$cmp_choice]]
      pid_info <- rv$pid_info
      if (!is.null(input$cats)) pid_info <- pid_info %>% filter(category %in% input$cats)

      cmpA <- cmp_sel$A
      cmpB <- cmp_sel$B

      dispA <- pid_info %>% filter(category == cmpA$type)
      if (!is.null(cmpA$FA)) dispA <- dispA %>% filter(FA_status==cmpA$FA)

      dispB <- pid_info %>% filter(category == cmpB$type)
      if (!is.null(cmpB$FA)) dispB <- dispB %>% filter(FA_status==cmpB$FA)

      rv$availA <- dispA
      rv$availB <- dispB

      if (input$nA > nrow(dispA) || input$nB > nrow(dispB)) {
        output$status <- renderText("⚠️ nA o nB mayor que disponibilidad.")
        return(NULL)
      }

      selA <- select_reproducible(dispA, input$nA, input$seedA)
      selB <- select_reproducible(dispB, input$nB, input$seedB)

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
      rv$subset_dir <- subset_dir
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
                       linetype = "dashed")
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
                       linetype = "dashed")
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
# === DEFINICIÓN: MÓDULO Conectividad (App 2) ===
# ============================================================
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

connectivity_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("threshold"),"Threshold (ε)", value=15, min=1),
        textInput(ns("outdir"),"Directorio donde buscar subsets (opcional)", value=DEFAULT_OUTPUT),
        actionButton(ns("run"),"Generar conectividades")
      ),
      mainPanel(
        h4("Estatus del proceso"),
        verbatimTextOutput(ns("status")),
        h4("Conectividad Grupo A"),
        imageOutput(ns("imgA"), height = "auto"),
        h4("Conectividad Grupo B"),
        imageOutput(ns("imgB"), height = "auto"),
        h4("Comparación A vs B"),
        imageOutput(ns("imgAB"), height = "auto")
      )
    )
  )
}

connectivity_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$status <- renderText("Listo para iniciar (Connectivity).")

    observeEvent(input$run, {

      withProgress(message="Procesando...", value = 0, {

        # ------------------------------------------
        # 1) Buscar el PID_seleccionados_final.csv más reciente
        # ------------------------------------------
        output$status <- renderText("Buscando PID_seleccionados_final.csv (último)...")

        search_dir <- input$outdir
        if (is.null(search_dir) || search_dir == "") search_dir <- DATA_DIR

        matches <- list.files(path = search_dir, pattern = "PID_seleccionados_final\\.csv$",
                              recursive = TRUE, full.names = TRUE)

        # fallback al DEFAULT_PID_FILE si no hay matches
        if (length(matches) == 0 && file.exists(DEFAULT_PID_FILE)) {
          pid_file <- DEFAULT_PID_FILE
        } else if (length(matches) == 0) {
          output$status <- renderText("ERROR: No se encontró ningún PID_seleccionados_final.csv.")
          return()
        } else {
          info <- file.info(matches)
          pid_file <- rownames(info[which.max(info$mtime), , drop = FALSE])
        }

        output$status <- renderText(paste("Usando PID CSV:", pid_file))
        pid_df <- read.csv(pid_file, stringsAsFactors = FALSE)

        incProgress(0.15)

        if (!"grupo" %in% colnames(pid_df)) {
          output$status <- renderText("ERROR: El archivo PID no tiene columna 'grupo'")
          return()
        }

        # Carpeta donde se guardarán las imágenes -> carpeta del PID encontrado
        pid_dir <- dirname(pid_file)

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
        if (!file.exists(DEFAULT_DATA_FILE)) {
          output$status <- renderText("ERROR: No se encontró datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv")
          return()
        }

        df <- read.csv(DEFAULT_DATA_FILE, stringsAsFactors = FALSE)
        df <- df %>% filter(PID %in% pid_df$PID)
        cell.types <- sort(unique(df$type))

        M_A <- matrix(0, length(cell.types), length(cell.types), dimnames=list(cell.types, cell.types))
        M_B <- M_A

        incProgress(0.25)

        output$status <- renderText("Generando conectividades por PID...")

        for (i in seq_len(nrow(pid_df))) {
          pid <- pid_df$PID[i]
          grp <- pid_df$grupo[i]

          sub <- df %>% filter(PID == pid)
          if (nrow(sub) == 0) next

          mat <- generate_connectivity(sub, input$threshold, "type", cell.types)

          if (grp == "A") M_A <- M_A + mat
          if (grp == "B") M_B <- M_B + mat
        }

        nA <- sum(pid_df$grupo == "A")
        nB <- sum(pid_df$grupo == "B")

        M_A <- M_A / max(nA,1)
        M_B <- M_B / max(nB,1)

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

        # combinamos lado a lado (si patchwork está disponible)
        if (requireNamespace("patchwork", quietly = TRUE)) {
          pAB <- pA + pB
        } else {
          pAB <- pA
        }

        incProgress(0.05)
        output$status <- renderText("Guardando imágenes en carpeta...")

        # Guardar con nombres basados en etiquetas reales y en la carpeta del subset
        outA  <- file.path(pid_dir, paste0("connectivity_", make.names(grupo_A_name), ".svg"))
        outB  <- file.path(pid_dir, paste0("connectivity_", make.names(grupo_B_name), ".svg"))
        outAB <- file.path(pid_dir, paste0("connectivity_", make.names(grupo_A_name), "_vs_", make.names(grupo_B_name), ".svg"))

        ggsave(outA, pA, width = 8, height = 6)
        ggsave(outB, pB, width = 8, height = 6)
        ggsave(outAB, pAB, width = 12, height = 6)

        # ------------------------------------------
        # 5) Mostrar imágenes en la app
        # ------------------------------------------
        output$imgA  <- renderImage({ list(src = outA,  contentType = "image/svg+xml") }, deleteFile = FALSE)
        output$imgB  <- renderImage({ list(src = outB,  contentType = "image/svg+xml") }, deleteFile = FALSE)
        output$imgAB <- renderImage({ list(src = outAB, contentType = "image/svg+xml") }, deleteFile = FALSE)

        output$status <- renderText(paste0("¡Listo! Imágenes guardadas en: ", pid_dir))

        incProgress(0.05)
      })
    })
  })
}

# ============================================================
# === UI PRINCIPAL: NAVBAR con dos pestañas (módulos) ===
# ============================================================
ui <- navbarPage(
  "Super-App: TopKAT + Connectivity",
  tabPanel("TopKAT Explorer", topkat_ui("topkat")),
  tabPanel("Conectividad Promedio", connectivity_ui("connect"))
)

# ============================================================
# === SERVER PRINCIPAL: Arranca los módulos ===
# ============================================================
server <- function(input, output, session) {
  topkat_server("topkat")
  connectivity_server("connect")
}

shinyApp(ui, server)
