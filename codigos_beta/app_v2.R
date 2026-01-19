# app.R ------------------------------------------------------------------
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
library(here)
library(base64enc)   # necesario para convertir PNG a dataURI

# Helpers
source("R/helpers_io.R")
source("R/helpers_logic.R")

# === RUTA BASE ===
BASE_DIR <- getwd()
DATA_DIR <- file.path(BASE_DIR, "data")

DEFAULT_DATA_FILE <- file.path(DATA_DIR, "datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv")
DEFAULT_OUT_DIR   <- DATA_DIR

# Comparaciones
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
# UI
# ============================================================

ui <- fluidPage(
  titlePanel("TopKAT Explorer — Shiny"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput("data_file", "Ruta CSV metadata", value = DEFAULT_DATA_FILE),
      textInput("out_dir", "Directorio con RDS", value = DEFAULT_OUT_DIR),
      actionButton("load_all", "📥 Cargar metadata y RDS"),
      hr(),
      uiOutput("available_categories_ui"),
      hr(),
      selectInput("cmp_choice", "Elegir comparación", choices = names(comparaciones_list)),
      numericInput("nA", "N muestras grupo A", value = 3, min = 1),
      numericInput("nB", "N muestras grupo B", value = 3, min = 1),
      numericInput("seedA", "Seed Grupo A", value = 123),
      numericInput("seedB", "Seed Grupo B", value = 456),
      actionButton("run_select", "🎯 Seleccionar muestras y generar subset"),
      hr(),
      actionButton("run_topkat", "⚙️ Ejecutar TopKAT y scale_importance"),
      hr(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Resumen",
                 h4("Tabla pid_info"),
                 DTOutput("pid_info_table") %>% withSpinner(),
                 hr(),
                 h4("Conteos por categoría y FA"),
                 DTOutput("counts_table") %>% withSpinner()
        ),
        tabPanel("Selección",
                 h4("Disponibilidad (Grupo A / Grupo B)"),
                 DTOutput("avail_A"),
                 DTOutput("avail_B"),
                 hr(),
                 h4("PIDs seleccionados"),
                 DTOutput("selected_pids_table")
        ),
        tabPanel("Resultados",
                 h4("TopKAT summary"),
                 verbatimTextOutput("topkat_summary"),
                 hr(),
                 h4("scale_importance (PNG generado)"),
                 uiOutput("scale_png_ui")
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

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
    scale_png_path_log = NULL
  )

  # -----------------------------------------------------------
  # LOAD ALL
  # -----------------------------------------------------------
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
    pickerInput("cats","Filtrar categorías:",
                choices=cats, selected=cats, multiple=TRUE,
                options=list(`actions-box`=TRUE))
  })

  # tablas
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

  # -----------------------------------------------------------
  # RUN SELECT
  # -----------------------------------------------------------

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

    pid_df_final <- bind_rows(
      selA %>% mutate(grupo="A"),
      selB %>% mutate(grupo="B")
    ) %>% mutate(PID_chr = as.character(PID))

    subset_dir <- file.path(input$out_dir,
                            paste0("subset_TOPKAT_",format(Sys.time(),"%Y%m%d_%H%M%S")))
    dir.create(subset_dir, recursive=TRUE)

    selected_pids_chr <- as.character(pid_df_final$PID)

    create_subset_files(
      selected_pids_chr,
      rv$rips_full, rv$K0_full, rv$K1_full, rv$PIDs_full,
      subset_dir
    )

    safe_write_csv(pid_df_final,file.path(subset_dir,"PID_seleccionados_final.csv"))

    rv$pid_df_final <- pid_df_final
    rv$subset_dir <- subset_dir

    output$status <- renderText(paste("✔ Subset creado en:", subset_dir))
  })

  output$avail_A <- renderDT({ req(rv$availA); datatable(rv$availA) })
  output$avail_B <- renderDT({ req(rv$availB); datatable(rv$availB) })
  output$selected_pids_table <- renderDT({
    req(rv$pid_df_final); datatable(rv$pid_df_final)
  })

  # -----------------------------------------------------------
  # RUN TOPKAT — AHORA GENERA 2 PNGs
  # -----------------------------------------------------------

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
      #   GENERAR PNGs — NORMAL Y LOG10
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


  # ---------------------------------------------------
  # MOSTRAR LAS FIGURAS PNG EN LA APP
  # ---------------------------------------------------
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

  # -----------------------------------------------------------
  # Summary TopKAT
  # -----------------------------------------------------------
  output$topkat_summary <- renderPrint({
    req(rv$topkat_res)
    cat("p-value global:", rv$topkat_res$overall.pval, "\n\n")
    print(rv$topkat_res$p.vals)
  })

}

shinyApp(ui, server)
