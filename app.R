# app.R
library(shiny)
library(shinyWidgets)
library(DT)
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(stringr)
library(shinycssloaders)

# helpers
source("R/helpers_io.R")
source("R/helpers_logic.R")

# DEFAULT PATHS (ajusta si hace falta)
DEFAULT_DATA_FILE <- "/home/jupyter-luisraul/APP/datos_para_TopKAT_TOTAL_con_FinalPhenotype.csv"
DEFAULT_OUT_DIR   <- "/home/jupyter-luisraul/APP"  # ahí es donde dijiste que guardaste los RDS

# Comparaciones predefinidas (en el mismo orden de tu script)
comparaciones_list <- list(
  "carcinoma FA vs carcinoma NOFA" = list(A = list(type = "carcinoma", FA = "FA"),
                                          B = list(type = "carcinoma", FA = "NOFA")),
  "dysplasia FA vs dysplasia NOFA" = list(A = list(type = "dysplasia", FA = "FA"),
                                          B = list(type = "dysplasia", FA = "NOFA")),
  "stroma_ad_carcinoma FA vs NOFA" = list(A = list(type = "stroma_ad_carcinoma", FA = "FA"),
                                          B = list(type = "stroma_ad_carcinoma", FA = "NOFA")),
  "stroma_ad_dysplasia FA vs NOFA" = list(A = list(type = "stroma_ad_dysplasia", FA = "FA"),
                                          B = list(type = "stroma_ad_dysplasia", FA = "NOFA")),
  "carcinoma vs dysplasia (ignorar FA)" = list(A = list(type = "carcinoma", FA = NULL),
                                               B = list(type = "dysplasia", FA = NULL)),
  "carcinoma FA vs dysplasia FA" = list(A = list(type = "carcinoma", FA = "FA"),
                                        B = list(type = "dysplasia", FA = "FA")),
  "carcinoma NOFA vs dysplasia NOFA" = list(A = list(type = "carcinoma", FA = "NOFA"),
                                            B = list(type = "dysplasia", FA = "NOFA")),
  "stroma_ad_carcinoma FA vs carcinoma FA" = list(A = list(type = "stroma_ad_carcinoma", FA = "FA"),
                                                  B = list(type = "carcinoma", FA = "FA")),
  "stroma_ad_carcinoma FA vs stroma_ad_dysplasia FA" = list(A = list(type = "stroma_ad_carcinoma", FA = "FA"),
                                                            B = list(type = "stroma_ad_dysplasia", FA = "FA"))
)

ui <- fluidPage(
  titlePanel("TopKAT Explorer — Shiny"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput("data_file", "Ruta CSV metadata (data_file)", value = DEFAULT_DATA_FILE),
      textInput("out_dir", "Directorio con RDS (out_dir)", value = DEFAULT_OUT_DIR),
      actionButton("load_all", "📥 Cargar metadata y RDS (rips, K0, K1, PIDs)"),
      hr(),
      uiOutput("available_categories_ui"),
      hr(),
      selectInput("cmp_choice", "Elegir comparación", choices = names(comparaciones_list)),
      numericInput("nA", "Número muestras grupo A", value = 3, min = 1, step = 1),
      numericInput("nB", "Número muestras grupo B", value = 3, min = 1, step = 1),
      numericInput("seedA", "Seed Grupo A", value = 123, min = 0),
      numericInput("seedB", "Seed Grupo B", value = 456, min = 0),
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
                 h4("Tabla pid_info (una fila por PID)"),
                 DTOutput("pid_info_table") %>% withSpinner(),
                 hr(),
                 h4("Conteos por categoría y FA"),
                 DTOutput("counts_table") %>% withSpinner()
        ),
        tabPanel("Selección",
                 h4("Disponibilidad (Grupo A / Grupo B)"),
                 DTOutput("avail_A") ,
                 DTOutput("avail_B"),
                 hr(),
                 h4("PIDs seleccionados"),
                 DTOutput("selected_pids_table")
        ),
        tabPanel("Resultados",
                 h4("TopKAT summary"),
                 verbatimTextOutput("topkat_summary"),
                 hr(),
                 h4("scale_importance (pval vs thresh)"),
                 plotOutput("scale_plot")
        )
      )
    )
  )
)

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
    selA = NULL,
    selB = NULL,
    pid_df_final = NULL,
    subset_dir = NULL,
    topkat_res = NULL,
    res_scale_import = NULL
  )

  observeEvent(input$load_all, {
    output$status <- renderText("Cargando archivos ...")
    tryCatch({
      data_file <- input$data_file
      out_dir <- input$out_dir

      # 1) metadata
      data1.df <- safe_read_csv(data_file)

      # 2) build pid_info
      pid_info <- build_pid_info(data1.df, archivo_col = "archivo_base")

      # 3) read rips & kernels & PIDs
      rips_path <- file.path(out_dir, "rips_list.rds")
      K0_path   <- file.path(out_dir, "K_dim0.rds")
      K1_path   <- file.path(out_dir, "K_dim1.rds")
      PIDs_path <- file.path(out_dir, "PIDs.rds")

      rips_full <- safe_read_rds(rips_path)
      K0_full <- safe_read_rds(K0_path)
      K1_full <- safe_read_rds(K1_path)
      PIDs_full <- safe_read_rds(PIDs_path)

      # assign
      rv$data1 <- data1.df
      rv$pid_info <- pid_info
      rv$rips_full <- rips_full
      rv$K0_full <- K0_full
      rv$K1_full <- K1_full
      rv$PIDs_full <- PIDs_full
      rv$data_loaded <- TRUE

      output$status <- renderText("✔ Carga completa")
    }, error = function(e){
      output$status <- renderText(paste("Error al cargar:", e$message))
    })
  })

  output$available_categories_ui <- renderUI({
    req(rv$data_loaded)
    cats <- sort(unique(rv$pid_info$category))
    pickerInput("cats", "Filtrar categorías (mostrar en tabla):", choices = cats, selected = cats, multiple = TRUE, options = list(`actions-box` = TRUE))
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
    colnames(tab) <- c("category", "FA_status", "n")
    datatable(tab)
  })

  # when user clicks run_select: compute availability, select samples reproducibly and save subset RDSs
  observeEvent(input$run_select, {
    req(rv$data_loaded)
    output$status <- renderText("Determinando disponibilidad y seleccionando muestras...")
    cmp_sel <- comparaciones_list[[input$cmp_choice]]

    # filter pid_info based on categories selected in sidebar (input$cats)
    pid_info <- rv$pid_info
    if (!is.null(input$cats)) pid_info <- pid_info %>% filter(category %in% input$cats)
    # prepare dispA & dispB
    cmpA <- list(type = cmp_sel$A$type, FA = cmp_sel$A$FA)
    cmpB <- list(type = cmp_sel$B$type, FA = cmp_sel$B$FA)

    dispA <- pid_info %>% filter(category == cmpA$type)
    if (!is.null(cmpA$FA)) dispA <- dispA %>% filter(FA_status == cmpA$FA)
    dispB <- pid_info %>% filter(category == cmpB$type)
    if (!is.null(cmpB$FA)) dispB <- dispB %>% filter(FA_status == cmpB$FA)

    rv$availA <- dispA
    rv$availB <- dispB

    # validate counts
    if (nrow(dispA) < 1 || nrow(dispB) < 1) {
      output$status <- renderText("⚠️ Alguno de los grupos no tiene muestras disponibles.")
      return(NULL)
    }
    if (input$nA > nrow(dispA) || input$nB > nrow(dispB)) {
      output$status <- renderText("⚠️ nA o nB mayor que disponibilidad. Ajústelos.")
      return(NULL)
    }

    # reproducible selection
    selA <- select_reproducible(dispA, input$nA, input$seedA)
    selB <- select_reproducible(dispB, input$nB, input$seedB)

    # combine final pid_df
    pid_df_final <- bind_rows(selA %>% mutate(grupo = "A"),
                              selB %>% mutate(grupo = "B")) %>%
      mutate(PID_chr = as.character(PID)) %>%
      select(grupo, PID, PID_chr, archivo_base, category, FA_status)

    # create subset folder
    subset_dir <- file.path(input$out_dir, paste0("subset_TOPKAT_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

    # create subsets using helpers
    tryCatch({
      selected_pids_chr <- as.character(pid_df_final$PID)
      create_subset_files(selected_pids_chr,
                          rv$rips_full, rv$K0_full, rv$K1_full, rv$PIDs_full,
                          subset_dir)
      # Save pid_df_final
      safe_write_csv(pid_df_final, file.path(subset_dir, "PID_seleccionados_final.csv"))
      rv$selA <- selA; rv$selB <- selB; rv$pid_df_final <- pid_df_final; rv$subset_dir <- subset_dir
      output$status <- renderText(paste("✔ Subset creado en:", subset_dir))
    }, error = function(e) {
      output$status <- renderText(paste("Error creando subset:", e$message))
    })
  })

  output$avail_A <- renderDT({
    req(rv$availA)
    datatable(rv$availA)
  })

  output$avail_B <- renderDT({
    req(rv$availB)
    datatable(rv$availB)
  })

  output$selected_pids_table <- renderDT({
    req(rv$pid_df_final)
    datatable(rv$pid_df_final)
  })

  # RUN TOPKAT and scale_importance on current subset
  observeEvent(input$run_topkat, {
    req(rv$subset_dir)
    output$status <- renderText("Ejecutando TopKAT (~ puede tardar)...")
    tryCatch({
      # load subset
      rips_sub  <- readRDS(file.path(rv$subset_dir, "rips_list_subset.rds"))
      K0_sub    <- readRDS(file.path(rv$subset_dir, "K_dim0_subset.rds"))
      K1_sub    <- readRDS(file.path(rv$subset_dir, "K_dim1_subset.rds"))
      PIDs_sub  <- readRDS(file.path(rv$subset_dir, "PIDs_subset.rds"))

      # prepare K.list for TopKAT
      K.list <- list(K0 = K0_sub, K1 = K1_sub)

      # create y vector from pid_df_final
      pid_df_final <- rv$pid_df_final
      # Ensure rows are ordered as PIDs_sub
      pid_df_final <- pid_df_final %>% arrange(match(as.character(PID), PIDs_sub))
      y <- ifelse(pid_df_final$grupo == "A", 1, 0)

      # run TopKAT
      library(TopKAT)
      res <- TopKAT(y = y, X = NULL, K.list = K.list, omega.list = c(0, 0.5, 1), outcome.type = "binary")
      rv$topkat_res <- res
      # save result
      saveRDS(res, file.path(rv$subset_dir, "TopKAT_result.rds"))

      # run scale_importance
      res_scale_import <- scale_importance(pd.list = rips_sub,
                                          y = y,
                                          omega.list = c(0, 0.5, 1),
                                          threshold = 1000,
                                          PIDs = seq_along(PIDs_sub),
                                          outcome.type = "binary")
      rv$res_scale_import <- res_scale_import
      saveRDS(res_scale_import, file.path(rv$subset_dir, "scale_importance_result.rds"))

      output$status <- renderText(paste0("✔ TopKAT (p=", signif(res$overall.pval, 3), ") y scale_importance ejecutados y guardados en: ", rv$subset_dir))
    }, error = function(e) {
      output$status <- renderText(paste("Error TopKAT/scale_importance:", e$message))
    })
  })

  output$topkat_summary <- renderPrint({
    req(rv$topkat_res)
    print(rv$topkat_res$overall.pval)
    print(rv$topkat_res$p.vals)
  })

  output$scale_plot <- renderPlot({
    req(rv$res_scale_import)
    dfp <- data.frame(thresh = rv$res_scale_import$threshold.seq, pval = rv$res_scale_import$pvals)
    ggplot(dfp, aes(x = thresh, y = pval)) + geom_point() + scale_y_log10() +
      geom_vline(xintercept = rv$res_scale_import$min.thresh, linetype = "dashed") +
      theme_bw() + xlab(expression(epsilon)) + ylab("p-value (log scale)")
  })
}

shinyApp(ui, server)
