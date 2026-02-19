# app.R --- MAPS Shiny Application (Full Pipeline, Production-Grade)

library(shiny)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(openxlsx)
library(tibble)

# ---------------------------------------------------------------------------
# EXTERNAL FILE MAPPING
# This allows Shiny to "see" your Archive folder as a web directory named /assets/
# ---------------------------------------------------------------------------
shiny::addResourcePath(prefix = "assets", directoryPath = "~/")

# ---------------------------------------------------------------------------
# Modification table
# ---------------------------------------------------------------------------
modification_db <- tibble::tribble(
  ~Modification, ~Mass.Change,
  "Isomeric", 0,
  "Dihydro", 2.016,
  "Dehydro", -2.016,
  "Hydroxylated", 15.9949,
  "Deoxygenated", -15.9949,
  "Dihydroxylated", 31.9898,
  "Dideoxygenated", -31.9898,
  "Hydration", 18.0106,
  "Dehydration", -18.0106,
  "Carboxy", 43.9898,
  "Decarboxy", -43.9898,
  "Methylated", 14.0157,
  "Demethylated", -14.0157,
  "Formylated", 27.9949,
  "Deformylated", -27.9949,
  "Acetylated", 42.0106,
  "Deacetylated", -42.0106,
  "Sulfated", 79.9568,
  "Desulfated", -79.9568,
  "Phosphorylated", 79.9663,
  "Dephosphorylated", -79.9663,
  "Glycosylated", 162.0528,
  "Deglycosylated", -162.0528,
  "Pentose addition", 132.0423,
  "Pentose loss", -132.0423,
  "Deoxyhexose addition", 146.0579,
  "Deoxyhexose loss", -146.0579,
  "Glucuronidation", 176.0321,
  "GlcNAc", 203.0794
)

# ---------------------------------------------------------------------------
# Package bootstrap
# ---------------------------------------------------------------------------
check_and_install <- function(packages, github_packages = list(),
                              github_pat = Sys.getenv("GITHUB_PAT")) {
  if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% names(github_packages)) {
        remotes::install_github(repo = github_packages[[pkg]], auth_token = github_pat)
      } else {
        utils::install.packages(pkg, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

#Round 1
github_packages <- ""
required_packages <- c(
  "dplyr", "tidyr", "stringr", "readr",
  "reshape2", "ggplot2", "svglite", "readxl", "data.table",
  "openxlsx", "tidyverse", "rvest", "jsonlite", "xml2", "progress",
  "DBI", "RSQLite", "httr", "DT", "shiny", "shinyjs", "shinycssloaders"
)
check_and_install(required_packages, github_packages)

#Round 2
github_packages <- list("MAPS.Package" = "michael-cowled/MAPS-Package-Public")
required_packages <- "MAPS.Package"
check_and_install(required_packages, github_packages)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),

  tags$head(
    tags$style(HTML("
            .shiny-notification-progress { width: 500px !important; font-size: 16px; }
            .progress-bar { height: 25px !important; }
            .shiny-notification-message { font-weight: bold; }
            /* Align the logo and title text */
            .title-container { display: flex; align-items: center; margin-bottom: 20px; }
            .title-logo { height: 60px; margin-right: 20px; }
        "))
  ),

  # Custom Title Panel with Logo
  tags$div(class = "title-container",
           tags$img(src = "assets/MAPS-logo.png", class = "title-logo"),
           tags$h2("MAPS – Metabolite Annotation Propagation & Synthesis")
  ),

  sidebarLayout(
    sidebarPanel(
      h4("Dataset configuration"),
      textInput("dataset_id", "Dataset ID", ""),
      textInput("gnps_task", "GNPS task ID", ""),
      textInput("folder", "Processed data folder", ""),
      textInput("cache", "CID Cache CSV path", ""),
      textInput("lipids", "Lipids CSV path", ""),
      textInput("cid_db", "PubChem SQLite path", ""),

      hr(),
      h4("Acceptance thresholds"),
      numericInput("gnps_prob", "GNPS probability", 0.7, min = 0, max = 1, step = 0.01),
      numericInput("canopus_prob", "CANOPUS probability", 0.7, min = 0, max = 1, step = 0.01),
      numericInput("csi_prob", "CSI:FingerID probability", 0.64, min = 0, max = 1, step = 0.01),
      numericInput("ms2_prob", "MS2Query probability", 0.7, min = 0, max = 1, step = 0.01),
      numericInput("rt_tol", "RT tolerance (min)", 0.1, min = 0, step = 0.01),
      numericInput("ppm_tol", "PPM tolerance", 5, min = 1, max = 20),

      hr(),
      h4("Optional configurations"),
      checkboxInput("plots", "Generate & save plots", TRUE),
      checkboxInput("standardisation", "PubChem standardisation", TRUE),
      checkboxInput("lv1_subclasses", "Use Level 1 subclasses", TRUE),
      checkboxInput("lv2_mzmine", "Use MZmine Level 2 annotations", TRUE),

      hr(),
      actionButton("run", "▶ Run MAPS", class = "btn-primary btn-lg"),
      br(), br(),
      downloadButton("download", "⬇ Download annotations"),
      width = 3
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Status", verbatimTextOutput("status")),
        tabPanel("Annotations", shinycssloaders::withSpinner(DTOutput("results_table"))),
        tabPanel("Plots",
                 shinycssloaders::withSpinner(plotOutput("plot1")),
                 shinycssloaders::withSpinner(plotOutput("plot2")),
                 shinycssloaders::withSpinner(plotOutput("plot3"))
        ),
        tabPanel("Run log", verbatimTextOutput("log"))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  status_msg <- reactiveVal("Idle")
  run_log <- reactiveVal(character())
  final_annotations <- reactiveVal(NULL)
  plot_objects <- reactiveVal(NULL)

  append_log <- function(msg) {
    run_log(c(run_log(), paste(Sys.time(), "-", msg)))
  }

  observeEvent(input$run, {
    disable("run")
    status_msg("Preparing to run MAPS...")
    run_log(character())

    progress <- shiny::Progress$new()
    progress$set(message = "MAPS Pipeline Initialization", value = 0)
    on.exit(progress$close())

    update_progress <- function(stage, value, detail = NULL) {
      progress$set(value = value, message = paste("MAPS Pipeline:", stage), detail = detail)
      status_msg(paste("Running:", stage))
      append_log(stage)
    }

    tryCatch({
      update_progress(stage = "Starting execution and data loading", value = 0.05)

      fa <- MAPS.Package::MAPS(
        dataset.id        = input$dataset_id,
        folder            = input$folder,
        cache.location    = input$cache,
        lipids.location   = input$lipids,
        gnps.task.id      = input$gnps_task,
        gnps.prob         = input$gnps_prob,
        canopus.prob      = input$canopus_prob,
        csi.prob          = input$csi_prob,
        ms2query.prob     = input$ms2_prob,
        ppm.tol           = input$ppm_tol,
        rt.tol            = input$rt_tol,
        cid_database_path = ifelse(input$cid_db == "", NULL, input$cid_db),
        standardisation   = input$standardisation,
        lv1.subclasses    = input$lv1_subclasses,
        lv2.mzmine        = input$lv2_mzmine,
        modification_db   = modification_db,
        updateProgress    = update_progress
      )

      final_annotations(fa)
      update_progress(stage = "MAPS annotation table generated", value = 0.90)

      if (input$plots) {
        update_progress(stage = "Generating plots", value = 0.95)
        results_list <- MAPS.Package::make_plots(fa, input$folder)
        plot_objects(results_list)
      }

      update_progress(stage = "Pipeline finished. Displaying results.", value = 1)
      status_msg("✔ SUCCESS: MAPS Pipeline Complete")

    }, error = function(e) {
      status_msg("✖ MAPS failed – see log")
      append_log(paste("ERROR:", e$message))
      progress$close()
    })

    enable("run")
  })

  output$status <- renderText(status_msg())
  output$log <- renderText(paste(run_log(), collapse = "\n"))

  output$results_table <- DT::renderDT({
    req(final_annotations())
    dt <- final_annotations()
    if ("Samples" %in% colnames(dt)) dt <- dt[, setdiff(colnames(dt), "Samples")]
    datatable(dt, filter = 'top', options = list(scrollX = TRUE, pageLength = 10))
  })

  output$plot1 <- renderPlot({ req(plot_objects()); plot_objects()$barchart_plot })
  output$plot2 <- renderPlot({ req(plot_objects()); plot_objects()$histogram_plot })
  output$plot3 <- renderPlot({ req(plot_objects()); plot_objects()$bubblechart_plot })

  output$download <- downloadHandler(
    filename = function() paste0(input$dataset_id, "_MAPS_annotations.xlsx"),
    content = function(file) openxlsx::write.xlsx(final_annotations(), file)
  )
}

shinyApp(ui, server)
