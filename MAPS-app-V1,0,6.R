# app.R --- MAPS Shiny Application v1.06

# ---------------------------------------------------------------------------
# USER DEFAULTS
# ---------------------------------------------------------------------------
DEFAULT_DATASET_ID <- ""
DEFAULT_GNPS_TASK  <- ""
DEFAULT_FOLDER     <- ""
DEFAULT_CACHE      <- "~/MAPS/cid_cache.csv"
DEFAULT_LIPIDS     <- "Y:/MA_BPA_Microbiome/Databases/PubChem/lipids_expanded.tsv"
DEFAULT_CID_DB     <- "Y:/MA_BPA_Microbiome/Databases/PubChem/PubChem_Indexed.sqlite"

# ---------------------------------------------------------------------------
# Dependency Management (Sequential: CRAN first, then GitHub)
# ---------------------------------------------------------------------------
cran_packages <- c(
  "shiny", "shinyjs", "shinycssloaders", "DT", "openxlsx", "tibble", "bslib",
  "dplyr", "tidyr", "stringr", "readr", "reshape2", "ggplot2", "svglite", 
  "readxl", "data.table", "tidyverse", "rvest", "jsonlite", "xml2", 
  "progress", "DBI", "RSQLite", "httr", "plotly"
)

# 1. Install and load CRAN packages sequentially
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing CRAN package:", pkg))
    utils::install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# 2. Install and load GitHub packages (MAPS.Package)
if (!requireNamespace("remotes", quietly = TRUE)) {
  utils::install.packages("remotes")
}

if (!requireNamespace("MAPS.Package", quietly = TRUE)) {
  message("Installing MAPS.Package from GitHub...")
  remotes::install_github("michael-cowled/MAPS-Package-Public") 
}
suppressPackageStartupMessages(library(MAPS.Package))

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
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  useShinyjs(),
  theme = bs_theme(version = 5, bootswatch = "zephyr", primary = "#2C3E50"),
  
  tags$head(
    tags$style(HTML("
            .shiny-notification-progress { width: 500px !important; font-size: 16px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
            .progress-bar { height: 25px !important; }
            .shiny-notification-message { font-weight: bold; padding-bottom: 5px;}
            .container-fluid { padding-top: 20px; padding-bottom: 40px; }
        "))
  ),
  
  div(class = "text-center my-4",
      h1("MAPS", class = "display-4 fw-bold text-primary"),
      h4("Metabolite Annotation Propagation & Synthesis", class = "text-muted")
  ),
  
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("1. Data Configuration", class = "card-title mb-3"),
        textInput("dataset_id", "Dataset ID (Output Name)", value = DEFAULT_DATASET_ID, placeholder = "e.g., Experiment_01"),
        textInput("gnps_task", "GNPS Task ID", value = DEFAULT_GNPS_TASK, placeholder = "32-character ID"),
        textInput("folder", "Processed Data Folder Path", value = DEFAULT_FOLDER, placeholder = "Use forward slashes '/'"),
        textInput("cache", "CID Cache CSV Path (Optional)", value = DEFAULT_CACHE, placeholder = "Use forward slashes '/'"),
        textInput("lipids", "Lipids CSV Path (Optional)", value = DEFAULT_LIPIDS, placeholder = "Use forward slashes '/'"),
        textInput("cid_db", "PubChem SQLite Path (Optional)", value = DEFAULT_CID_DB, placeholder = "Use forward slashes '/'")
      ),
      wellPanel(
        h4("2. Acceptance Thresholds", class = "card-title mb-3"),
        helpText("Minimum scores required to accept an annotation."),
        fluidRow(
          column(6, numericInput("gnps_prob", "GNPS Prob.", 0.7, min=0, max=1, step=0.01)),
          column(6, numericInput("canopus_prob", "CANOPUS Prob.", 0.7, min=0, max=1, step=0.01))
        ),
        fluidRow(
          column(6, numericInput("csi_prob", "CSI:FingerID", 0.64, min=0, max=1, step=0.01)),
          column(6, numericInput("ms2_prob", "MS2Query", 0.7, min=0, max=1, step=0.01))
        ),
        fluidRow(
          column(6, numericInput("rt_tol", "RT Tol. (min)", 0.1, min=0, step=0.01)),
          column(6, numericInput("ppm_tol", "PPM Tol.", 5, min=1, max=20))
        )
      ),
      wellPanel(
        h4("3. Options", class = "card-title mb-3"),
        checkboxInput("plots", "Generate & save plots", TRUE),
        checkboxInput("standardisation", "PubChem standardisation", TRUE),
        checkboxInput("enable_local_db", "Local PubChem DB Lookup", TRUE),
        checkboxInput("enable_api", "API PubChem Enrichment", TRUE),
        checkboxInput("lv1_subclasses", "Use Level 1 subclasses", TRUE),
        checkboxInput("lv2_mzmine", "Use MZmine Level 2 annotations", TRUE),
        checkboxInput("msnovelist", "Use MSNovelist annotations", FALSE)
      ),
      hr(),
      actionButton("run", "▶ Run MAPS Pipeline", class = "btn-primary btn-lg w-100 my-2"),
      downloadButton("download", "⬇ Download Results", class = "btn-outline-primary w-100")
    ),
    
    mainPanel(
      tabsetPanel(type = "pills", id = "main-tabs",
                  tabPanel("Status", icon = icon("info-circle"),
                           card(class = "mt-3 p-3 bg-light", verbatimTextOutput("status"))
                  ),
                  tabPanel("Annotations", icon = icon("table"),
                           class = "mt-3",
                           shinycssloaders::withSpinner(DTOutput("results_table"), color = "#2C3E50")
                  ),
                  tabPanel("Plots", icon = icon("chart-bar"),
                           class = "mt-3",
                           fluidRow(
                             column(6, card(card_header("Feature Counts per Sample"), card_body(shinycssloaders::withSpinner(plotOutput("plot1"))))),
                             column(6, card(card_header("Cumulative Annotations"), card_body(shinycssloaders::withSpinner(plotOutput("plot2")))))
                           ),
                           fluidRow(
                             column(6, card(class = "mt-3", card_header("Class Starburst"), card_body(shinycssloaders::withSpinner(plotOutput("plot4"))))),
                             column(6, card(class = "mt-3", card_header("Metabolite Class Bubble Chart"), card_body(shinycssloaders::withSpinner(plotOutput("plot3")))))
                           )
                  ),
                  tabPanel("Run Log", icon = icon("terminal"),
                           card(class = "mt-3 p-3 bg-light", verbatimTextOutput("log"))
                  )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  status_msg <- reactiveVal("Idle – Ready for configuration.")
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
    updateTabsetPanel(session, "main-tabs", selected = "Run Log")
    
    progress <- shiny::Progress$new()
    progress$set(message = "MAPS Pipeline Initialization", value = 0)
    on.exit(progress$close())
    
    update_progress <- function(stage, value, detail = NULL) {
      progress$set(value = value, message = paste("MAPS Pipeline:", stage), detail = detail)
      status_msg(paste("Running:", stage))
      append_log(stage)
    }
    
    # -----------------------------------------------------------------------
    # LOGGING BLOCK START
    # -----------------------------------------------------------------------
    params_log <- c(
      "========================================",
      "MAPS Pipeline Run Parameters",
      "========================================",
      paste("Timestamp:", Sys.time()),
      paste("Dataset ID:", input$dataset_id),
      paste("GNPS Task ID:", input$gnps_task),
      paste("Processed Data Folder:", input$folder),
      paste("CID Cache Path:", input$cache),
      paste("Lipids File Path:", input$lipids),
      paste("PubChem SQLite DB:", input$cid_db),
      "----------------------------------------",
      "Acceptance Thresholds",
      "----------------------------------------",
      paste("GNPS Probability:", input$gnps_prob),
      paste("CANOPUS Probability:", input$canopus_prob),
      paste("CSI:FingerID Score:", input$csi_prob),
      paste("MS2Query Probability:", input$ms2_prob),
      paste("Retention Time Tolerance (min):", input$rt_tol),
      paste("Mass Tolerance (PPM):", input$ppm_tol),
      "----------------------------------------",
      "Pipeline Options",
      "----------------------------------------",
      paste("Generate Plots:", input$plots),
      paste("PubChem Standardisation:", input$standardisation),
      paste("Use Level 1 Subclasses:", input$lv1_subclasses),
      paste("Use MZmine Level 2:", input$lv2_mzmine),
      paste("Use MSNovelist:", input$msnovelist),
      "========================================"
    )
    
    file_prefix <- ifelse(input$dataset_id == "", "MAPS", input$dataset_id)
    log_file_path <- file.path(input$folder, paste0(file_prefix, "_run_parameters.log"))
    
    tryCatch({
      if (dir.exists(input$folder)) {
        writeLines(params_log, log_file_path)
        append_log(paste("Saved parameter log to:", log_file_path))
      } else {
        append_log("Warning: Folder does not exist yet. Parameter log not saved.")
      }
    }, error = function(e) {
      append_log(paste("Warning: Could not save parameter log:", e$message))
    })
    # -----------------------------------------------------------------------
    # LOGGING BLOCK END
    # -----------------------------------------------------------------------
    
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
        msnovelist        = input$msnovelist,
        modification_db   = modification_db,
        updateProgress    = update_progress,
        enable_local_db   = input$enable_local_db,
        enable_api        = input$enable_api
      )
      
      final_annotations(fa)
      update_progress(stage = "MAPS annotation table generated", value = 0.90)
      
      if (input$plots) {
        update_progress(stage = "Generating plots", value = 0.95)
        results_list <- MAPS.Package::make_plots(fa, input$folder)
        plot_objects(results_list)
      }
      
      update_progress(stage = "Pipeline finished.", value = 1)
      status_msg("✔ SUCCESS: MAPS Pipeline Complete. Check Annotations and Plots tabs.")
      # Note: Leaving selected = Annotations so the user sees the table immediately, 
      # but they can click Plots manually.
      updateTabsetPanel(session, "main-tabs", selected = "Plots") 
      
    }, error = function(e) {
      status_msg("✖ MAPS failed – see log for details.")
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
    datatable(dt, filter = 'top', style = 'bootstrap5', options = list(scrollX = TRUE, pageLength = 15))
  })
  
  # Render outputs securely using print() and req() to avoid errors on NULLs
  output$plot1 <- renderPlot({ 
    req(plot_objects(), plot_objects()$barchart_plot)
    print(plot_objects()$barchart_plot) 
  })
  
  output$plot2 <- renderPlot({ 
    req(plot_objects(), plot_objects()$histogram_plot)
    print(plot_objects()$histogram_plot) 
  })
  
  output$plot3 <- renderPlot({ 
    req(plot_objects(), plot_objects()$bubblechart_plot)
    print(plot_objects()$bubblechart_plot) 
  })
  
  output$plot4 <- renderPlot({ 
    req(plot_objects(), plot_objects()$starburst_plot)
    print(plot_objects()$starburst_plot) 
  })
  
  output$download <- downloadHandler(
    filename = function() paste0(input$dataset_id, "_MAPS_annotations.xlsx"),
    content = function(file) openxlsx::write.xlsx(final_annotations(), file)
  )
}

shinyApp(ui, server)