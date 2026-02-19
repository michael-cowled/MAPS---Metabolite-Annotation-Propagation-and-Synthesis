# app.R --- MAPS Shiny Application (Full Pipeline, Production-Grade)

# ---------------------------------------------------------------------------
# USER DEFAULTS
# ---------------------------------------------------------------------------
# Edit these variables to set default paths and IDs when the app launches.
# Use forward slashes (/) for file paths, even on Windows.
# Leave as "" for a blank default.

DEFAULT_DATASET_ID <- ""
DEFAULT_GNPS_TASK  <- ""
DEFAULT_FOLDER     <- ""
DEFAULT_CACHE      <- ""
DEFAULT_LIPIDS     <- ""
DEFAULT_CID_DB     <- ""

# ---------------------------------------------------------------------------
# Dependency Management (Sequential: CRAN first, then GitHub)
# ---------------------------------------------------------------------------
cran_packages <- c(
  "shiny", "shinyjs", "shinycssloaders", "DT", "openxlsx", "tibble", "bslib",
  "dplyr", "tidyr", "stringr", "readr", "reshape2", "ggplot2", "svglite", 
  "readxl", "data.table", "tidyverse", "rvest", "jsonlite", "xml2", 
  "progress", "DBI", "RSQLite", "httr"
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
  
  # Modern Theme
  theme = bs_theme(version = 5, bootswatch = "zephyr", primary = "#2C3E50"),
  
  tags$head(
    tags$style(HTML("
            .shiny-notification-progress { width: 500px !important; font-size: 16px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
            .progress-bar { height: 25px !important; }
            .shiny-notification-message { font-weight: bold; padding-bottom: 5px;}
            .container-fluid { padding-top: 20px; padding-bottom: 40px; }
        "))
  ),
  
  # Clean, centered title
  div(class = "text-center my-4",
      h1("MAPS", class = "display-4 fw-bold text-primary"),
      h4("Metabolite Annotation Propagation & Synthesis", class = "text-muted")
  ),
  
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("1. Data Configuration", class = "card-title mb-3"),
        # The 'value' argument is now tied to the global variables at the top of the script
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
        checkboxInput("lv1_subclasses", "Use Level 1 subclasses", TRUE),
        checkboxInput("lv2_mzmine", "Use MZmine Level 2 annotations", TRUE)
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
                             column(6, card(card_header("Annotation Counts"), card_body(shinycssloaders::withSpinner(plotOutput("plot1"))))),
                             column(6, card(card_header("Mass Error Histogram"), card_body(shinycssloaders::withSpinner(plotOutput("plot2")))))
                           ),
                           fluidRow(
                             column(12, card(class = "mt-3", card_header("Retention Time vs. M/Z"), card_body(shinycssloaders::withSpinner(plotOutput("plot3")))))
                           )
                  ),
                  tabPanel("Run Log", icon = icon("terminal"),
                           card(class = "mt-3 p-3 bg-light",
                                verbatimTextOutput("log")
                           )
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
      
      update_progress(stage = "Pipeline finished.", value = 1)
      status_msg("✔ SUCCESS: MAPS Pipeline Complete. Check Annotations and Plots tabs.")
      updateTabsetPanel(session, "main-tabs", selected = "Annotations")
      
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
  
  output$plot1 <- renderPlot({ req(plot_objects()); plot_objects()$barchart_plot })
  output$plot2 <- renderPlot({ req(plot_objects()); plot_objects()$histogram_plot })
  output$plot3 <- renderPlot({ req(plot_objects()); plot_objects()$bubblechart_plot })
  
  output$download <- downloadHandler(
    filename = function() paste0(input$dataset_id, "_MAPS_annotations.xlsx"),
    content = function(file) openxlsx::write.xlsx(final_annotations(), file)
  )
}

shinyApp(ui, server)