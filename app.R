
source('global.R')

ui <- dashboardPage(
  dashboardHeader(title = "scATACseq Analysis"),
  dashboardSidebar(
    tags$head(
      tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
    ),
    sidebarMenu(
      id = "tab",
      useShinyjs(),
      menuItem("Home Page", tabName = "home", icon = icon("list")),
      # scATACseq Analyzer
      menuItem("scATACseq Analyzer", tabName = "scATAC_input", icon = icon("edit")),
      conditionalPanel(
        condition = "input.tab == 'scATAC_input'",
        div(
          fileInput("atac_file", "Upload ATAC File", multiple = FALSE, accept = c('.rds')),
          actionButton("atac_reset", "Reset", icon = icon("undo"),
                       style = "color: #fff; background-color: #dc3545; width: 87.25%"),
          actionButton("atac_run", "Run", icon = icon("play"),
                       style = "color: #fff; background-color: #28a745; width: 87.25%"),
          selectizeInput("region", "Open Chromatin Region", choices = NULL, multiple = FALSE)
        )
      ),
      # scRNAseq Analyzer
      menuItem("scRNAseq Analyzer", tabName = "scRNA_input", icon = icon("edit")),
      conditionalPanel(
        condition = "input.tab == 'scRNA_input'",
        div(
          fileInput("rna_file", "Upload RNA File", multiple = FALSE, accept = c('.rds')),
          actionButton("rna_reset", "Reset", icon = icon("undo"),
                       style = "color: #fff; background-color: #dc3545; width: 87.25%"),
          actionButton("rna_run", "Run", icon = icon("play"),
                       style = "color: #fff; background-color: #28a745; width: 87.25%")
        )
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              tabPanel("Instruction", uiOutput("instruction_content"))
      ),
      tabItem(tabName = "scATAC_input",
              tabsetPanel(id = "atac_tabs",
                          tabPanel("Instruction", includeMarkdown("markdown/scATACseq_instructions.md"))
              )
      ),
      tabItem(tabName = "scRNA_input",
              tabsetPanel(id = "rna_tabs",
                          tabPanel("Instruction", includeMarkdown("markdown/scRNAseq_instructions.md"))
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$instruction_content <- renderUI({
    md_text <- paste(readLines("markdown/landing.md"), collapse = "\n")
    HTML(markdown::markdownToHTML(text = md_text, fragment.only = TRUE))
  })
  
  options(shiny.maxRequestSize = 1000 * 1024^2)
  
  # Store Seurat object to Reactive
  atac_obj <- reactiveVal(NULL)
  rna_obj <- reactiveVal(NULL)
  
  observe({
    shinyjs::disable("region")
  })
  
  observe({
    if (!is.null(input$atac_file)) {
      shinyjs::enable("atac_run")
    } else {
      shinyjs::disable("atac_run")
    }
  })
  
  observe({
    if (!is.null(input$rna_file)) {
      shinyjs::enable("rna_run")
    } else {
      shinyjs::disable("rna_run")
    }
  })
  
  observeEvent(input$atac_reset, {
    shinyjs::reset("atac_file")
    shinyjs::disable("atac_run")
    shinyjs::disable("region")
    removeTab("atac_tabs", "DA_Table")
    removeTab("atac_tabs", "DA_Plot")
    removeTab("atac_tabs", "Coverage_Plot")
      

  })
  observeEvent(input$rna_reset, {
    shinyjs::reset("rna_file")
    shinyjs::disable("rna_run")
    removeTab("rna_tabs", "UMAP")
    removeTab("rna_tabs", "Gene Expression")
  })
  # Load scATAC
  observeEvent(input$atac_file, {
    seurat_atac <- load_seurat_obj(input$atac_file$datapath)
    
    if (!is.vector(seurat_atac)) {
      atac_obj(seurat_atac)
      shinyjs::enable("atac_run")
    } else {
      shinyjs::disable("atac_run")
      
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See details below:</h5><br>",
             paste(unlist(seurat_atac), collapse = "<br></br>"))
      ))
    }
  })
  # Load scRNA
  observeEvent(input$rna_file, {
    seurat_rna <- load_seurat_obj(input$rna_file$datapath)
    
    if (!is.vector(seurat_rna)) {
      rna_obj(seurat_rna)
      shinyjs::enable("rna_run")
    } else {
      shinyjs::disable("rna_run")
      
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See details below:</h5><br>",
             paste(unlist(seurat_rna), collapse = "<br></br>"))
      ))
    }
  })
  
  # Compute atac results
  observeEvent(input$atac_run, {
    shinyjs::disable("atac_run")
    show_modal_spinner(text = paste("Finding Markers ..."))
    
    if (is.vector(atac_obj())) {
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded.</h5><br>",
             paste(unlist(atac_obj()), collapse = "<br>"))
      ))
      shinyjs::enable("atac_run")
    } else {
      da_peaks <- create_da_table(atac_obj())
      shinyjs::enable("region")
      
      output$da_table <- renderDataTable({
        DT::datatable(da_peaks)
      })
      
      output$accessibility_region_plot <- renderPlot({
        req(input$region)
        selected_region <- input$region[1]
        create_accessibility_plot(atac_obj(), region = selected_region)
      })
      
      output$coverage_plot <- renderPlot({
        req(input$region)
        selected_region <- input$region[1]
        create_coverage_plot(atac_obj(), region = selected_region)
      })

      updateSelectizeInput(session, "region",
                           choices = rownames(da_peaks),
                           selected = if (length(rownames(da_peaks)) > 0) rownames(da_peaks)[1] else NULL,
                           server = TRUE)
      
      output$downloadDAPlot <-downloadHandler(
        filename = function(){
          paste0(input$region,"_DAPlot.png")
        },
        content =  function(file){
          plot <- create_accessibility_plot(atac_obj(), input$region)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
      output$downloadCoveragePlot <-downloadHandler(
        filename = function(){
          paste0(input$region,"_CoveragePlot.png")
        },
        content =  function(file){
          plot <- create_coverage_plot(atac_obj(), input$region)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
      
            
      # show the da peaks table
      insertTab(
        inputId = "atac_tabs",
        tabPanel(
          "DA_Table",
          fluidRow(
            column(
              width = 12,
              bslib::card(
                full_screen = TRUE,
                card_header("Differential Accessibility Table"),
                card_body(
                  DT::dataTableOutput(outputId = "da_table")
                )
              )
            )
          )
        )
      )
      insertTab(
        inputId = "atac_tabs",
        tabPanel(
          "DA_Plot",
          fluidRow(
            column(
              width = 12,
              bslib::card(
                full_screen = TRUE,
                card_header("Cell Expression Level"),
                card_body(
                  plotOutput(outputId = "accessibility_region_plot"),
                  downloadButton("downloadDAPlot","Download Differeitial Accessibility Plot")
                )
              )
            )
          )
        )
      )
      insertTab(
        inputId = "atac_tabs",
        tabPanel(
          "Coverage_Plot",
          fluidRow(
            column(
              width = 12,
              bslib::card(
                full_screen = TRUE,
                card_header("Coverage Plot"),
                card_body(
                  plotOutput(outputId = "coverage_plot"),
                  downloadButton("downloadCoveragePlot","Download Covergae Plot")
                )
              )
            )
          )
        )
      )
      remove_modal_spinner()
    }
  })
  
  # Compute rna results
  observeEvent(input$rna_run, {
    shinyjs::disable("rna_run")
    show_modal_spinner(text = paste("Finding Markers ..."))
    
    if (is.vector(rna_obj())) {
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded.</h5><br>",
             paste(unlist(rna_obj()), collapse = "<br>"))
      ))
      shinyjs::enable("rna_run")
    } else {

      output$umap <- renderPlot({
        if (!is.null(input$metadata_col)) {
          create_metadata_UMAP(rna_obj(), input$metadata_col)
        }
      })
      
      output$featurePlot <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot(rna_obj(), input$gene)
        }
      })
      
      output$download_umap <-downloadHandler(
        filename = function(){
          paste0(input$metadata_col,"_UMAP.png")
        },
        content =  function(file){
          plot <- create_metadata_UMAP(rna_obj(), input$metadata_col)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
      output$downloadFeaturePlot <-downloadHandler(
        filename = function(){
          paste0(input$gene,"_FeaturePlot.png")
        },
        content =  function(file){
          plot <- create_feature_plot(rna_obj(), input$gene)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
      
      # show the da peaks table
      insertTab(
        inputId = "rna_tabs",
        tabPanel(
          "UMAP",
          fluidRow(
            column(
              width = 8,
              bslib::card(
                full_screen = TRUE,
                card_header("Cells Clustering - UMAP"),
                card_body(
                  plotOutput(outputId = "umap"),
                  downloadButton("download_umap","Download UMAP")
                )
              )
            ),
            column(
              width = 4,
              selectizeInput("metadata_col","Metadata Column",colnames(rna_obj()@meta.data))
            )
          )
        )
      )
      insertTab(
        inputId = "rna_tabs",
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8,
              bslib::card(
                full_screen = TRUE,
                card_header("Cells Clustering - UMAP"),
                card_body(
                  plotOutput(outputId = "featurePlot"),
                  downloadButton("downloadFeaturePlot","Download Feature Plot")
                )
              )
            ),
            column(
              width = 4,
              selectizeInput("gene","Genes", rownames(rna_obj()))
            )
          )
        )
      )
      remove_modal_spinner()
    }
  })
}

shinyApp(ui, server)