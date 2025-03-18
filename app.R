
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
      menuItem("scATACseq Analyzer", tabName = "input", icon = icon("edit")),
      conditionalPanel(
        condition = "input.tab=='input'",
        div(
          fileInput("file", "Upload File", multiple = FALSE, accept = c('.rds')),
          actionButton("reset", "Reset", icon = icon("undo"),
                       style = "color: #fff; background-color: #dc3545; width: 87.25%"),
          actionButton("run", "Run", icon = icon("play"),
                       style = "color: #fff; background-color: #28a745; width: 87.25%"),
          selectizeInput("region", "Open Chromatin Region", choices = NULL, multiple = FALSE)
          
          )
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              tabPanel("Instruction", uiOutput("instruction_content"))
      ),
      tabItem(tabName = "input",
              tabsetPanel(id = "main_tabs",
                          tabPanel("Instruction", includeMarkdown("markdown/instructions.md"))
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
  
  obj <- reactiveVal(NULL)  # Store Seurat object to Reactive
  
  observe({
    shinyjs::disable("region")
  })
  
  observe({
    if (!is.null(input$file)) {
      shinyjs::enable("run")
      print(input$file)
    } else {
      shinyjs::disable("run")
    }
  })
  
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("run")
    shinyjs::disable("region")
    removeTab("main_tabs", "DA_Table")
    removeTab("main_tabs", "DA_Plot")
    removeTab("main_tabs", "Coverage_Plot")
      

  })
  
  observeEvent(input$file, {
    seurat_obj <- load_seurat_obj(input$file$datapath)
    
    if (!is.vector(seurat_obj)) {
      obj(seurat_obj)
      shinyjs::enable("run")
    } else {
      shinyjs::disable("run")
      
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See details below:</h5><br>",
             paste(unlist(seurat_obj), collapse = "<br></br>"))
      ))
    }
  })
  
  
  # Store computed results
  observeEvent(input$run, {
    shinyjs::disable("run")
    show_modal_spinner(text = paste("Finding Markers ..."))
    
    if (is.vector(obj())) {
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded.</h5><br>",
             paste(unlist(obj()), collapse = "<br>"))
      ))
      shinyjs::enable("run")
    } else {
      da_peaks <- create_da_table(obj())
      shinyjs::enable("region")
      
      output$da_table <- renderDataTable({
        DT::datatable(da_peaks)
      })
      
      output$accessibility_region_plot <- renderPlot({
        req(input$region)
        selected_region <- input$region[1]
        create_accessibility_plot(obj(), region = selected_region)
      })
      
      output$coverage_plot <- renderPlot({
        req(input$region)
        selected_region <- input$region[1]
        create_coverage_plot(obj(), region = selected_region)
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
          plot <- create_accessibility_plot(obj(), input$region)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
      output$downloadCoveragePlot <-downloadHandler(
        filename = function(){
          paste0(input$region,"_CoveragePlot.png")
        },
        content =  function(file){
          plot <- create_coverage_plot(obj(), input$region)
          ggsave(filename = file, width=10, height = 5, type = "cairo")
        }
      )
      
            
      # show the da peaks table
      insertTab(
        inputId = "main_tabs",
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
        inputId = "main_tabs",
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
        inputId = "main_tabs",
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
}

shinyApp(ui, server)