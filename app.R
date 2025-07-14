#UI based calculation of gibbs free energy calculations

library(shiny)
library(ggplot2)
library(plotly)
library(tidyr)
library(dplyr)

#import calculation file
source("free_energy_calcs.R")

ui <- fluidPage(

    titlePanel("Gibbs Energy Estimations"),

    sidebarLayout(
      
        sidebarPanel(
          
          #download physicochemical data
          h4(style = "margin-bottom: 0px","Download file with physicochemical data"),
          tags$hr(style = "border-top: 1px solid #333; margin-top: 0px; margin-bottom: 25px;"),

          HTML('<strong>Upload csv (see '),
          downloadLink("download_template", label = "template"),HTML('):</strong>'),
          
          fileInput("physicochemistry_file", NULL, accept = ".csv"),
          
          #inputs for calculations
          h4(style = "margin-bottom: 0px","Setup for calculations"),
          tags$hr(style = "border-top: 1px solid #333; margin-top: 0px; margin-bottom: 25px;"),
          
          h5(style = "margin-bottom: 20px",
             HTML('<strong>IMPORTANT:</strong> All species must be in <a href="https://chnosz.net/vignettes/anintro.html" target=>CHNOSZ</a> nomenclature')),
          
          selectInput("ionic_strength","Select ionic strength",
                      choices = c(0.001,0.01,0.1,0.7),selected = "" ),
          
          
          textAreaInput("reactions","Enter reactions", 
                      placeholder ="Format: name = species1,species2,...|coefs1,coefs2,...
                      \nExample:\nFeOx_NitrateRed = Fe+2,NO3-,H+,Fe+3,NO2-,H2O|-2,-1,-2,2,1,1",
                      rows = 7),
          
          textInput("minerals","List of minerals (comma separate, no spaces!)",value =""),
          
          # run calculation
          actionButton("calc_button","Calculate"),
          ),
          
        # show a summary of the table (top 10 rows)
        mainPanel(
          uiOutput("results_ui")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #download template when link is clicked
    output$download_template <- downloadHandler(
      filename = function() {
        "input_file_template.csv"
      },
      content = function(file) {
        file.copy("input_file_template.csv", file)
      }
    )
    
    #perform calculations
    input_data <- reactiveVal(NULL)
    results <- reactiveVal(NULL)
    
    observeEvent(input$calc_button, {
      
      # there needs to be an input file and reactions
      req(input$physicochemistry_file)
      req(input$reactions)
      
      #load the input file
      redox <- read.csv(input$physicochemistry_file$datapath, check.names = FALSE)
      input_data(redox)
      
      #constant setup
      minerals <- unlist(strsplit(input$minerals, ","))
      
      I <- as.numeric(input$ionic_strength)
      
      reaction_data <- strsplit( trimws(input$reactions),"\n")[[1]]
      reactions_list <- list()
      
      for (num in reaction_data){
        parse <- strsplit(num,"=")[[1]]
        name <- trimws(parse[1])
        chem <- strsplit(trimws(parse[2]),"\\|")[[1]]
        species <- strsplit(chem[1], ",")[[1]]
        coefs <- as.numeric(strsplit(chem[2],",")[[1]])
        reactions_list[[name]] <- list(species,coefs)
      }
      
      reaction_names <- names(reactions_list)
      
      # call calculations
      RealGibbs <- free_energy_calcs(redox,reactions_list,reaction_names,minerals,I)
      
      # because R is SO insistent on renaming my columns ...
      colnames(RealGibbs)[1:4] <- c("Sample", "Temperature (C)", "Pressure (bar)", "pH")
      
      return(results(RealGibbs))
      }
    )
    
    # output table preview of results
    output$preview_table <- renderTable({
      req(results())
      head(results())
    })
    
    #download results if button pressed (and computation completed)
    output$download_table <- downloadHandler(
      filename = function() {
        req(results())
        cleaned <- gsub("\\s+", "_", input$table_output_name)
        paste0(cleaned, ".csv")
      },
      content = function(file) {
        req(results())
        write.csv(results(), file, row.names = FALSE)
      }
    )
    
    #generate reactive plot (so it can be saved without redoing the code)
    plot_result <- reactive({
      req(results(), input$y_axis_input)
      
      ycol <- input$y_axis_input
      
      consistents <- c("Sample", "Temperature (C)", "Pressure (bar)", "pH")
      
      #pull out the reactions
      rxns <- setdiff(names(results()), c(consistents))
  
      df <- merge(results(),input_data(),by = consistents)
      
      melted <- df %>%
        pivot_longer(cols = all_of(rxns), names_to = "Reaction", values_to = "DeltaG")
      
      plot <- ggplot(melted, aes(x= DeltaG, y = .data[[ycol]], color = Reaction)) + 
        (if (is.numeric(melted[[ycol]])) geom_path(size = 1) else geom_point(size = 2)) +
        theme_minimal() +
        labs(title = paste("ΔG vs.", ycol), x = "ΔG (kJ/mol)", y = ycol) +
        expand_limits(x = -10)+
        geom_vline(xintercept = 0,linetype = "dashed", color = "grey")
      
      return(plot)
    })
    
    #output the plot
    output$gibbs_plot <- renderPlotly({
      req(plot_result())
      ggplotly(plot_result())
    })
    
    #save the plot
    #output$download_plot <- downloadHandler(
    #  filename = function() {
    #    cleaned <- gsub("\\s+", "_", input$plot_output_name)
    #    paste0(cleaned, ".png")
    #  },
    #  content = function(file) {
    #    ggsave(file, plot = plot_result(), width = 8, height = 6, dpi = 300)
    #  }
    #)
    
    # generate active ui for table and plot results and downloads
    output$results_ui <- renderUI({
      req(results())
      
      tabsetPanel(
        
        tabPanel("Table", h3(tags$strong("Tabular results (preview)"),
               tags$hr(style = "border-top: 1px solid #333; margin-top: 0px; margin-bottom: 25px;"),
               div(style = "height: 300px; overflow-y: auto; border: 1px solid #ccc; padding: 5px;",
                   tableOutput("preview_table")),
               textInput("table_output_name","Name of output file:",value = "tabular_results"),
               downloadButton("download_table","Download tabular results"))),
        tabPanel("Plot", h3(tags$strong("Plot variable vs. \u0394G")),
               tags$hr(style = "border-top: 1px solid #333; margin-top: 0px; margin-bottom: 25px;"),
               #select y-axis for plot
               selectInput("y_axis_input","Select y-axis (options from input data):",choices = names(input_data())),
               h5("You can download as png by clicking the camera icon on the plotly plot."),
               #textInput("plot_output_name","Name of output file:",value = "plotted_results"),
               #downloadButton("download_plot","Download plotted results"),
               plotlyOutput("gibbs_plot",height = "600px"))
      )
    })

}

shinyApp(ui, server)
