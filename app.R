#UI based calculation of gibbs free energy calculations

library(shiny)

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
             HTML('<strong>IMPORTANT:</strong> All species must be in <a href="https://chnosz.net/vignettes/anintro.html" target=>CHNOSZ</a> format')),
          
          selectInput("ionic_strength","Select ionic strength",
                      choices = c(0.001,0.01,0.1,0.7),selected = "" ),
          
          
          textAreaInput("reactions","Enter reactions", 
                      placeholder ="Format: name = species1,species2,...|coefs1,coefs2,...
                      \nExample:\nFeOx_NitrateRed = Fe+2,NO3-,H+,Fe+3,NO2-,H2O|-2,-1,-2,2,1,1",
                      rows = 7),
          
          textInput("minerals","List of minerals (comma separate, no spaces!)",value =""),
          
          textInput("output_name","Name of output file:",value = ""),
          
          # run calculation
          actionButton("calc_button","Calculate"),
          ),
          
        # show a summary of the table (top 10 rows)
        mainPanel(
          
          fluidRow(
            
            column(width = 6, h4("Tabular results (preview)"),
                   tags$hr(style = "border-top: 1px solid #333; margin-top: 0px; margin-bottom: 25px;"),
                   tableOutput("preview_table"),
                   downloadButton("download_table","Download tabular results"))
          )
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
    results <- reactiveVal(NULL)
    
    observeEvent(input$calc_button, {
      
      # there needs to be an input file and reactions
      req(input$physicochemistry_file)
      req(input$reactions)
      
      #load the input file
      redox <- read.csv(input$physicochemistry_file$datapath, check.names = FALSE)
      
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
      
      return(results(RealGibbs))
      }
    )
    
    # output preview of results
    output$preview_table <- renderTable({
      head(results())
    })
    
    #download results if button pressed
    output$download_table <- downloadHandler(
      filename = function() {
        cleaned <- gsub("\\s+", "_", input$output_filename)
        paste0(cleaned, ".csv")
      },
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      }
    )

}

shinyApp(ui, server)
