#UI based calculation of gibbs free energy calculations

library(shiny)

#import calculation file
source("free_energy_calcs.R")

ui <- fluidPage(

    titlePanel("Gibbs Energy Estimations"),

    sidebarLayout(
      
        sidebarPanel(
          h4(style = "margin-bottom: 20px","IMPORTANT:All chemistry must be in CHNOSZ format"),
          
          fileInput("physicochemistry_file","Upload physicochemisty csv",accept=".csv"),
          
          selectInput("ionic_strength","Select ionic strength",
                      choices = c(0.001,0.01,0.1,0.7),selected = "" ),
          
          textAreaInput("reactions","Enter reactions", 
                      placeholder ="Format: name = species1,species2,...|coefs1,coefs2,...
                      \nExample:\nFeOx_NitrateRed = Fe+2,NO3-,H+,Fe+3,NO2-,H2O|-2,-1,-2,2,1,1",
                      rows = 7),
          
          textInput("minerals","List of minerals (comma separate, no spaces!)",value =""),
          
          textInput("output_name","Name of output file:",value = ""),
          
          actionButton("calc_button","Calculate"),
          
          downloadButton("download_results","Download results")
          ),

        # show a summary of the table (top 10 rows)
        mainPanel(
          
           h4("Results (preview)"),
           
           tableOutput("preview_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    results <- reactiveVal(NULL)
    
    observeEvent(input$calc_button, {
      
      # there needs to be an input file
      req(input$physicochemistry_file)
      
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
        reactions[[name]] <- list(species,coefs)
      }
      
      reaction_names <- names(reactions_list)
      
      # call calculations
      RealGibbs <- free_energy_calcs(redox,reactions_list,reaction_names,minerals,I)
      
      return(results(RealGibbs))
      }
    )
    
    output$preview_table <- renderTable({
      head(results())
    })
    
    output$download_results <- downloadHandler(
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
