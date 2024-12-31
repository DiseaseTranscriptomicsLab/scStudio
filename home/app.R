#INFORMATION-----------------------------


#LOAD LIBRARIES ------------------------
library(shiny)
library(shinyjs)
print("Sucessfully loaded libraries.")


#LOAD TABS-------------------------------
source("UI.R")
print("Successfully loaded tabs.")

#USER INFERFACE ----------------------------

ui <- fluidPage(
  ## theme ----------------------------------------
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  
## navbarPage --------------------------------------- 
  navbarPage("",
             tabPanel("SCS",
                      tab_HOME    
             )#close HOME
  )#close navbarPage 
  
) #close UI 

##-------------------------------------##  


server <- function(input, output, session) {
  
  
}#close server

shinyApp(ui = ui, server = server)








