##-------------------------------------##
#####            REPORT             #####
##-------------------------------------##

tab_REPORT <- tabItem(
  tabName = "Report",
  textOutput(outputId = "session_id"),
  tabsetPanel(type = "tabs",
              tabPanel("Report",
                       sidebarLayout(
                         sidebarPanel(width = 4,
                                      selectInput(inputId = "select_analysis_type", label = "Type of analysis:", 
                                                  choices = c("-"), multiple = FALSE, width = "100%"),
                                      selectInput(inputId = "select_analysis_id", label = "ID:", 
                                                  choices = c("-"), multiple = FALSE, width = "100%"),
                                      
                         ), #close sidebarPanel
                         
                         mainPanel(             
                           fluidRow(
                             verbatimTextOutput(outputId = "report")
                             
                           ) #close fluidrow
                           
                         ) #close mainPanel    
                       ) #close sidebarLayout
              ) #close report tab
  ) # close tabpanel
) #close tabItem




