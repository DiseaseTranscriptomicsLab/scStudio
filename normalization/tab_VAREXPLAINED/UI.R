##-------------------------------------##
##              NORM TAB               ##
##-------------------------------------##

tab_VAREXPLAINED <- tabItem(
  tabName = "Variance Explained",
  textOutput(outputId = "session_id"),
  br(),br(),
  tabsetPanel(type = "tabs",
     
                       sidebarLayout(
                         sidebarPanel(width = 3,
                          selectInput(inputId = "select_matrix_varexp", label = "Select count matrix", 
                                                  choices = c("-"), multiple = FALSE, width = "100%"),  
                   
                          selectInput(inputId = "variance_vars", label = "Select variables to inspect", 
                                                  choices = c("-"), multiple = TRUE, width = "100%"),  
                          
                          actionButton(inputId = "run_variance", "Run"),
                         ), #close sidebarPanel
                         mainPanel(width = 9,
                              box(withSpinner(plotlyOutput(outputId = "var_exp", height = 500), type = 8),
                                  width = NULL), #close box
                         ) #close mainPanel
                       ) #close sidebarLayout
  ) # close tabsetPanel
) #close tab_norm
















