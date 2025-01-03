##-------------------------------------##
##            DIMRED TAB               ##
##-------------------------------------##

tab_FEATURE_SELECTION <- tabItem(
  tabName = "Feature Selection",
  textOutput(outputId = "session_id"),
    sidebarLayout(
      sidebarPanel(width = 3,
                   selectInput(inputId = "select_matrix_fselect", 
                               label = "Select count matrix", 
                               choices = c("-"), 
                               multiple = FALSE, 
                               width = "100%"),  
                   textInput(inputId = "hvgs_id", 
                             label = "Analysis ID:",
                             value = "top_hvgs"),
                   
                   actionButton(inputId = "run_fs", "Run"),
                   
                          br(),br(),
                   
                 
                   
                          textInput(inputId = "top_hvg", 
                                    label = "Top highly variable genes:",
                                    value = 1000),
                   
                          textInput(inputId = "max_pval", 
                                    label = "P-value:",
                                    value = 0.05),
                   
                 
                          selectizeInput(inputId = "gene_highlight", 
                                  label = "Highlight a gene:",
                                  choices = NULL, 
                                  selected = NULL, 
                                  multiple = FALSE, 
                                  options = NULL),
                          
                          actionButton(inputId = "subset_features", 
                                       "Add subset of features")
                         
                         ), #close sidebarPanel
      
                         mainPanel(width = 9,
                              fluidRow(
                               column(width = 6,
                                box(withSpinner(plotlyOutput(outputId = "varvsmean", height = 500), type = 8),
                                  width = NULL)), #close column
                               column(width = 6,
                                box(withSpinner(plotlyOutput(outputId = "pvalvsbio", height = 500), type = 8),
                                  width = NULL)), #close column
                           ) #close fluidRow
                           #fluidRow(
                          #   column(width = 7,
                          #          box(div(withSpinner(DT::dataTableOutput(
                          #            outputId = "hvgs_table"),
                          #            type = 7),
                          #            style = "font-size: 60%; width: 80%"), 
                          #            width = NULL)
                          #   ) #close column  
               
                           #) #close fluidrow
                         ) #close mainPanel
                       ) #close sidebarLayout
  ) #close tab_DIMRED
















