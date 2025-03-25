
##-------------------------------------##
##              DEA TAB                ##
##-------------------------------------##

tab_DEA<- tabItem(
  tabName = "Gene Ranking",
  textOutput(outputId = "session_id"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h3("Differential Expression Analysis"),
                 
                 selectInput(inputId = "select_matrix_dea", label = "Select matrix", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 selectInput(inputId = "dea_condition", label = "Select conditions to test:", 
                             choices = c("-"), multiple = TRUE, width = "100%"),
                 
                 selectInput(inputId = "dea_group1", label = "Group 1:", 
                             choices = c("-"), multiple = TRUE, width = "100%"),
                 
                 radioButtons( inputId = "dea_grouping1",
                               label = "Grouping type:",
                               choices = c("Union", "Intersection"),
                               selected = "Union",
                               inline = TRUE),
                 
                 selectInput(inputId = "dea_group2", label = "Group 2:", 
                             choices = c("-"), multiple = TRUE, width = "100%"),
                 
                 radioButtons( inputId = "dea_grouping2",
                               label = "Grouping type:",
                               choices = c("Union", "Intersection"),
                               selected = "Union",
                               inline = TRUE),
                 
                 selectInput(inputId = "select_method_dea", label = "Select method", 
                             choices = c("MAST", "t", "bimod", "wilcox", "roc"), 
                             multiple = TRUE, width = "100%",
                             selected = "MAST"),
                 
                 textInput(inputId = "dea_id", label = "Identifier:",
                           value = "dea_1"),
                 
                 actionButton(inputId = "run_dea", "Run DEA"),
                 
                  
                 
                 h2("Plotting"),
                 
                 selectInput(inputId = "select_dea", label = "Select DEA", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                  
                 
                 h3("Scatter Plot"),
                 
                 textInput(inputId = "dea_cutoff", label = "Expressing cells cutoff (%)", value = 90),
                 
                  
                 
                 h3("Heatmap"),
                 
                 selectInput(inputId = "heatmap_metric", label = "Order by:", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 textInput(inputId = "heatmap_topGenes", label = "Top genes:",
                           value = 5),
                 
                 radioButtons( inputId = "dea_scale_heatmap",
                               label = "Min-Max scaling:",
                               choices = c("row", "column", "none"),
                               selected = "none",
                               inline = TRUE),
                 
                 radioButtons( inputId = "cluster_heatmap",
                               label = "Cluster:",
                               choices = c("row","column", "both", "none"),
                               selected = "none",
                               inline = TRUE),
                
                 
                 h3("Volcano plot"),
                 
                 selectInput(inputId = "volcano_X", label = "x-axis:", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 selectInput(inputId = "volcano_Y", label = "y-axis:", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 uiOutput("volcano_X_limits"),
                 uiOutput("volcano_Y_limits")
                 
        
  ), #close sidebarPanel
  
  mainPanel(
    fluidRow(
      column(width = 4,
             box(div(withSpinner(
               DT::dataTableOutput(outputId = "table_dea"),
               type = 7),
               style = "font-size: 60%; width: 80%"), width = NULL) 
      ), #close column  
      column(width = 4),
      column(width = 4,
             actionButton("get_dea_plot", "Generate Scatter Plot"),
             box(withSpinner(
               plotlyOutput(outputId = "dea_plot"), 
               type = 5),
               width = NULL)
             ), #close column
      ), #close fluidrow
    
    br(),br(),br(),
    
    fluidRow(column(width = 12,
                    actionButton("get_heatmap", "Generate Heatmap"),
                    box(withSpinner(
                      plotOutput(outputId = "heatmap_dea"), 
                      type = 1), 
                      width = NULL)
                    ) #close column
            
             ), #close fluidrow
    br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),
    fluidRow(
    column(width = 12,
           actionButton("get_volcano", "Generate Volcano Plot"),
           box(withSpinner(
             plotlyOutput(outputId = "volcano_dea", 
                          height = 800), 
             type = 5),
             width = NULL)
    ) #close column
    
    ) #close fluidrow
    ) #close mainPanel    
  ) #close sidebarLayout
) #close DEA
















































