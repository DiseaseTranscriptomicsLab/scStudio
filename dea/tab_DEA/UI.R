
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
                 
                 h3("Plotting"),
                 
                 selectInput(inputId = "select_dea", label = "Select DEA", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 h4("Scatter Plot"),
                 
                 textInput(inputId = "dea_cutoff", label = "Expressing cells cutoff (%)", value = 90),
                 
                 h4("Volcano plot"),
                 
                 sliderInput(inputId = "volcano_pvalue", "-log10(Adjusted p-value)", 
                             min = 0, max = 100, step = 0.01 , value = 1.30103, ticks = TRUE),
                 
                 sliderInput(inputId = "volcano_fc", "Average Log2FC", 
                             min = 0, max = 10, step = 0.05 , value = 1, ticks = TRUE),
                 
                 h4("Heatmap"),
                 
                 textInput(inputId = "heatmap_topGenes", label = "Top genes (average Log2FC):",
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
                               inline = TRUE)
                 
                 #actionButton(inputId = "go_plot_dea", "Plot results")
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
                          height = 500), 
             type = 5),
             width = NULL)
    ) #close column
    
    ) #close fluidrow
    ) #close mainPanel    
  ) #close sidebarLayout
) #close DEA
















































