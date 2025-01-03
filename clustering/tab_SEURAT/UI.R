
##-------------------------------------##
##              SEURAT TAB             ##
##-------------------------------------##

tab_SEURAT <- tabItem(
  tabName = "Seurat",
  textOutput(outputId = "session_id"),
   sidebarLayout(
     sidebarPanel(width = 2,
                  h4("Clustering"),
                  p("Run a new clustering analysis."),
                
                  radioButtons( inputId = "genesvspcs_cgraph",
                               label = "Use genes or PCs:",
                               choices = c("Genes", "PCs"),
                               selected = "PCs",
                               inline = TRUE),
                  
                 uiOutput("select_matrix_seurat"),
                 
                 uiOutput("select_features_cgraph"),
                 
                 uiOutput("select_pcs_cgraph"),
                 
                 uiOutput("ncomponents_cgraph"),
                 
                 textInput(inputId = "cgraph_id", label = "Identifier:",
                           value = "clust_1"),
                 
                 actionButton(inputId = "run_cgraph", "Run clustering"),
                 
                 h4("Marker genes"),
                 
                 p("Get the marker genes for the currently selected analysis."),
                 
                 selectInput(inputId = "select_matrix_run_mkgenes", 
                             label = "Select matrix", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 selectInput(inputId = "select_method_dea_cgraph", 
                             label = "Select method", 
                             choices = c("MAST", "t", "bimod", "wilcox", "roc"), 
                             multiple = TRUE, 
                             width = "100%",
                             selected = "MAST"),
                 
                 textInput(inputId = "minLogFC_cgraph", 
                           label = "Min Log2FC:",
                           value = 0.50),
                 
                 actionButton(inputId = "run_mkgenes", "Get marker genes"),
                 
                 br(),br(),
                 
                 h4("Plotting"),
                 
                 selectInput(inputId = "select_clustering_analysis", 
                             label = "Select clustering analysis", 
                             choices = c(""), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 selectInput(inputId = "select_dimred_cgraph", 
                             label = "Select dimred", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 sliderInput(inputId = "resolution_cgraph", "Resolution", min = 0.1, max = 2, step = 0.1 , value = 0.1, ticks = FALSE),
                 p("Increased values of resolution lead to a greater number of clusters."),
    
                 br(),
                 
                 h4("Table options:"),
                 
                 selectInput(inputId = "select_dea_method_cgraph", 
                             label = "Method", 
                             choices = "", 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 textInput(inputId = "selected_cluster_cgraph", label = "Cluster", value = ""),
                 
                 textInput(inputId = "adjP_cgraph", label = "Bonferroni-adjusted p-value", value = "0.01")
                 ), #close sidebarPanel
     mainPanel(             
       fluidRow(
         column(width = 6,
                box(withSpinner(
                  plotlyOutput(
                    outputId = "cgraph_plot"), 
                  type = 5),
                  width = NULL),
                ), #close column
         column(width = 6,
                box(withSpinner(plotOutput(
                  outputId = "cgraph_clustree", 
                  height = 600), 
                  type = 5),
                  width = NULL),
                ), #close column
         ), #close fluidrow
       
       br(),br(),br(),
       
       fluidRow(
         column(width = 7,
                
                h4("Marker genes"), #save list option
                p("Select a gene to view the violin plot."),
                
                box(div(withSpinner(DT::dataTableOutput(
                  outputId = "mkgs_cgraph"),
                  type = 7),
                  style = "font-size: 60%; width: 80%"), 
                  width = NULL)
                ), #close column  
         column(width = 5,
                box(withSpinner(plotlyOutput(
                  outputId = "violins_cgraph"), 
                  type = 1), 
                  width = NULL)
                ) #close column
         ) #close fluidrow
       ) #close mainPanel    
     ) #close sidebarLayout
) #close CLUSTERS
















































