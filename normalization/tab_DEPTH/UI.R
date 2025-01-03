##-------------------------------------##
##              NORM TAB               ##
##-------------------------------------##

tab_DEPTH <- tabItem(
  tabName = "Depth Normalization",
  sidebarLayout(
    sidebarPanel(width = 3,
                 selectInput("norm_method", 
                             label = "Select method ", 
                             choices = c("scran", "SCTransform"),
                             selected = "scran"),
                 
                 textInput(inputId = "norm_mat_name", 
                           label = "Normalization ID:", 
                           value = "normalization1"),
                 
                 actionButton(inputId = "run_normalization", "Run normalization"),     
                 
                 br(),br(),
                 
                 h4("Plotting"), 
                 
                 selectInput(inputId = "select_matrix_depth", 
                             label = "Select count matrix", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 selectInput("var_depth_boxplot",
                             label = "Select condition",
                             choices = list(),
                             selected = NULL),
                 
                 sliderInput(inputId = "select_boxplot_range1",
                             label = "Plot limits:",
                             min = 0, max = 10, value = c(1,6)),
                 
                 actionButton(inputId = "generate_norm_boxplots", label = "Generate plots")
                 
    ), #close sidebarPanel
    mainPanel(width = 9,
              
              box(withSpinner(plotlyOutput(outputId = "boxplot_raw_depth", height = 500), type = 8)),
              box(withSpinner(plotlyOutput(outputId = "boxplot_norm_depth", height = 500), type = 8))
              
              
    ) #close mainPanel
  ) #close sidebarLayout
) #close tab_norm
















