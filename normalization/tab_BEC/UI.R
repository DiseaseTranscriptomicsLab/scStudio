##-------------------------------------##
##               BEC TAB               ##
##-------------------------------------##

tab_BEC <- tabItem(
  tabName = "Batch Effect Correction",
  sidebarLayout(
    sidebarPanel(width = 3, 
                 selectInput(inputId = "select_matrix_bec", label = "Select count matrix", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 selectInput("be_method", 
                             label = "Select method", 
                             choices = list("limma"),
                             selected = "limma"),
                 selectInput("batch_label1",
                             label = "Select variable indicating batches",
                             choices = list(),
                             selected = NULL),
                 selectInput("batch_label2",
                             label = "Select variable indicating second series of batches (optional)",
                             choices = list(),
                             selected = NULL),
                 
                 selectInput("adjust_covariate",
                             label = "Select a covariate to be adjusted for (optional)",
                             choices = list(),
                             selected = NULL),
                 
                 selectInput("keep_biological",
                             label = "Select condition to be preserved (optional)",
                             choices = list(),
                             selected = NULL),
                 textInput(inputId = "mat_name", label = "New matrix name:", value = "bec_mat1"),
                 
                 actionButton(inputId = "correct_be", "Correct batch effect"),
                 
                 h4("Plotting"),
                 
                 selectInput(inputId = "select_matrix_plot_bec", label = "Select count matrix", 
                             choices = c("-"), multiple = FALSE, width = "100%"),
                 
                 selectInput("var_be_boxplot",
                             label = "Select condition",
                             choices = list(),
                             selected = NULL),
                 
                 sliderInput(inputId = "select_boxplot_range3", label = "Plot limits:",
                             min = 0, max = 10, value = c(1,6)),
                 
                 actionButton(inputId = "generate_be_boxplots", label = "Generate plots")  
                 
                 
                 
    ), #close sidebarPanel
    mainPanel(width = 9,
              fluidRow(
                box(withSpinner(plotlyOutput(outputId = "boxplot_be", height = 500), type = 8))
              ) #close fluidrow
    ) #close mainPanel
  ) #close sidebarLayout
) #close tab_bec
















