
##-------------------------------------##
##                QC TAB               ##
##-------------------------------------##

tab_QC_DOUBLETS <- tabItem(
  tabName = "Doublets",
  
  actionButton(inputId = "calculate_doublets", "Calculate doublets"),
  actionButton(inputId = "remove_doublets", "Remove doublets"),
  br(), br(), 
  sidebarLayout(
    sidebarPanel(width = 2, 
       selectInput("doublet_method", 
                   label = "Select method ", 
                   choices = list("scDblFinder"),
                   selected = "scDblFinder"),
       
               selectInput("doublet_channel",
                           label = "Sample IDs:",
                           choices = list(),
                           selected = NULL),
       p("If you have multiple cell captures, 
         it is preferable to look for doublets separately for each sample"),
       
       sliderInput(inputId = "selectDoubletCutoff", label = "Doublet score:",
                   min = 0, max = 1,
                   value = 0.9, step = 0.01),
       p("The higher the score, the more likely that the cell is a doublet.")
  ), #close sidebarPanel
  
  mainPanel(width = 10,
            fluidRow(
                   box(withSpinner(plotlyOutput(outputId = "doublet_plots"), type = 8),
                       width = NULL)), #close fluidrow
            fluidRow(div(style = "height:200px;")), #close fluidrow
            fluidRow(
                   box(withSpinner(dataTableOutput(outputId = "doublet_metrics", width = "50%"), type = 8), width = "NULL") 
            ), #close fluidrow
         ) #close mainPanel
       ) #close sidebarLayout
) #close tab_Doublets





