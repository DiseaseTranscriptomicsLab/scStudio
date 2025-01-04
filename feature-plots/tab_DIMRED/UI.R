##-------------------------------------##
##             DIMRED TAB              ##
##-------------------------------------##
tab_DIMRED <- tabItem(
  tabName = "Dimensionality reduction plots",
  textOutput(outputId = "session_id"),
   sidebarLayout(
     sidebarPanel(width = 3,
                  
                  selectInput(inputId = "select_featplot_type_dimred", 
                              label = "Select type of plot", 
                              choices = c("pca", "tsne", "umap"), 
                              multiple = FALSE, 
                              width = "100%",
                              selected = "tsne"),
                                      
                  selectInput(inputId = "select_featplot_dimred", 
                              label = "Select reduction", 
                              choices = c(""), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectInput(inputId = "select_featplot_xx", 
                              label = "Select xx", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectInput(inputId = "select_featplot_yy", 
                              label = "Select yy", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectInput(inputId = "select_matrix_dimred", 
                              label = "Select matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectizeInput(inputId = "featplot_genes", 
                                 label = "Genes",
                                 choices = NULL, 
                                 selected = NULL, 
                                 multiple = TRUE, 
                                 options = NULL),
                  
                  selectizeInput(inputId = "featplot_gsets", 
                                 label = "Gene set scores",
                                 choices = NULL, 
                                 selected = NULL, 
                                 multiple = TRUE, options = NULL),
                  
                  radioButtons( inputId = "featplot_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                                     
                  sliderInput(inputId = "featplot_dimred_width", 
                              label = "Plot width:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  sliderInput(inputId = "featplot_dimred_height", 
                              label = "Plot heigth:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  actionButton(inputId = "run_featplot", 
                               "Plot"),
                                      
                  br(),br()
                  ), #close sidebarPanel
     mainPanel(             
       fluidRow(
         column(width = 6,
                box(
                  withSpinner(
                    plotOutput(outputId = "grid_featplots"),
                    color = "#C91004", 
                    type = 5), 
                  width = "100%")
                ) #close column
         ) #close fluidrow
       ) #close mainPanel    
     ) #close sidebarLayout
) #close FEATPLOTS













