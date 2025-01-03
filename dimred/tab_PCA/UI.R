##-------------------------------------##
##              PCA TAB                ##
##-------------------------------------##

tab_PCA <- tabItem(
  tabName = "PCA",
   sidebarLayout(
     sidebarPanel(width = 3,
                  selectInput(inputId = "select_matrix_PCA", 
                              label = "Select count matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "select_features", 
                              label = "Select features", 
                              choices = c("All"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  radioButtons( inputId = "do_scaling",
                                label = "Scale to unit variance:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                  
                  sliderInput(inputId = "ncomponents", 
                              label = "Number of components:",
                              min = 2, 
                              max = 500, 
                              step = 1,
                              value = c(100)),
                                      
                  textInput(inputId = "pca_id", 
                            label = "Identifier:",
                            value = "pca_1"),
                  
                  actionButton(inputId = "run_pca", "Run PCA"),
                                      
                  br(),br(),
                  
                  h4("Plotting"),
                  
                  selectInput(inputId = "choose_pca", 
                              label = "Choose PCA:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "pca_x", 
                              label = "X-axis:", 
                              choices = c(1), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "pca_y", 
                              label = "Y-axis:", 
                              choices = c(2), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "pca_color", 
                              label = "Colour by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  radioButtons( inputId = "pca_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                  
                  textInput(inputId = "pca_title", 
                            label = "Title:",
                            value = "PCA"),
                                  
                  sliderInput(inputId = "pca_pointSize", 
                              label = "Point size",
                              min = 0.1, max = 5, step = 0.1,
                              value = c(1)),
                  
                  actionButton(inputId = "plot_pca", "Plot PCA"),
        ), #close sidebarPanel
     
        mainPanel(width = 9,
                  fluidRow(
                    box(withSpinner(
                      plotlyOutput(outputId = "pca_plot", 
                                   height = 800, 
                                   width = 950), 
                                   type = 5))), #close fluidrow
                  
                                br(),br(),
                  
                                fluidRow(
                                      box(withSpinner(
                                            plotlyOutput(outputId = "scree_plot"), 
                                            type = 5), 
                                          width = NULL)) #close fluidrow 

                    ) #close mainPanel
                ) #close sidebarLayout
) #close tab_PCA
















