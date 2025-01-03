##-------------------------------------##
##             UMAP TAB                ##
##-------------------------------------##

tab_UMAP <- tabItem(
  tabName = "UMAP",
   sidebarLayout(
     sidebarPanel(width = 3, 
                  radioButtons( inputId = "genesvspcs_umap",
                                label = "Use genes or PCs:",
                                choices = c("Genes", "PCs"),
                                selected = "PCs",
                                inline = TRUE),
                                                 
                  uiOutput("select_matrix_umap"),
                                                 
                  uiOutput("select_features_umap"),
                                                 
                  uiOutput("select_umap_pcs"),
                                                 
                  uiOutput("ncomponents_umap"),
                  
                  sliderInput(inputId = "n_neighbors_umap", 
                              label = "# Neighbors",
                              min = 3, 
                              max = 1000, 
                              step = 1, #need to update max value
                              value = c(15)),
                  
                  sliderInput(inputId = "min_dist_umap", 
                              label = "Minimum distance",
                              min = 0, max = 1, step = 0.01,
                              value = c(0.01)),
                  
                  textInput(inputId = "umap_id", 
                            label = "Identifier:", 
                            value = "umap_1"),
                  
                  actionButton(inputId = "run_umap", "Run UMAP"),
                                                 
                  br(), br(),
                                                 
                  h4("Plot 1"),
                                                 
                  selectInput(inputId = "choose_umap1", 
                              label = "Choose UMAP:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "umap_color1", 
                              label = "Colour by:", 
                              choices = c("-"), 
                              multiple = FALSE,
                              width = "100%"),
                                                 
                  textInput(inputId = "umap_title1", 
                            label = "Title:",
                            value = "UMAP"),
                                                
                  radioButtons( inputId = "umap1_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                                                 
                  sliderInput(inputId = "umap1_pointSize", 
                              label = "Point size",
                              min = 0.1, 
                              max = 5, 
                              step = 0.1,
                              value = c(1)),
                                                 
                  actionButton(inputId = "plot_umap1", 
                               "Plot UMAP"),
                                                 
                  br(), br(),
                                                 
                  h4("Plot 2"),
                                                 
                  selectInput(inputId = "choose_umap2", 
                              label = "Choose UMAP:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "umap_color2", 
                              label = "Colour by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  textInput(inputId = "umap_title2", 
                            label = "Title:",
                            value = "UMAP"), 
                                                
                  radioButtons( inputId = "umap2_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                                                 
                  sliderInput(inputId = "umap2_pointSize", 
                              label = "Point size",
                              min = 0.1, 
                              max = 5, 
                              step = 0.1,
                              value = c(1)),
                                                 
                  actionButton(inputId = "plot_umap2", "Plot UMAP"),
                  ), #close sidebarPanel
     mainPanel(width = 9,
               fluidRow(
                 box(withSpinner(
                   plotlyOutput(outputId = "umap_plot1", 
                                height = 800, 
                                width = 950), 
                   type = 5))), #close fluidrow
                                  
                br(),br(),
                                   
               fluidRow(
                 box(withSpinner(
                   plotlyOutput(outputId = "umap_plot2", 
                                height = 800, 
                                width = 950), 
                   type = 5), width = NULL)) #close fluidrow           
               ) #close mainPanel
     ) #close sidebarLayout
) #close tab_UMAP
















