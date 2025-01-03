##-------------------------------------##
##              TSNE TAB               ##
##-------------------------------------##

tab_TSNE<- tabItem(
  tabName = "t-SNE",
   sidebarLayout(
     sidebarPanel(width = 3, 
                  radioButtons( inputId = "genesvspcs",
                                label = "Select genes or PCs:",
                                choices = c("Genes", "PCs"),
                                selected = "PCs",
                                inline = TRUE),
                  
                  uiOutput("select_matrix_tsne"),
                  uiOutput("select_features_tsne"),
                  uiOutput("select_pcs_tsne"),
                  uiOutput("ncomponents_tsne"),
                  
                  sliderInput(inputId = "perplexity", label = "Perplexity:",
                              min = 5, max = 100, 
                              step = 5,
                              value = c(30)),
                                      
                  textInput(inputId = "tsne_id", label = "Identifier:", value = "tsne_1"),
                                      
                  actionButton(inputId = "run_tsne", 
                               "Run TSNE"),
                  br(), br(),
                                    
                  h4("Plot 1"),
                                      
                  selectInput(inputId = "choose_tsne1", 
                              label = "Choose TSNE:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  selectInput(inputId = "tsne_color1", 
                              label = "Colour by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  textInput(inputId = "tsne_title1", 
                            label = "Title:",
                            value = "TSNE"),
                  
                  radioButtons( inputId = "tsne1_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                  
                  sliderInput(inputId = "tsne1_pointSize", 
                              label = "Point size",
                              min = 0.1, 
                              max = 5, 
                              step = 0.1,
                              value = c(1)),
                  
                  actionButton(inputId = "plot_tsne1", 
                               "Plot TSNE"),
                  br(), br(),
                                     
                  h4("Plot 2"),
                                      
                  selectInput(inputId = "choose_tsne2", 
                              label = "Choose TSNE:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),

                  selectInput(inputId = "tsne_color2", 
                              label = "Colour by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  textInput(inputId = "tsne_title2", 
                            label = "Title:",
                            value = "TSNE"), 
                                     
                  radioButtons( inputId = "tsne2_density",
                                label = "Density lines:",
                                choices = c("TRUE", "FALSE"),
                                selected = "FALSE",
                                inline = TRUE),
                  
                   sliderInput(inputId = "tsne2_pointSize", 
                               label = "Point size",
                               min = 0.1, 
                               max = 5, 
                               step = 0.1,
                               value = c(1)),
                                     
                  actionButton(inputId = "plot_tsne2", 
                               "Plot TSNE"),
                         ), #close sidebarPanel
      mainPanel(width = 9,
                fluidRow(
                  box(withSpinner(
                    plotlyOutput(outputId = "tsne_plot1", 
                                 height = 800, 
                                 width = 950), 
                                 type = 5))), #close fluidrow
                                   
                br(),br(),
                                   
                fluidRow(
                  box(withSpinner(
                    plotlyOutput(outputId = "tsne_plot2", 
                                 height = 800, 
                                 width = 950), 
                    type = 5), 
                    width = NULL)) #close fluidrow           
                ) #close mainPanel
     ) #close sidebarlayout
) #close tab_TSNE
















