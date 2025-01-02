
##-------------------------------------##
##              GENES QC               ##
##-------------------------------------##

tab_QC_GENES <- tabItem(
  tabName = "Quality Control", 
  tabPanel("Features",
         actionButton(inputId = "remove_features", "Remove features"),
         br(),br(),
         sidebarLayout(
           sidebarPanel(width = 2, 
                        h5("Filtering options"),
                        textInput(inputId = "select_nr_counts", label = "Total counts:",
                                    value = 0),
                        textInput(inputId = "select_nr_cells", label = "Number of expressing cells:",
                                    value = 0),
                        h5("Plotting options"),
                        sliderInput(inputId = "top_genes_n", min = 5, max = 50, value = 20, step = 5,
                                    label = "Top most expressed genes:"),
                        
                        #selectInput(inputId = "top_genes_sample", 
                        #            label = "Subset:", 
                        #            choices = "orig.ident"),
                        
                        radioButtons(inputId = "top_genes_scale",
                                     label = "Scale",
                                     choices = c("Linear", "Log2"),
                                     selected = "Linear")
                        
                       
           ), #close sidebarPanel
           mainPanel(width = 10,
                     
                     
                     column(width = 6,
                       box(withSpinner(plotlyOutput(outputId = "discarded_features"), type = 8))
                       ),
                     column(width = 6,
                            box(withSpinner(plotOutput(outputId = "topGenes", width = "500px", 
                                                       height = "1000px"), type = 8))
                     ) #close column
           ) #close mainPanel
         ) #close sidebarLayout
      ) # close tabPanel Feature inspection  
  ) #close tab_QC





