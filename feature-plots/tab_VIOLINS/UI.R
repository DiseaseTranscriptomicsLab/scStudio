##-------------------------------------##
##            VIOLINS TAB              ##
##-------------------------------------##
tab_VIOLINS <- tabItem(
  tabName = "Violin plots",
   sidebarLayout(
     sidebarPanel(width = 3,

                  selectInput(inputId = "select_matrix_violins", 
                              label = "Select matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectizeInput(inputId = "featplot_violin_genes", 
                                 label = "Genes",
                                 choices = NULL, 
                                 selected = NULL, 
                                 multiple = TRUE, 
                                 options = NULL),
                                      
                  selectizeInput(inputId = "featplot_violins_gsets", 
                                 label = "Gene set scores",
                                 choices = NULL, 
                                 selected = NULL,
                                 multiple = TRUE, 
                                 options = NULL),
                                      
                  radioButtons(inputId = "featplot_violin_scale",
                               label = "Scale:",
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE",
                               inline = TRUE),
                                      
                  selectInput(inputId = "select_featplot_violin_var", 
                              label = "Group by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectInput(inputId = "select_featplot_violin_subvar", 
                              label = "Sub-group by:", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  sliderInput(inputId = "featplot_violin_ptsize", 
                              label = "Point size:",
                              min = 0, 
                              max = 2, 
                              step = 0.1, 
                              value = 0.1),
                                      
                  sliderInput(inputId = "featplot_violin_width", 
                              label = "Plot width:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  sliderInput(inputId = "featplot_violin_height", 
                              label = "Plot heigth:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  actionButton(inputId = "run_featplot_violin", "Plot"),
                                      
                  br(),br()
                  ), #close sidebarPanel
     mainPanel(
       fluidRow(
         column(width = 6,
                box(withSpinner(
                  plotOutput(outputId = "grid_violinplots"),
                  color = "#C91004", 
                  type = 5),
                  width = "100%")
                ) #close column
         ) #close fluidrow
       ) #close mainPanel    
     ) #close sidebarLayout
  ) #close Violins














