##-------------------------------------##
##           DOTPLOTS TAB             ##
##-------------------------------------##
tab_DOTPLOTS <- tabItem(
  tabName = "Dot Plots",
   sidebarLayout(
     sidebarPanel(width = 2,
                                      
                  selectInput(inputId = "select_matrix_dotplot", 
                              label = "Select matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  selectizeInput(inputId = "featplot_dotplot_genes", 
                                 label = "Genes",
                                 choices = NULL, 
                                 selected = NULL, 
                                 multiple = TRUE, 
                                 options = NULL),
                                      
                                      
                  #selectizeInput(inputId = "featplot_dotplot_gsets", 
                  #               label = "Gene set scores",
                  #               choices = NULL, 
                  #               selected = NULL, 
                  #               multiple = TRUE, 
                  #               options = NULL),
                                      
                 # radioButtons( inputId = "featplot_dotplot_scale",
               #                 label = "Scale:",
               #                 choices = c("TRUE", "FALSE"),
               #                 selected = "TRUE",
               #                 inline = TRUE),
                                      
                  selectInput(inputId = "select_featplot_dotplot_var", 
                              label = "Group by:",
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                                      
                  #selectInput(inputId = "select_featplot_dotplot_subvar", 
                  #            label = "Sub-group by:", 
                  #            choices = c("-"), 
                  #            multiple = FALSE, 
                  #            width = "100%"),
                                      
                  #radioButtons( inputId = "featplot_dotplot_reverse",
                  #              label = "Reverse plot:",
                  #              choices = c("TRUE", "FALSE"),
                  #              selected = "FALSE",
                  #              inline = TRUE),
                                      
                  sliderInput(inputId = "featplot_dotplot_width", 
                              label = "Plot width:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  sliderInput(inputId = "featplot_dotplot_height", 
                              label = "Plot heigth:",
                              min = 500, 
                              max = 3000, 
                              step = 100, 
                              value = 500),
                                      
                  actionButton(inputId = "run_featplot_dotplot", "Plot"),
                                      
                  br(),br()
                  ), #close sidebarPanel
     mainPanel(
       fluidRow(
         column(width = 6,
                box(withSpinner(
                  plotOutput(outputId = "dotplot"), 
                  color = "#C91004", 
                  type = 5), 
                  width = "100%")
                ) #close column
         ) #close fluidrow
       ) #close mainPanel    
     ) #close sidebarLayout
) #close DOTPLOTS













