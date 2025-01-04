
##-------------------------------------##
##             SUBSET TAB              ##
##-------------------------------------##

tab_SUBSET<- tabItem(
  tabName = "Metadata",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectInput(inputId = "subset_var", 
                                      label = "Select conditions:", 
                                      choices = c("-"), 
                                      multiple = TRUE, 
                                      width = "100%"),
                                
                          selectInput(inputId = "subset_group", 
                                      label = "Select groups:", 
                                            choices = c("-"), 
                                      multiple = TRUE, 
                                      width = "100%"),
                          
                          selectInput(inputId = "choose_plot_subset", 
                                      label = "Select plot type:", 
                                      choices = c("-"), 
                                      multiple = FALSE, 
                                      width = "100%"),
                          
                          radioButtons(inputId = "subset_grouping",
                                       label = "Grouping type:",
                                       choices = c("Union", "Intersection"),
                                       selected = "Union",
                                       inline = TRUE),
                                
                          actionButton(inputId = "subset_dataset", 
                                       label = "Subset dataset")
                                
                   ), #close sidebarPanel
                   mainPanel(width = 9,
                             withSpinner(plotlyOutput(outputId = "subset_tsne", 
                                                      height = 800, 
                                                      width = 950),
                                         type = 5)
                   ) #close mainPanel
                 ) #close sidebarLayout


  
  ) # close tabItem               






