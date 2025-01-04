
##-------------------------------------##
##            METADATA TAB             ##
##-------------------------------------##

tab_METADATA <- tabItem(
  tabName = "Metadata",
   sidebarLayout(
     sidebarPanel(width = 4,
                  #h4("Create new annotation from selection"),
                #  
                #  textInput(inputId = "add_var_name", label = "New variable name:"),
                #                
                #  actionButton(inputId = "save_selection", "Add selection"),
                #               
                #   br(), br(), 
                                
                  h4("Alter pre-existing annotation"),
                                
                  selectInput(inputId = "select_var_annotation",
                              label = "Select variable",
                              choices = c()),
                                
                  textInput(inputId = "change_var_name", label = "Change variable name to:"),
           
                  selectInput(inputId = "select_group_annotation",
                              label = "Select group",
                              choices = list(),
                              selected = NULL),
                  
                  textInput(inputId = "select_new_annotation", label = "Change group name to:"),
                   
                  actionButton(inputId = "save_new_annotation", "Save changes"),
                                
                  br(),  br(),
                                
                  h4("Add cluster information"),
                                
                  selectInput(inputId = "select_clustering_analysis",
                              label = "Select clustering analysis",
                              choices = ""),
                
                  selectInput(inputId = "select_clustering_annotation",
                              label = "Select resolution",
                              choices = ""),
                         
                  actionButton(inputId = "add_clustering_metadata", label = "Add to metadata")
                  
                  ), #close sidebarPanel 
     mainPanel(  
       ) #close mainPanel    
     ) #close sidebarLayout
  
) # close tabItem               






