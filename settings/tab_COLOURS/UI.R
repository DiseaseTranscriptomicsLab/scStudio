##-------------------------------------##
##             COLOURS TAB             ##
##-------------------------------------##

tab_COLOURS <- tabItem(
  tabName = "Colours",
  textOutput(outputId = "session_id"),
   sidebarLayout(
     sidebarPanel(width = 4,
                  
                  selectInput(inputId = "select_discrete_col_var",
                              label = "Select variable",
                              choices = NULL),
              
                  selectInput(inputId = "select_discrete_col_group",
                              label = "Select group",
                              choices = list(),
                              selected = NULL),
           
                  colourpicker::colourInput("select_discrete_col","", "#FFA9FA"),
             
                  br(), 
                  
                  actionButton(inputId = "save_discrete_col", "Save"),
                  
                  textOutput(outputId = "show_discrete_cols")
                  ), #close sidebarPanel 
    
     mainPanel(  
      
       fluidRow(
         plotlyOutput(outputId = "colors_plot")  
         )#close fluidrow
       ) #close mainPanel    
     ) #close sidebarLayout
) #close tabItem               






