##-------------------------------------##
#####            REPORT             #####
##-------------------------------------##

tab_SESSIONINFO <- tabItem(
  tabName = "Session Information",
  
           #sidebarLayout(
             #sidebarPanel(width = 2,
                          
                          
             #), #close sidebarPanel
             
             #mainPanel(             
            column(width = 6,
           verbatimTextOutput(outputId = "session_info")
            ) #close column
                 
                 
                 
               #) #close fluidrow
               
             #) #close mainPanel    
           #) #close sidebarLayout
) #close tabItem




