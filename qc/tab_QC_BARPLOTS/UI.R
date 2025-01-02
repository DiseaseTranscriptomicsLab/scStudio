
##-------------------------------------##
##                QC TAB               ##
##-------------------------------------##



tab_QC_BARPLOTS <- tabItem(
  tabName = "Quality Control",
  tabPanel("BARPLOTS",
          sidebarLayout(
          sidebarPanel(width = 2,
            selectInput(inputId = "summary_var",
                         label = "Select variables",
                         choices = list(),
                         multiple = TRUE,
                         selected = NULL),  
                         p("Select two variables to view the barplot."),
                         actionButton(inputId = "go_summary", "Plot"),
                       ), #close sidebarPanel
           mainPanel(width = 10,
                    fluidRow(
                      column(width = 4,
                         box(withSpinner(DT::dataTableOutput(outputId = "summary_table"), 
                                    type = 5, color = "#C91004"),
                         width = NULL),
                      ), #close column
                      
                    column(width = 8,
                           box(withSpinner(
                             plotlyOutput(outputId = "summary_barplot", height = 500), 
                             type = 5, color = "#C91004"),
                               width = NULL),
                    ) #close column
                    ) #close fluidrow
           )#close mainPanel
         )#close sidebarlayout
       ) #close summary
) #close tab_QC





