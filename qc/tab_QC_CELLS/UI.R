
##-------------------------------------##
##            CELL QC TAB              ##
##-------------------------------------##

tab_QC_CELLS <- tabItem(
  tabName = "Quality Control",
  textOutput(outputId = "session_id"),
  br(),br(),
  actionButton(inputId = "save_session", "Save Session"),
  actionButton(inputId = "flag_cells", "Flag cells"),
  actionButton(inputId = "remove_cells", "Remove cells"),
  br(),br(),
  sidebarLayout(
   sidebarPanel(width = 3,

########## Library Size and Features ########### 
      h5("Scater plots:"),
      radioButtons( inputId = "selectScale",
                        label = "Select scale",
                        choices = c("Linear", "Log10", "Log2"),
                        selected = "Linear",
                        inline = TRUE),

      sliderInput(inputId = "selectLibSizes", label = "Library sizes:",
            min = -1, max = 1000,
            value = c(-1,1000)),

      sliderInput(inputId = "selectLibBins", label = "Histogram bins:",
            min = 5, max = 100,
            value = 50),

      sliderInput(inputId = "selectNrFeatures", label = "Number of features:",
            min = 1, max = 1000,
            value = c(1,1000)),
      
      sliderInput(inputId = "selectFeatBins", label = "Histogram bins:",
            min = 5, max = 100,
            value = 50),
      sliderInput(inputId = "selectMtCutoff", label = "% MT-genes cutoff:",
            min = 0, max = 100,
            value = 50),

      h4("Subset cells:"),
      selectInput("features_sample",
            label = "Select condition",
            choices = list(),
            selected = NULL),
     
      selectInput("features_subsample",
                  label = "",
                  choices = list(),
                  selected = NULL,
                  multiple = TRUE),

      radioButtons(inputId = "violin_qc_showDots", choices = c("TRUE", "FALSE"), label = "Show dots on violins:", selected ="FALSE"),

               ), #close sidebarPanel

navset_card_underline(
  nav_panel("Scatter plots",
  mainPanel(width = 12,
            fluidRow(
              column(width = 6,
                     box(withSpinner(plotlyOutput(outputId = "histLibSize", height = 500), type = 8),
                         width = NULL)), #close column
                     column(width = 6,
                            box(withSpinner(plotlyOutput(outputId = "histFeatures", height = 500), type = 8),
                                width = NULL)) #close column
                     
                     
              
              ), #close fluidRow 
            
            fluidRow(
              column(width = 6,
                     box(withSpinner(plotlyOutput(outputId = "scatterLibFeat", height = 500), type = 8),
                         width = NULL)), #close column
              column(width = 6,
                     box(withSpinner(plotlyOutput(outputId = "scatterLibMt", height = 500), type = 8),
                         width = NULL)) #close column
              
            ), #close fluidRow 
            
          
            

  ) #close mainpanel
  ), #close nav_panel scatter plots
  nav_panel("Violin plots",
   mainPanel(width = 12,
             actionButton(inputId = "update_violins", "Update"),
  fluidRow(
    column(width = 6,
           box(withSpinner(plotlyOutput(outputId = "boxplot_libraries"), type = 8))   
    ), #close column
    #column(width = 1),
    column(width = 6,
           box(withSpinner(plotlyOutput(outputId = "boxplot_features"), type = 8))
    ) #close column
  ), #close fluidRow
    #column(width = 1),
  br(),
  br(),
  br(),
  br(),
  fluidRow(
    column(width = 6,
           box(withSpinner(plotlyOutput(outputId = "boxplot_mtgenes"), type = 8))
    ) #close fluidRow
    ) #close column
  
   ) #close mainpanel
  ), #close nav_panel violin plots  
) #close navset_card_underline

  ) #close sidebarLayout
) #close tab_QC





