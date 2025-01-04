
##-------------------------------------##
##              GSEA TAB               ##
##-------------------------------------##

tab_GSEA <- tabItem(
  tabName = "GSEA",
  sidebarLayout(
    sidebarPanel(width = 2,
                 
                 h4("Gene Set Enrichment Analysis"),
                 
                 selectInput(inputId = "select_dea_gsea", 
                             label = "Select DEA", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 selectInput(inputId = "gsea_ordering", 
                             label = "Select ordering:", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 selectInput(inputId = "gsea_organism", 
                             label = "Organism:", 
                             choices = c("Homo sapiens",
                                         "Mus musculus",
                                         "Danio rerio",
                                         "Rattus norvegicus",
                                         "Macaca mulatta",
                                         "Drosophila melanogaster",
                                         "Caenorhabditis elegans",
                                         "Saccharomyces cerevisiae",
                                         "Anolis carolinensis",
                                         "Bos taurus",
                                         "Canis lupus familiaris",
                                         "Equus caballus",
                                         "Felis catus",
                                         "Gallus gallus",
                                         "Monodelphis domestica",
                                         "Ornithorhynchus anatinus",
                                         "Pan troglodytes",
                                         "Schizosaccharomyces pombe 972h-",
                                         "Sus scrofa",
                                         "Xenopus tropicalis"
                                         ), multiple = FALSE, width = "100%"),
                 
                 selectInput(inputId = "gsea_pathways", 
                             label = "Pathways:", 
                             choices = c("H",
                                         "C1",
                                         "CGP",
                                         "CP",
                                         "CP:BIOCARTA",
                                         "CP:KEGG",
                                         "CP:PID",
                                         "CP:REACTOME",
                                         "CP:WIKIPATHWAYS",
                                         "MIR:MIRDB",
                                         "MIR:MIR_Legacy",
                                         "TFT:GTRD",
                                         "TFT:TFT_Legacy",
                                         "CGN",
                                         "CM",
                                         "GO:BP",
                                         "GO:CC",
                                         "GO:MF",
                                         "HPO",
                                         "C6",
                                         "IMMUNESIGDB",
                                         "VAX",
                                         "C8"), 
                             multiple = TRUE, 
                             width = "100%"),
                 
                 textInput(inputId = "gsea_id", 
                           label = "Identifier:",
                           value = "gsea_1"),
                 
                 actionButton(inputId = "run_gsea", "Run GSEA"),
                 
                 h4("Plotting"),
                 
                 selectInput(inputId = "select_gsea", 
                             label = "Select GSEA", 
                             choices = c("-"), 
                             multiple = FALSE, 
                             width = "100%"),
                 
                 textInput(inputId = "gsea_pval", 
                           label = "Adjusted p-value cutoff ", 
                           value = 0.01),
                 
                 textInput(inputId = "gsea_NES", 
                           label = "NES cutoff ", 
                           value = 1.5),
                 
                 actionButton(inputId = "go_plot_gsea", 
                              label = "Plot results")
                 ), #close sidebarPanel
    mainPanel(
      fluidRow(
        column(width = 12,
               box(withSpinner(
                 plotOutput(outputId = "summary_gsea", 
                            height = 500), 
                 type = 5),
                 width = NULL)) #close column
      ), #close fluidRow
      br(),br(), br(), br(), br(), br(), br(), br(),
        fluidRow(
        column(width = 9,
               box(withSpinner(
                 DT::dataTableOutput(outputId = "table_gsea"), 
                 type = 5),
                 width = NULL)
               ) #close column
        ), #close fluidrow
      
      br(),br(),br(),
      
      fluidRow(
        column(width = 6,
               box(withSpinner(
                 plotOutput(outputId = "pathway_gsea"), 
                 type = 1), 
                 width = NULL)
               ), #close column
        column(width = 6,
               box(withSpinner(
                 textOutput(outputId = "leading_edge"),
                 type = 7), 
                 width = NULL)
               ) #close column  
        ) #close fluidrow
      ) #close mainPanel    
    ) #close sidebarLayout
) #close FEA
















































