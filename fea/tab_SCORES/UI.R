
##-------------------------------------##
##              SCORES TAB             ##
##-------------------------------------##

tab_SCORES <- tabItem(
  tabName = "Scores",
  textOutput(outputId = "session_id"),
   sidebarLayout(
     sidebarPanel(width = 4,
                  
                  h4("Calculate new score:"),
                  
                  selectInput(inputId = "select_matrix_scores", 
                              label = "Select matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  textInput(inputId = "gset_id", label = "Gene set ID:"),
                  
                  selectInput(inputId = "gset_scores_organism", 
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
                  
                  selectInput(inputId = "select_database", 
                              label = "Pathways:", 
                              choices = c("Custom", 
                                          "H",
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
                              selected = "H"), 
                                     
                  uiOutput("gset_pathway"),
                                    
                  uiOutput("gene_set"),
                  
                  actionButton(inputId = "calculate_gset_score", 
                               label = "Calculate"),
                  
                  h4("Plotting"),
                  
                  selectInput(inputId = "select_gset_scores", 
                              label = "Gene set scores:", 
                              choices = c("-"), 
                              multiple = TRUE), 
                                     
                  p("Choose at least 2 scores to view heatmap."),
                                     
                  selectInput(inputId = "gset_var", 
                              label = "Colour by:", 
                              choices = c("-")),
                                     
                  radioButtons( inputId = "scores_scale_heatmap",
                                label = "Scale:",
                                choices = c("row", "column", "none"),
                                selected = "none",
                                inline = TRUE),
                                     
                  radioButtons( inputId = "cluster_scores",
                                label = "Cluster:",
                                choices = c("row", "column", "both", "none"),
                                selected = "none",
                                inline = TRUE),
                  
                  actionButton(inputId = "plot_gset_scores", label = "Plot heatmap"),
                  ), #close sidebarPanel
                        
     mainPanel(
       fluidRow(
         column(width = 10,
                box(withSpinner(
                  plotOutput(outputId = "gset_heatmap"), 
                  type = 5),
                  width = NULL) #, height = 150
                ) #close column
         ) #close fluidrow
       ) #close mainPanel
     ) #close sidebarlayout
) #close FEA
















































