
##-------------------------------------##
##           Pseudobulk TAB            ##
##-------------------------------------##

tab_PSEUDO <- tabItem(
  tabName = "Pseudobulk",
  textOutput(outputId = "session_id"),
   sidebarLayout(
     sidebarPanel(width = 4,
                  
                  h3("Convert to pseudobulk:"),
                  
                  selectInput(inputId = "select_matrix_pseudo", 
                              label = "Select matrix", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  selectInput(inputId = "select_group_pseudo", 
                              label = "Select groups", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  selectInput(inputId = "select_condition_pseudo", 
                              label = "Select condition", 
                              choices = c("-"), 
                              multiple = FALSE, 
                              width = "100%"),
                  
                  h3("Select gene sets:"),
                  selectInput(inputId = "gset_scores_organism_pseudo", 
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
                  
                  selectInput(inputId = "select_database_pseudo", 
                              label = "Pathways:", 
                              choices = c(#"Custom", 
                                          "H",
                                          "C1",
                                          "CGP",
                                          "CP",
                                          "CP:BIOCARTA",
                                          "CP:KEGG_LEGACY",
                                          "CP:KEGG_MEDICUS",
                                          "CP:PID",
                                          "CP:REACTOME",
                                          "CP:WIKIPATHWAYS",
                                          "MIR:MIRDB",
                                          "MIR:MIR_Legacy",
                                          "TFT:GTRD",
                                          "TFT:TFT_Legacy",
                                          "3CA",
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
                              selected = "Custom"), 
                                     
                  uiOutput("gset_pathway_pseudo"),
                                    
                  #uiOutput("gene_set_pseudo"),
                  
                  h3("Select scoring method:"),
                  
                  selectInput(inputId = "method_pseudobulk", 
                              label = "Method:", 
                              choices = c("logmedian", 
                                          "ranking",
                                          "ssGSEA"), 
                              multiple = FALSE,
                              selected = "ranking"),
                  
                  actionButton(inputId = "calculate_pseudobulk", 
                               label = "Calculate")
                  
                                     
                  
), #close sidebarPanel
                        
     mainPanel(
       fluidRow(
         column(width = 10,
                box(withSpinner(
                  plotOutput(outputId = "marker"), 
                  type = 5),
                  width = NULL,
                  height = "700px") #, height = 150
                ) #close column
         ), #close fluidrow
       fluidRow(
         column(width = 10,
                box(withSpinner(
                  plotOutput(outputId = "marker_roc"), 
                  type = 5),
                  width = NULL,
                  height = "700px") #, height = 150
         ) #close column
       ) #close fluidrow
       ) #close mainPanel
     ) #close sidebarlayout
) #close FEA

