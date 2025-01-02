
##-------------------------------------##
##         SELECTDATA TAB              ##
##-------------------------------------##
get_dataset_names <- function(){
  files <- list.dirs("./public_datasets")
  files <- gsub("./public_datasets/", "", files)
  files <- files[-1]
  return(files)
}


tab_SELECTDATA<- tabItem(
  tabName = "Select Data",
  textOutput(outputId = "session_id"),
  br(),br(),
  actionButton(inputId = "upload_session", "Upload Session"),
  actionButton(inputId = "save_session", "Save Session"),
  actionButton(inputId = "publish_analysis", "Publish Analysis"),
  
  
  br(),br(),
  
sidebarLayout(
  sidebarPanel(width = 4,
               ########## Upload GEO ###########                
               h3("Gene Expression Omnibus:"), 
               p("Explore single-cell RNA sequencing datasets that \
                 have been openly shared on the NCBI's GEO database, accessible at https://www.ncbi.nlm.nih.gov/geo/."),
               textInput(inputId = "GEO_id", label = "GEO ID", value = "GSE181878"),
               actionButton(inputId = "view_GEO", "Load"),
               br(),br(),
########## Pre-processed datasets ###########                 
               h4("Public analyses:"),  
               p("Explore datasets and analyses shared by fellow scStudio users."), 
               selectInput(inputId = "selectData", label = "Dataset", 
                           choices = get_dataset_names()),
               actionButton(inputId = "view_preprocessed", "Load"),

br(),br(),

########## Upload dataset ###########                 
               h4("Upload data:"), 

               p("Upload your own dataset."),

               textInput(inputId = "upload_id", label = "Name:", value = "matrix1"),

               radioButtons(inputId = "select_file_type", 
                            label = "Select file type: ", 
                            choices = c("table (.csv/.tsv/.txt)", 
                                        "Seurat (.rds)",
                                        "SingleCellExperiment (.rds)",
                                        "10X Genomics (triplet format)",
                                        "10X Genomics (.h5)"),
                            selected = "Seurat (.rds)"),

               uiOutput("matrix_upload"),

               

               uiOutput("seurat_upload"),
            
               uiOutput("sce_upload"),

               uiOutput("matrix_10X"),
               uiOutput("barcodes_10X"),
uiOutput("features_10X"),

               uiOutput("h5_10X"),

uiOutput("metadata_upload"),

               actionButton(inputId = "view_upload", "Load"), 

               br(),br(),br(),
               actionButton(inputId = "delete_session", "Delete Session",
                            style="border-color: #C91004"),
               
               ), #close sidebarPanel
  mainPanel(width = 8,
            
            selectInput(inputId = "availableCountMatrices", label = "Available count matrices in session:", 
                        choices = c("-"), multiple = TRUE, width = "100%"),
            
            actionButton(inputId = "delete_countMatrix", "Delete", 
                         style="color: #fff; background-color: #C91004; border-color: #C91004"),
            br(), br(),
            dataTableOutput(outputId = "countMatricesProperties"),
            br(), br(),
            
            selectInput(inputId = "selectCountMatrices", label = "Add samples:", 
                        choices = c("-"), multiple = TRUE, width = "100%"),
            p("If you have multiple samples or data files you wish to analyze collectively,\
              simply choose the desired ones and merge them into a single raw count matrix."),
            actionButton(inputId = "kickoff", "Add samples", 
                         style="color: #fff; background-color: #C91004; border-color: #C91004"),
            br(), br(),
            dataTableOutput(outputId = "sampleProperties"),
            
            selectInput(inputId = "subsetMetadata", label = "Sample information:", 
                        choices = c("-"), multiple = TRUE, width = "100%"),
            p("Use this filter to select the specific cell annotations of interest.\
              Keep in mind that once you save the session, any excluded information will be irretrievable."),
            box(dataTableOutput(outputId = "metadata_preview"), width = 12)
            
            ) #close mainPanel
    ) #close sidebarLayout
) #close tab_SELECTDATA





