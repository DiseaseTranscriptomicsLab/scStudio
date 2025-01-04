#INFORMATION-----------------------------

#LOAD LIBRARIES ------------------------
library(data.table)
library(dplyr)
library(DT)
library(ff)
library(fgsea)
library(GEOquery)
library(ggheatmap) #install dev github version
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(gridExtra)
library(limma)
library(Matrix)
library(matrixStats)
library(magrittr)
library(msigdbr)
library(plotly)
library(reshape2)
library(scater)
library(scran)
library(scuttle)
library(Seurat)
library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(shinyjs)
library(shinythemes)
library(stringi)
library(stringr)
library(tools)
library(vroom)

print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_REPORT/UI.R")

source("tab_SESSIONINFO/UI.R")

print("Successfully loaded tabs.")

#USER INFERFACE ----------------------------

ui <- fluidPage(
## theme ----------------------------------------

  useShinyjs(),
  theme = shinytheme("flatly"),
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  tags$script(src="https://cdn.plot.ly/plotly-latest.min.js"),
  
tags$head(
  tags$style(HTML(
    'body {
        font-size: 20px;
     }
     .irs-grid-text {
        font-size: 10pt;
     }
      button, .btn {
       font-size: 16px;
    }'
  ))
),

## navbarPage --------------------------------------- 
  navbarPage("scStudio",
          tabPanel("Report",
                     tab_REPORT),
          tabPanel("Session Information",
                   tab_SESSIONINFO)
      )#close navbarPage 
  
) #close UI 

# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  #output$dotplot <- renderDataTable(data.frame())

  }, silent = FALSE)
  
# OVERALL SESSION VARIABLES -----------------------------------
  overall_vars <- reactiveValues(metadata = data.frame())
  
# timer  ---------------------------------------------
autoInvalidate <- reactiveTimer(30000)  
# page configuration -----------------------------------  
output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)})

# UPLOAD SESSION ------------------------------------
## main modal -----------------------------------  
   uploadModal <- function(failed = FALSE) {
     modalDialog(
       h1("Welcome to scStudio - Report"),
       p("Please load a session to begin."),
       h5("Option 1"),
       textInput("upload_token", "Insert token",
                 placeholder = 'E.g.: exzy8sl9zc'),
       h5("Option 2"),
       fileInput(inputId = 'upload_session_zip', label = 'Upload token.zip', multiple = FALSE),
       span(""),
       if (failed)
         div(tags$b("Invalid token or file", style = "color: red;")),
       
       footer = tagList(
         actionButton("ok_upload_session", "Upload session"))
       )} #close dataModal
   
   observe({
     showModal(uploadModal())
   }) # close upload session

## ok_upload_session --------------------------------   
   observeEvent(input$ok_upload_session, {
### zip modality -------------------------------------    
     if (!is.null(input$upload_session_zip)) {
       
       removeModal()
       print(input$upload_session_zip)
       
       token <- gsub("\\.zip","",input$upload_session_zip$name)
       unzip(zipfile = input$upload_session_zip$datapath,
             exdir = paste0(getwd(),"/tokens/", "temp"))
       old <- paste0(getwd(),"tokens/temp","/tokens/",token)
    
       new <- paste0(getwd(),"/tokens/",token)

       file.move(old, new)
       unlink(getwd(),"/tokens/temp", recursive = TRUE)
       
       overall_vars$session_token <- token
       
       withProgress(message = 'Uploading data', value = 0, {
         
       files <- list.files(new)
         
       if ("countMatrices.rds" %in% files & "metadata.rds" %in% files){
         
      incProgress(0.1, detail = "")
            
             incProgress(0.6, detail = "Finishing...")
             
             if ("report.rds" %in% files){
               overall_vars$report<- upload_session(overall_vars$session_token, "report")
             } #close if
             else {overall_vars$report <- list()}
             overall_vars$md5$report <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/report.rds"))
             
        } #close if 
       
       else{
         showNotification("Session is not valid.")
       } #close else
         incProgress(1, detail = "Done.")
    }) #close progress
} #close upload zip
     
## token modality -------------------------------------      
     # Check that token exists
     else if (!is.null(input$upload_token)) {
       removeModal()
       
       upload_token_in_dir <- input$upload_token %in% list.dirs(path = "./tokens/", full.names = FALSE)[-1]
       
       if (!upload_token_in_dir ) {showModal(uploadModal(failed = TRUE))}
       
       else{
         removeModal()
         
         overall_vars$session_token <- input$upload_token
           
         files <- list.files(paste0(getwd(),"/tokens/",input$upload_token))
         
         if ("countMatrices.rds" %in% files & "metadata.rds" %in% files){
           
           withProgress(message = 'Uploading data', value = 0, {
             
           incProgress(0.5, detail = "Retrieving metadata...")
           
           
           incProgress(0.6, detail = "Finishing...")

           
           if ("report.rds" %in% files){
             overall_vars$report <- upload_session(overall_vars$session_token, "report")
           } #close if 
          
            else{overall_vars$report <- list()}
           
           overall_vars$md5$report <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/report.rds"))
           
          incProgress(1, detail = "Done.")
         }) #close progress
      } #close if 
           
         else{
           showNotification("Session is not valid.")
         } #close else
         
       } #close else
  } #close first if 
  
}) #close upload session
   

# UPDATE INPUT OPTIONS/PLOTS ---------------------------------------------------
  observe({
    
## md5sum checks --------------------------------------------------------------
    
    try({
      # check if files changed
      check_md5sum_report <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/report.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_report != overall_vars$md5$report){
        overall_vars$report <- upload_session(overall_vars$session_token, "report")
        overall_vars$md5$report <- check_md5sum_report
       } # close if
      
    }) #close try
  
}) #close observe
  
## update inputs ----------------------------------------------------   
observe({
  
updateSelectInput(session, "select_analysis_type",
                   choices = c("QC", "NORMALIZATION", "BEC",
                               "PCA","TSNE","UMAP","Clustering", "DEA", "SCORES", "GSEA"), 
                   selected = "")
  
  output$session_info <- renderText({
    paste(capture.output(sessionInfo()), collapse = "\n")  # Capture and format
  })

})
   
   observeEvent(input$select_analysis_type, {
     choices <- c()
     updateSelectInput(session, "select_analysis_id",
                       choices = "")
     
     #if (input$select_analysis_type == "QC" ) {choices = "Main"}
     if (input$select_analysis_type == "NORMALIZATION" ){choices = names(overall_vars$report$normalization)}
     if (input$select_analysis_type == "BEC" ) {choices = names(overall_vars$report$bec)}
     if (input$select_analysis_type == "PCA" ) {choices = names(overall_vars$report$pca)}
     if (input$select_analysis_type == "TSNE" ) {choices = names(overall_vars$report$tsne)}
     if (input$select_analysis_type == "UMAP" ) {choices = names(overall_vars$report$umap)}
     if (input$select_analysis_type == "Clustering" ) {choices = names(overall_vars$report$clustering)}
     if (input$select_analysis_type == "DEA" ) {choices = names(overall_vars$report$DEA)}
     if (input$select_analysis_type == "SCORES" ) {choices = names(overall_vars$report$scores)}
     if (input$select_analysis_type == "GSEA" ) {choices = names(overall_vars$report$gsea)}
     
     updateSelectInput(session, "select_analysis_id",
                       choices = choices)
     
     if (input$select_analysis_type == "QC" ){
       
       print("Here")
       
       text <- paste0("Total number of discarded cells: ", 
               as.numeric(table(overall_vars$report$QC$filter_cells))[1], "\n",
               "Scale: ", overall_vars$report$QC$scale, "\n",
               "Library size >= ", overall_vars$report$QC$library_size_min, "\n",
               "Library size <= ", overall_vars$report$QC$library_size_max, "\n",
               "Total unique genes >= ", overall_vars$report$QC$total_features_size_min, "\n",
               "Total unique genes <= ", overall_vars$report$QC$total_features_size_max, "\n",
               "% counts attributed to mitochondrial genes >= ", 
               overall_vars$report$QC$percentage_mt_genes, "\n",
               "\n",
               "Total number of discarded genes: ", 
               as.numeric(table(overall_vars$report$QC$filter_genes))[1], "\n",
               "Total counts across cells > ", overall_vars$report$QC$total_gene_expression, "\n",
               "Number of cells expressing the gene > ", overall_vars$report$QC$number_expressing_cells) 
     }
     print(text)
     output$report <- renderText({ text })
   })
   
   observeEvent(input$select_analysis_id, {
     
     if (input$select_analysis_type == "DEA" ){
       
       print("selected DEA analyses")
       
       id <- gsub("_MAST$","",input$select_analysis_id)
       
       print(id)
       
       id <- gsub("_t$","",id)
       id <- gsub("_wilcox$","",id)
       id <- gsub("_roc$","",id)
       id <- gsub("_bimod$","",id)
       print(id)
       print(overall_vars$report$DEA[[id]]$matrix)
       
       text <- paste0("Matrix: ",overall_vars$report$DEA[[id]]$matrix, "\n",
                      "Methods: ", overall_vars$report$DEA[[id]]$method,"\n",
                      "Groups: ", overall_vars$report$DEA[[id]]$group1," (", 
                      overall_vars$report$DEA[[id]]$group1_type,")",
                      " vs ", 
                      overall_vars$report$DEA[[id]]$group2," (",
                      overall_vars$report$DEA[[id]]$group2_type,")","\n")
     }
     
     if (input$select_analysis_type == "SCORES" ){
       id <- input$select_analysis_id
       text <- paste0("Matrix: ", overall_vars$report$scores[[id]]$matrix, "\n",
                      "Organism: ", overall_vars$report$scores[[id]]$organism, "\n",
                      "Category: ", overall_vars$report$scores[[id]]$category, "\n",
                      "Genes: ", "\n", toString(overall_vars$report$scores[[id]]$genes), "\n") 
     }
     
     if (input$select_analysis_type == "NORMALIZATION" ){
       id <- input$select_analysis_id
       text <- paste0("Method: ", overall_vars$report$normalization[[id]]$method, "\n") 
     }
     
     if (input$select_analysis_type == "BEC" ){
       id <- input$select_analysis_id
       text <- paste0("Matrix: ", overall_vars$report$bec[[id]]$bec, "\n",
                      "Method: ", overall_vars$report$bec[[id]]$method, "\n",
                      "Batch label: ", overall_vars$report$scores[[id]]$batch_label1, "\n",
                      "Batch label 2: ", overall_vars$report$scores[[id]]$batch_label2, "\n",
                      "Biological condition: ", overall_vars$report$scores[[id]]$biological_condition, "\n") 
     }
     
     if (input$select_analysis_type == "PCA" ){
       id <- input$select_analysis_id
       print(overall_vars$report$pca[[id]])
       text <- paste0("Matrix: ", overall_vars$report$pca[[id]]$matrix, "\n",
                      "Components: ", overall_vars$report$pca[[id]]$ncomponents, "\n",
                      "Features: ", overall_vars$report$pca[[id]]$features, "\n",
                      "Scaling: ", overall_vars$report$pca[[id]]$scaling, "\n") 
     }
     
     if (input$select_analysis_type == "TSNE" ){
       id <- input$select_analysis_id
       if (overall_vars$report$tsne[[id]]$method == "PCA"){
         print(overall_vars$report$tsne[[id]]$matrix)
         text <- paste0("PCA: ", overall_vars$report$tsne[[id]]$matrix, "\n",
                        "Components: ", overall_vars$report$tsne[[id]]$features, "\n",
                        "Perplexity: ", overall_vars$report$tsne[[id]]$perplexity, "\n") 
       }
       else {
         text <- paste0("Matrix: ", overall_vars$report$tsne[[id]]$matrix, "\n",
                        "Features: ", overall_vars$report$tsne[[id]]$features, "\n",
                        "Perplexity: ", overall_vars$report$tsne[[id]]$perplexity, "\n")
       }
     }
     
     if (input$select_analysis_type == "UMAP" ){
       id <- input$select_analysis_id
       if (overall_vars$report$umap[[id]]$method == "PCA"){
         text <- paste0("PCA: ", overall_vars$report$umap[[id]]$matrix, "\n",
                        "Components: ", overall_vars$report$umap[[id]]$features, "\n",
                        "Minimum distance: ", overall_vars$report$umap[[id]]$min_dist, "\n",
                        "#neighbors: ", overall_vars$report$umap[[id]]$nneigh, "\n") 
       }
       else {
         text <- paste0("Matrix: ", overall_vars$report$umap[[id]]$matrix, "\n",
                        "Features: ", overall_vars$report$umap[[id]]$features, "\n",
                        "Minimum distance: ", overall_vars$report$umap[[id]]$min_dist, "\n",
                        "#neighbors: ", overall_vars$report$umap[[id]]$nneigh, "\n")
       }
     }
     
     
     if (input$select_analysis_type == "Clustering" ){

       id <- input$select_analysis_id
       if (overall_vars$report$clustering[[id]]$method == "PCA"){
         text <- paste0("PCA: ", overall_vars$report$clustering[[id]]$matrix, "\n",
                        "Components: ", overall_vars$report$clustering[[id]]$features, "\n") 
       }
       else {
    
         text <- paste0("Matrix: ", overall_vars$report$tsne[[id]]$matrix, "\n",
                        "Features: ", overall_vars$report$tsne[[id]]$features, "\n")
       }
       output$report <- renderText({ text })
     }
     
     if (input$select_analysis_type == "GSEA" ){
       id <- input$select_analysis_id
       print(overall_vars$report$gsea[[id]])
       text <- paste0("Rank: ", overall_vars$report$gsea[[id]]$rank, "\n",
                      "Ordering: ", overall_vars$report$gsea[[id]]$ordering, "\n",
                      "Organism: ", overall_vars$report$gsea[[id]]$organism, "\n",
                      "Categories: ", overall_vars$report$gsea[[id]]$subcategories, "\n") 
     }
     
     
     output$report <- renderText({ text })
   })
   


}#close server
  
  
 
 shinyApp(ui = ui, server = server)


