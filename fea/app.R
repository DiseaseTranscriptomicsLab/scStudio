#INFORMATION-----------------------------


#LOAD LIBRARIES ------------------------
library(data.table)
library(DT)
library(dplyr)
library(ff)
library(fgsea)
library(ggheatmap) #install dev github version
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(gridExtra)
library(Matrix)
library(matrixStats)
library(magrittr)
library(msigdbr)
#library(msigdbdf)
library(plotly)
library(reshape2)
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


print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_SCORES/UI.R")
source("tab_SCORES/server.R")

source("tab_GSEA/UI.R")
source("tab_GSEA/server.R")

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
  navbarPage("scStudio: FEA",
          tabPanel("Scores",
                     tab_SCORES),
          tabPanel("GSEA",
                   tab_GSEA)
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$gset_heatmap <- renderDataTable(data.frame())
  output$table_gsea<- renderDataTable(data.frame())
  output$summary_gsea <- renderDataTable(data.frame())
  output$pathway_gsea<- renderDataTable(data.frame())



  }, silent = TRUE)
  
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
       h1("Welcome to scStudio - Functional Enrichment Analysis"),
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
     print("BUTTON PUSHED")
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
             
             incProgress(0.5, detail = "Retrieving metadata...")
             
             overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
             overall_vars$md5$metadata <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
             
             overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
             overall_vars$md5$mat_names <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
             
             overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
             overall_vars$md5$genes <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
             
             incProgress(0.6, detail = "Finishing...")
             
             if ("colours.rds" %in% files){
               overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
             } #close if 
             
             else{overall_vars$colours <- list()}
             
             if ("report.rds" %in% files){
               overall_vars$report<- upload_session(overall_vars$session_token, "report")
             } #close if
             else {overall_vars$report <- list()}
             
             if ("dea.rds" %in% files){
               overall_vars$dea <- upload_session(overall_vars$session_token, "dea")
             } #close if 
             else {overall_vars$dea <- list()
             save_session(overall_vars$session_token ,overall_vars, "dea")}
             
             
             overall_vars$md5$dea <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
             
             if ("scores.rds" %in% files){
               overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
             } #close if 
             else {overall_vars$scores <- list()
             save_session(overall_vars$session_token ,overall_vars, "scores")}
             
             if ("gsea.rds" %in% files){
               overall_vars$gsea <- upload_session(overall_vars$session_token, "gsea")
             } #close if 
             else {overall_vars$gsea <- list()
             save_session(overall_vars$session_token ,overall_vars, "gsea")}
  
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
           
           overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
           
           overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
           
           overall_vars$md5$metadata <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
           
           overall_vars$md5$mat_names <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
           
           overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
           overall_vars$md5$genes <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
           
           incProgress(0.6, detail = "Finishing...")
           
           if ("colours.rds" %in% files){
             overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
           } #close if 
          
            else{overall_vars$colours <- list()}
           
           if ("report.rds" %in% files){
             overall_vars$report <- upload_session(overall_vars$session_token, "report")
           } #close if 
          
            else{overall_vars$report <- list()}
           
           if ("dea.rds" %in% files){
             overall_vars$dea <- upload_session(overall_vars$session_token, "dea")
           } #close if 
           else {overall_vars$dea <- list()
           save_session(overall_vars$session_token ,overall_vars, "dea")}
           
           overall_vars$md5$dea <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
           
           
           if ("scores.rds" %in% files){
             overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
           } #close if 
           else {overall_vars$scores <- list()
           save_session(overall_vars$session_token ,overall_vars, "scores")}
           
           if ("gsea.rds" %in% files){
             overall_vars$gsea <- upload_session(overall_vars$session_token, "gsea")
           } #close if 
           else {overall_vars$gsea <- list()
           save_session(overall_vars$session_token ,overall_vars, "gsea")}
        
          incProgress(1, detail = "Done.")
         }) #close progress
      } #close if 
           
         else{
           showNotification("Session is not valid.")
         } #close else
         
       } #close else
  } #close first if 
  
}) #close upload session
   

# SAVE SESSION -------------------------------------------
## main modal --------------------------------------------   
   saveModal <- function(failed = FALSE) {
     modalDialog(
       h3("Session: ", overall_vars$session_token),
       actionButton("ok_current_session", "Save current session"),
       actionButton("ok_create_session", "Save new session"),
       downloadButton("ok_download", "Download"),
       
       span(""),
       if (failed)
         div(tags$b("Invalid token name", style = "color: red;")),
       
       footer = tagList(
         modalButton("Cancel")
       )
     )} #close dataModal
   
   observeEvent(input$save_session,{
     showModal(saveModal())
   }) # close save session
   
## ok_download ------------------------------------------
   output$ok_download <- downloadHandler(
     filename = function() {
       paste0(overall_vars$session_token, ".zip")
     },
     content = function(file) {
       withProgress(message = 'Downloading latest saved session...', value = 0, {
         incProgress(0.4, detail = "")
         on.exit(removeModal())
         # Zip the selected files
         zip(file, paste0("tokens/", overall_vars$session_token))
       }) # close progress bar
     }
   )
   
## ok_create_session -------------------------------------      
   observeEvent(input$ok_create_session, {
     
     new_token <- stri_rand_strings(1, 6, pattern = "[a-z0-9]")
     
     removeModal()
     
     showModal(modalDialog(
       title = "New token generated", #if token is null, we create a new one
       paste("Use the following token to retrieve your session:",new_token)
     )) #close showModal
     
     withProgress(message = 'Saving data', value = 0, {
       incProgress(0.1, detail = "")
       
       save_session(new_token,overall_vars, c("metadata"))
       overall_vars$md5$metadata <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       
       save_session(new_token,overall_vars, c("gene_qc"))
       
       #COPY old TOKEN folder
       incProgress(0.9, detail = "Finishing...")
       
     }) #close progress bar
   }) # close create new session
   
## ok_current_session --------------------------------------   
   observeEvent(input$ok_current_session, {
     
     removeModal()
     
     withProgress(message = 'Saving data', value = 0, {
       
       token <- overall_vars$session_token
       
       incProgress(0.1, detail = "")
       
       save_session(token, overall_vars, c("metadata"))
       overall_vars$md5$metadata <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       
       save_session(token,overall_vars, c("gene_qc"))
       
       incProgress(0.9, detail = "Finishing...")
       
       showNotification("Session saved successfully.")
       
     }) #close progress bar
   }) # close create new session   
   
# UPDATE INPUT OPTIONS/PLOTS ---------------------------------------------------
  observeEvent(input$select_dea_gsea, {

      if (length(grep("roc",input$select_dea_gsea)) == 1){
        
        opts <- c("AUC", "Power", "Log2FC (mean)", "Log2FC (median)", "SNR", "|SNR|")
      }
      else { opts <- c("Log2FC (mean)", "Log2FC (median)", "Adjusted p-value", "SNR", "|SNR|")}
      
      updateSelectInput(session, "gsea_ordering",choices = opts)
  })
   
   
   observe({
    
## running jobs --------------------------------------------      
    autoInvalidate()
    print(Sys.time())
    print("Running jobs")
    print(overall_vars$jobs)
    print("Job result")
    
    for (job in names(overall_vars$jobs)){
      print(job)
      try(print(overall_vars$jobs[[job]]$get_result()), silent = FALSE)
    }
    

    
## md5sum checks --------------------------------------------------------------
    
    try({
      # check if files changed
      check_md5sum_metadata <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
      
      check_md5sum_mat_names <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
      
      check_md5sum_dea <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
      } # close if
      
      if (check_md5sum_dea != overall_vars$md5$dea){
        overall_vars$dea <- upload_session(overall_vars$session_token, "dea")
        
        
      } # close if
      
    }) #close try

}) #close observe
  

   
## update inputs ----------------------------------------------------   
observe({

  if (input$select_database != "Custom"){
    
    output$gene_set <- NULL
    
    get_categories <-  readRDS(paste0(getwd(), "/","/msigdbr_collections.rds"))

    if (input$select_database %in% c("H", "C1", "C6", "C8")){  
      options <- msigdbr(species = input$gset_scores_organism, category = input$select_database)
      options <- unique(options$gs_name)
    }
    else{
      cat <- get_categories[get_categories$gs_subcat == input$select_database, "gs_cat"]
      options <- msigdbr(species = input$gset_scores_organism, category = cat, subcategory = input$select_database)
      options <- unique(options$gs_name)
    }
    
    output$gset_pathway <- renderUI(selectInput(inputId = "gset_pathway", label = "Gene Sets:", 
                                                choices = options, multiple = TRUE))

  } #close if 
  
  else {output$gene_set <- renderUI(
    selectizeInput(inputId = "gene_set", 
                   label = "Custom Gene Set:",
                   choices = NULL, 
                   selected = NULL, 
                   multiple = TRUE, 
                   options = NULL))
  
  updateSelectizeInput(session, 
                       "gene_set", 
                       choices =  overall_vars$genes, 
                       server = TRUE)
  
  output$gset_pathway <- NULL}
  
  updateSelectInput(session, "select_matrix_scores",
                    choices = overall_vars$mat_names[-1])
  
  updateSelectInput(session, "gset_var",
                    choices = names(identify_discrete(overall_vars$metadata)), 
                    selected = "orig.ident")
  
  updateSelectInput(session, "select_gset_scores",
                    choices = names(overall_vars$scores))
  
  updateSelectizeInput(session, 
                       "featplot_gsets", 
                       choices =  names(overall_vars$scores), 
                       server = TRUE)
  
  updateSelectizeInput(session, 
                       "featplot_violins_gsets", 
                       choices =  names(overall_vars$scores), 
                       server = TRUE)
  
  updateSelectizeInput(session, 
                       "featplot_dotplot_gsets", 
                       choices =  names(overall_vars$scores), 
                       server = TRUE)
  
  
  
  updateSelectInput(session, "select_dea_gsea",choices = names(overall_vars$dea))
  
  updateSelectInput(session, "select_gsea",choices = names(overall_vars$gsea))
  

}) #close observe
 
   
# RUN SCORES -----------------------------------------------  

observeEvent(input$calculate_gset_score,{
  
  withProgress(message = 'Calculating scores...', value = 0, {
    
    incProgress(0.2, detail = "")
    
    if(input$select_database == "Custom"){  
      
      print("Pushed button to calculate custom gene set score")
      
      incProgress(0.4, detail = "Preparing count matrix")
      overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
      
      try({
        
        gset_score <- get_score(mat = overall_vars$countMatrices[[input$select_matrix_scores]], 
                                gene_set = input$gene_set, 
                                name = input$gset_id)
        
        overall_vars$countMatrices <- NULL
                        
        overall_vars$scores[[input$gset_id]] <- gset_score
        
        save_session(overall_vars$session_token, overall_vars, "scores")
                        
        overall_vars$report$scores[[input$gset_id]][["genes"]] <- input$gene_set
        overall_vars$report$scores[[input$gset_id]][["matrix"]] <- input$select_matrix_scores
        overall_vars$report$scores[[input$gset_id]][["organism"]] <- "-"
        overall_vars$report$scores[[input$gset_id]][["category"]] <- "Custom"
        }, silent = FALSE) # close try 
      } #close if 
    else{
      incProgress(0.4, detail = "Preparing count matrix")
      overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
      
      print("Pushed button to calculate pathway gene set score")
      try({
        gset_score_list <- get_score_pathways(
          mat_name = input$select_matrix_scores,
          mat = overall_vars$countMatrices[[input$select_matrix_scores]], 
          organism = input$gset_scores_organism,
          database = input$select_database,
          pathways = input$gset_pathway,
          name = input$gset_id)
        
        overall_vars$countMatrices <- NULL
        
        overall_vars$scores <- c(overall_vars$scores, gset_score_list[[1]])
        save_session(overall_vars$session_token, overall_vars, "scores")
                        
        overall_vars$report$scores <- c(overall_vars$report$scores, gset_score_list[[2]])
        }) # close try 
      }  # close else 
    
    updateSelectInput(session, "select_gset_scores",
                      choices = names(overall_vars$scores),
                      selected = NULL)
    
    updateSelectizeInput(session, 
                         "gene_set", 
                         choices =  overall_vars$genes, 
                         server = TRUE)
    }) # close progress bar
  
  save_session(overall_vars$session_token, overall_vars, "report")                
  showNotification("Scores ready.")
  })
   
observeEvent(input$plot_gset_scores ,{
     
  if (length(input$select_gset_scores) > 1 ){
    

        
    output$gset_heatmap <- bindEvent(renderPlot({
      
      p <- plot_gset_heatmap(scores = overall_vars$scores, 
                        ids = input$select_gset_scores, 
                        cluster = input$cluster_scores, 
                        scale = input$scores_scale_heatmap, 
                        selected_cols = overall_vars$colours[[input$gset_var]], 
                        groups = overall_vars$metadata[[input$gset_var]])
      p
    }, height = 800, width = 1500), input$plot_gset_scores)
    }
  })
   
# RUN GSEA --------------------------------------------------------------------   
observeEvent(input$run_gsea, {
  
  print("Running GSEA")
  
  withProgress(message = 'Running GSEA...', value = 0, {
    
    incProgress(0.4, detail = "")
    
    try({
         overall_vars$gsea[[input$gsea_id]] <- get_gsea(
           dea = overall_vars$dea[[input$select_dea_gsea]], 
           ordering = input$gsea_ordering, 
           organism = input$gsea_organism, 
           subcategories = input$gsea_pathways)
         
         updateSelectInput(session, "select_gsea",choices = names(overall_vars$gsea))
         
         overall_vars$report$gsea[[input$gsea_id]]$rank <- input$select_dea_gsea
         overall_vars$report$gsea[[input$gsea_id]]$ordering <- input$gsea_ordering
         overall_vars$report$gsea[[input$gsea_id]]$organism <- input$gsea_organism
         overall_vars$report$gsea[[input$gsea_id]]$subcategories <- input$gsea_pathways
         
         save_session(overall_vars$session_token, overall_vars, "report")
         
         save_session(overall_vars$session_token, overall_vars, "gsea")
         
       })
     }) # close progress bar
     
   })
   
   observeEvent(input$table_gsea_rows_selected, {
     
     output$pathway_gsea <- renderPlot({
       
       table <- overall_vars$gsea[[input$select_gsea]][[1]]
       
       pathway  <-  table[input$table_gsea_rows_selected,]$pathway
       print(pathway)
       
       p <- plot_pathway(selected_pathway = pathway, 
                         nrank = overall_vars$gsea[[input$select_gsea]][[2]], 
                         species = overall_vars$report$gsea[[input$select_gsea]]$organism)
       p + theme(text = element_text(size = 18))
     })
     
     tab <- overall_vars$gsea[[input$select_gsea]][[1]]
     leadingEdge <- tab$leadingEdge[[input$table_gsea_rows_selected]]
     output$leading_edge <- renderText({ c("Leading edge: ",leadingEdge) })
   })
   
   observeEvent(input$go_plot_gsea,{
     
     output$summary_gsea <- renderPlot({
       
       p <- plot_gsea(tab = overall_vars$gsea[[input$select_gsea]][[1]], 
                      pval_cut = as.numeric(input$gsea_pval), 
                      NES_cut = as.numeric(input$gsea_NES),
                      ID = input$select_gsea)
       p
     }, height = 700, width = 1200)
     
output$table_gsea <- DT::renderDataTable({ 
  table <- as.data.frame(overall_vars$gsea[[input$select_gsea]][[1]][,c(1,2,3,4,5,6,7)])
  table <- table[order(table[["NES"]], decreasing = TRUE),]
       
  colnames(table) <- c("Pathway",
                       "P-value", 
                        "Adjusted P-value", 
                        "Log2err", 
                        "Enrichment Score", 
                        "Normalized Enrichment Score",
                        "Size")
       DT::datatable(
         { table },
         
         extensions = 'Buttons',
         rownames= FALSE,
         
         options = list(
           paging = TRUE,
           searching = TRUE,
           fixedColumns = FALSE,
           autoWidth = FALSE,
           ordering = TRUE,
           dom = 'Bfrtip',
           lengthMenu=c(10,20,50),
           buttons = c('copy', 'csv', 'excel')),
         
         class = "display")}, 
       selection = 'single', server = FALSE)
   }) #close observeEvent
   

}#close server
  
  
 
 shinyApp(ui = ui, server = server)


