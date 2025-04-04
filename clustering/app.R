#INFORMATION-----------------------------

#LOAD LIBRARIES ------------------------
library(clustree)
library(data.table)
library(dplyr)
library(DT)
library(ff)
library(ggfun)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(gridExtra)
library(Matrix)
library(matrixStats)
library(plotly)
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

source("tab_SEURAT/UI.R")
source("tab_SEURAT/server.R")

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
  navbarPage("scStudio: Clustering",
          tabPanel("Seurat",
                     tab_SEURAT)
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$cgraph_plot <- renderDataTable(data.frame())
  output$cgraph_clustree <- renderDataTable(data.frame())
  output$mkgs_cgraph <- renderDataTable(data.frame())
  output$violins_cgraph <- renderDataTable(data.frame())



  }, silent = TRUE)
  
# OVERALL SESSION VARIABLES -----------------------------------
  overall_vars <- reactiveValues(metadata = data.frame())
  
# timer  ---------------------------------------------
autoInvalidate <- reactiveTimer(30000)  

## page configuration -----------------------------------  
output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)})

# UPLOAD SESSION ------------------------------------
## main modal -----------------------------------  
   uploadModal <- function(failed = FALSE) {
     modalDialog(
       h1("Welcome to scStudio - Clustering"),
       p("In this tab, you will be able to perform clustering analysis of your dataset.
         Please load a session to begin."),
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
         
      incProgress(0.1, detail = "Retrieving count matrix...")
         
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
       
       overall_vars$md5$countMatrices <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
             
             incProgress(0.5, detail = "Retrieving metadata...")
             
             overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
             overall_vars$md5$metadata <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
             
             overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
             
             incProgress(0.6, detail = "Finishing...")
             
             if ("doublets.rds" %in% files){
               overall_vars$doublets <- upload_session(overall_vars$session_token, "doublets")
             } #close if 
             
             else{overall_vars$doublets <- list()}
             
             
             if ("colours.rds" %in% files){
               overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
             } #close if 
             
             else{overall_vars$colours <- list()}
             
             if ("report.rds" %in% files){
               overall_vars$report<- upload_session(overall_vars$session_token, "report")
             } #close if
             else {overall_vars$report <- list()}
  
             
             if ("hvgs.rds" %in% files){
               overall_vars$hvgs <- upload_session(overall_vars$session_token, "hvgs")
             } #close if 
             else {overall_vars$hvgs <- list()}
             
             if ("dimred.rds" %in% files){
               overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
             } #close if 
             else {
               overall_vars$dimred <- list()
               save_session(overall_vars$session_token ,overall_vars, "dimred")
             }
             
             overall_vars$md5$dimred <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
             
             updateSelectInput(session, "select_dimred_cgraph",
                               choices = c(names(overall_vars$dimred$tsne),
                                           names(overall_vars$dimred$pca),
                                           names(overall_vars$dimred$umap)))
             
             if ("clusters.rds" %in% files){
               overall_vars$clusters <- upload_session(overall_vars$session_token, "clusters")
             } #close if 
             else {
               overall_vars$clusters <- list()
               save_session(overall_vars$session_token ,overall_vars, "clusters")
             }
             overall_vars$md5$clusters <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
             
             if ("mkgs.rds" %in% files){
               overall_vars$mkgs <- upload_session(overall_vars$session_token, "mkgs")
             } #close if 
             else {
               overall_vars$mkgs <- list()
               save_session(overall_vars$session_token ,overall_vars, "mkgs")
             }
             overall_vars$md5$mkgs <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/mkgs.rds"))
             
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
             
            incProgress(0.1, detail = "Retrieving count matrix...")
          
           overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
           
           overall_vars$md5$countMatrices <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
          
           incProgress(0.5, detail = "Retrieving metadata...")
           
           overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
           
           overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
           
           overall_vars$md5$metadata <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
           
           overall_vars$md5$mat_names <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
           
           incProgress(0.6, detail = "Finishing...")
          
           
           if ("doublets.rds" %in% files){
             overall_vars$doublets <- upload_session(overall_vars$session_token, "doublets")
           } #close if 
           
           else{overall_vars$doublets <- list()}
           
           if ("colours.rds" %in% files){
             overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
           } #close if 
          
            else{overall_vars$colours <- list()}
           
           if ("report.rds" %in% files){
             overall_vars$report <- upload_session(overall_vars$session_token, "report")
           } #close if 
          
            else{overall_vars$report <- list()}
        
           if ("hvgs.rds" %in% files){
             overall_vars$hvgs <- upload_session(overall_vars$session_token, "hvgs")
           } #close if 
           else {overall_vars$hvgs <- list()
                 save_session(overall_vars$session_token ,overall_vars, "hvgs")}
           
           overall_vars$md5$hvgs <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/hvgs.rds"))
           
           if ("dimred.rds" %in% files){
             overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
           } #close if 
           else {
             overall_vars$dimred <- list()
             save_session(overall_vars$session_token ,overall_vars, "dimred")
           }
           
           overall_vars$md5$dimred <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
           
           updateSelectInput(session, "select_dimred_cgraph",
                             choices = c(names(overall_vars$dimred$tsne),
                                         names(overall_vars$dimred$pca),
                                         names(overall_vars$dimred$umap)))
           
           if ("clusters.rds" %in% files){
             overall_vars$clusters <- upload_session(overall_vars$session_token, "clusters")
           } #close if 
           else {
             overall_vars$clusters <- list()
             save_session(overall_vars$session_token ,overall_vars, "clusters")
           }
           overall_vars$md5$clusters <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
           
           if ("mkgs.rds" %in% files){
             overall_vars$mkgs <- upload_session(overall_vars$session_token, "mkgs")
           } #close if 
           else {
             overall_vars$mkgs <- list()
             save_session(overall_vars$session_token ,overall_vars, "mkgs")
           }
           overall_vars$md5$mkgs <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/mkgs.rds"))
           
           
           
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
      
      check_md5sum_countMatrices <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
      } # close if
      
      if (check_md5sum_countMatrices != overall_vars$md5$countMatrices){
        overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
      } # close if
      
    }) #close try
    
    
  md5_dimred <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
  
    try({
      
      if (md5_dimred != overall_vars$md5$dimred){
        
        overall_vars$dimred <- readRDS(
          paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
        overall_vars$md5$dimred <- md5_dimred
        
        updateSelectInput(session, "select_pcs_cgraph",
                          choices = names(overall_vars$dimred$pca))
        
        updateSelectInput(session, "select_dimred_cgraph",
                          choices = 
                            c(names(overall_vars$dimred$tsne),
                              names(overall_vars$dimred$pca),
                              names(overall_vars$dimred$umap)))
        
        showNotification("New dimensionality reduction plot available.")
      }
    }, silent = TRUE)
  
  md5_hvgs <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/hvgs.rds"))
  
  try({
    
    if (md5_hvgs != overall_vars$md5$hvgs){
      
      overall_vars$hvgs <- readRDS(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/hvgs.rds"))
      overall_vars$md5$hvgs <- md5_hvgs
      
      updateSelectInput(session, "select_features_cgraph",
                        choices = c("All",names(overall_vars$hvgs)))
      
    }
    
    md5_mkgs <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/mkgs.rds"))

      
      if (md5_mkgs != overall_vars$md5$mkgs){
        
        overall_vars$mkgs <- readRDS(
          paste0(getwd(),"/tokens/", overall_vars$session_token, "/mkgs.rds"))
        overall_vars$md5$mkgs <- md5_mkgs
        
        options <- grep(input$select_clustering_analysis, names(overall_vars$mkgs), value = TRUE)
        options <- gsub(input$select_clustering_analysis, "", options)
        options <- gsub("[0-9]\\.[0-9]", "", options)
        options <- gsub("[0-9]", "", options)
        options <- unique(gsub("_", "", options))
        
        updateSelectInput(session, "select_dea_method_cgraph",
                          choices = options)
        
        showNotification("New marker gene analysis available.")
        
        
      } #close if 
      
    
  }, silent = FALSE)
  
## clustering update ---------------------------------------           
  try({
    
    md5_clusters <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
    
    if (md5_clusters != overall_vars$md5$clusters){
      
      overall_vars$clusters <- readRDS(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
      
      overall_vars$md5$clusters <- md5_clusters
      
    }
  }, silent = FALSE)
  

  }) #close observe
   
## update inputs ----------------------------------------------------   
observe({
  updateSelectInput(session, 
                    "select_clustering_analysis", 
                    choices = names(overall_vars$clusters))
                    
  if (input$genesvspcs_cgraph == "Genes"){
    
    output$select_features_cgraph <- renderUI(
      selectInput(inputId = "select_features_cgraph", 
                  label = "Select features", 
                  choices = c("All",names(overall_vars$hvgs)), 
                  multiple = FALSE, 
                  width = "100%"))
    
    output$select_matrix_seurat <- renderUI(  
      selectInput(inputId = "select_matrix_seurat", 
                  label = "Select count matrix", 
                  choices = overall_vars$mat_names, 
                  multiple = FALSE, width = "100%"))
    
    output$select_pcs_cgraph <- NULL
    
    output$ncomponents_cgraph <- NULL
    }
  
  else {
    output$select_pcs_cgraph <- renderUI(
      selectInput(inputId = "select_pcs_cgraph", 
                  label = "Select PCA", 
                  choices = names(overall_vars$dimred$pca), 
                  multiple = FALSE, 
                  width = "100%"))
    
    output$ncomponents_cgraph <- renderUI(
      textInput(inputId = "ncomponents_cgraph", 
                label = "Number of PCs:",value = 25))
    
    output$select_features_cgraph <- NULL
    
    output$select_matrix_seurat <- NULL
  }
  
  updateSelectInput(session, "select_matrix_run_mkgenes",
                    choices = overall_vars$mat_names[-1])
  
  options <- grep(input$select_clustering_analysis, names(overall_vars$mkgs), value = TRUE)
  options <- gsub(input$select_clustering_analysis, "", options)
  options <- gsub("[0-9]\\.[0-9]", "", options)
  options <- gsub("[0-9]", "", options)
  options <- unique(gsub("_", "", options))
  
  updateSelectInput(session, "select_dea_method_cgraph",
                    choices = options)
  
  if(input$select_dimred_cgraph %in% names(overall_vars$dimred$pca)){
    coordinates <- overall_vars$dimred$pca[[input$select_dimred_cgraph]][["x"]]
  }
  else if (input$select_dimred_cgraph %in% names(overall_vars$dimred$tsne)){
    coordinates <- overall_vars$dimred$tsne[[input$select_dimred_cgraph]]
  }
  else {
    coordinates <- overall_vars$dimred$umap[[input$select_dimred_cgraph]]
  }
  
  output$cgraph_plot <- renderPlotly({
    
    req(input$select_clustering_analysis)
    
    p <- plot_dimRed_clusters(coord = coordinates, 
                              x = 1, 
                              y = 2, 
                              metadata = overall_vars$clusters[[input$select_clustering_analysis]], 
                              sel_var = paste0("res_",input$resolution_cgraph), 
                              title = "",
                              cols = extra_cols)
    toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
              layout(dragmode = "lasso", height = 550, 
                     legend = list(orientation = "h", x = 0, y = -0.2))) 
  }) #close cgraph_plot
  
  output$mkgs_cgraph <- DT::renderDataTable({ 
    
    print("Plotting marker genes table")
    
    name <- paste(input$select_clustering_analysis, 
                  input$select_dea_method_cgraph, 
                  input$resolution_cgraph, 
                  sep ="_")

    table <- overall_vars$mkgs[[name]]
    
    filter <- table[,"Adjusted p-value"] <= as.numeric(input$adjP_cgraph)
    
    if (input$selected_cluster_cgraph ==""){
      
      print("No specific cluster selected.")
      
      DT::datatable(
        { table[filter,]},
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
          buttons = c('copy', 'csv', 'excel')
        ),
        class = "display"
      )
    }
    else {
      DT::datatable(
        { table[table$Cluster == input$selected_cluster_cgraph & filter, ] },
        
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
          buttons = c('copy', 'csv', 'excel')
        ),
        class = "display"
      )
    }
  }, 
  selection = 'single', server = FALSE)
  
  output$cgraph_clustree <- renderPlot({
    req(input$select_clustering_analysis)
    p <- plot_clustree(overall_vars$clusters[[input$select_clustering_analysis]], 
                       input$resolution_cgraph)
    p
    
    #toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
    #          layout(dragmode = "lasso", height = 800, 
    #                legend = list(orientation = "h", x = 0, y = -0.2))) 
  })  
  
  #updateSelectInput(session, "select_dimred_cgraph",
  #                  choices = c(names(overall_vars$dimred$tsne),
  #                              names(overall_vars$dimred$pca),
  #                              names(overall_vars$dimred$umap)))

}) #close observe
 
   
# RUN SEURAT -----------------------------------------------  

   observeEvent(input$run_cgraph, {
     
     withProgress(message = 'Calculating clusters...', value = 0, {
       incProgress(0.4, detail = "")
       
       subset_row <- c()  
       
       try({ if (input$select_features_cgraph == "All"){
         
         subset_row <- rownames(overall_vars$countMatrices[[input$select_matrix_seurat]])} 
         
         else {subset_row <- overall_vars$hvgs[[input$select_features_cgraph]]}
       }) #close try
       
       print("Calculating clusters using Seurat's Graph-based approach.")
       
       if(input$genesvspcs_cgraph == "PCs"){use_mat <- overall_vars$countMatrices[[1]]}
       else{use_mat <- overall_vars$countMatrices[[input$select_matrix_seurat]]}
       
       clusters <- run_graph(mat = use_mat, 
                             pca = overall_vars$dimred$pca[[input$select_pcs_cgraph]][["x"]], 
                             genesvspcs = input$genesvspcs_cgraph, 
                             features = subset_row, 
                             dims = c(1:as.numeric(input$ncomponents_cgraph)), 
                             token = overall_vars$session_token, 
                             session_obj =  overall_vars$clusters, 
                             ID = input$cgraph_id)
       
       overall_vars$clusters[[input$cgraph_id]] <- clusters 
       
       print("Clustering job finished successfully.")
       
       save_session(overall_vars$session_token, overall_vars, "clusters")
   
   
   
       if(input$genesvspcs_cgraph == "PCs"){
      
         overall_vars$report$clustering[[input$cgraph_id]]$method <- "PCA"

         overall_vars$report$clustering[[input$cgraph_id]]$matrix <- input$select_pcs_cgraph

         overall_vars$report$clustering[[input$cgraph_id]]$features <- input$ncomponents_cgraph

 
       }
       else {
         overall_vars$report$clustering[[input$cgraph_id]]$method <- "Genes"
         overall_vars$report$clustering[[input$cgraph_id]]$matrix <- input$select_matrix_seurat
         overall_vars$report$clustering[[input$cgraph_id]]$features <- input$select_features_cgraph
         
       }

       save_session(overall_vars$session_token, overall_vars, "report")


       print("Report on clustering assay was updated.")
     }) # close progress bar
    })
   
   #**************************************************************************************    
   
   observeEvent(input$mkgs_cgraph_rows_selected, 
                
                output$violins_cgraph <- renderPlotly({
                  
                name <- paste(input$select_clustering_analysis, 
                              input$select_dea_method_cgraph, 
                              input$resolution_cgraph, 
                              sep ="_")

                table <- overall_vars$mkgs[[name]]
                  
                filter <- table[,"Adjusted p-value"] <= as.numeric(input$adjP_cgraph)
                  
       
                if (input$selected_cluster_cgraph == ""){
            
                    
                    table <- table[filter,]
                }
                
                else {
                    table <- table[table$Cluster == input$selected_cluster_cgraph & filter, ]
                    
                }
                  
                gene  <-  table[input$mkgs_cgraph_rows_selected,]$Gene
                  
                clusters <- overall_vars$clusters[[input$select_clustering_analysis]]
                clusters <- clusters[[paste0("res_", input$resolution_cgraph)]]
                  
                mat <- overall_vars$report$mkgs[[name]]$matrix
                  
                  
                g <- plot_violin_clusters(mat = overall_vars$countMatrices[[mat]], 
                                          clusters = clusters, 
                                          gene = gene) 
                  
                  toWebGL(ggplotly(g, source = "main", key = "pointNumber", tooltip = "text") %>% layout(dragmode = "lasso", 
                                                                                                         height = 400,
                                                                                                         width = 500)) 
                  })) #close plot_violin_clusters
   
   observeEvent(input$run_mkgenes,{ 
     
     withProgress(message = 'Preparing data for marker gene analysis...', value = 0, {
       incProgress(0.4, detail = "")
       
       mat <- overall_vars$countMatrices[[input$select_matrix_run_mkgenes]]
       norm_mat <- (2**mat)-1
       ln_mat <- log(norm_mat + 1)
       
       incProgress(0.5, detail = "")
       
       srt <- Seurat::CreateSeuratObject(counts = ln_mat, min.cells = 0, min.features = 0)
     })
     
     print("Sending marker genes job")
     print("-----------------------------------------------------------------")
     
     showNotification("Running analysis in background.")
     overall_vars$jobs[[paste(input$select_clustering_analysis, 
                        input$select_method_dea_cgraph[1] ,
                        input$resolution_cgraph, 
                        sep = "_")]]  <- callr::r_bg(run_marker_genes, 
                                                     args = list(srt = srt, 
                                                                 mat = overall_vars$countMatrices[[input$select_matrix_run_mkgenes]],
                                                                 methods = input$select_method_dea_cgraph, 
                                                                 ths = input$minLogFC_cgraph, 
                                                                 all_clusters = overall_vars$clusters[[input$select_clustering_analysis]], 
                                                                 session_obj = overall_vars$mkgs, 
                                                                 token = overall_vars$session_token, 
                                                                 ID = input$select_clustering_analysis,
                                                                 sres = input$resolution_cgraph))

     for (method in input$select_method_dea_cgraph){
       name <- paste(input$select_clustering_analysis, method, input$resolution_cgraph, sep ="_")
       overall_vars$report$mkgs[[name]]$matrix <- input$select_matrix_run_mkgenes
       overall_vars$report$mkgs[[name]]$ths <- input$minLogFC_cgraph
     } #close for 
     
     save_session(overall_vars$session_token, overall_vars, "report")
     
   }) #close observeEvent run_mkgenes
   
   # Annotate cell selection
   
   observeEvent(event_data('plotly_selected', source = 'main'), {
     cells <- event_data('plotly_selected', source = 'main')
     
     showNotification(paste0("Selected ", nrow(cells) , " cells"), duration = 20)
     
     showModal(modalDialog(
       h1("Cell Selection"),
       p("Select how you wish to annotate your cells:"),
       
       h3("Option 1 - new variable"),
       textInput("new_variable_name", "Name:",
                 placeholder = 'E.g.: cell_type'),
       textInput("new_variable_group1", "Selection annotation",
                 placeholder = 'E.g.: Cell A'),
       textInput("new_variable_group2", "Other cells",
                 placeholder = 'E.g.: Cell B'),
       
       actionButton("go_annotate_new", "OK"),
       
       h3("Option 2 - existing variable"),
       
       selectInput("existing_variable_name", 
                   label = "", 
                   choices = colnames(overall_vars$metadata),
                   selected = NULL),
       
       textInput("replace_annotation", "Selection annotation",
                 placeholder = 'E.g.: Cell A'),
       
       radioButtons(inputId = "var_annotation_type",
                    label = "",
                    choices = c("Merge", "Replace"),
                    selected = "Replace",
                    inline = TRUE),
       
       actionButton("go_annotate_existing", "OK"),
       
       footer = tagList(
         
         modalButton("Cancel"))
     ) #close dataModal
     )
     
   }) # close selection  
   
   observeEvent(input$go_annotate_new, {
     removeModal()
     
     
     selected_cells <- as.numeric(event_data('plotly_selected', source = 'main')$key)
     
     overall_vars$metadata[[input$new_variable_name]] <- input$new_variable_group2
     
     overall_vars$metadata[[input$new_variable_name]][selected_cells] <- input$new_variable_group1
     
     save_session(overall_vars$session_token, overall_vars, "metadata") 
     
     showNotification("New annotation available.")
     
   }) #close observe Annotate new 
   
   
   observeEvent(input$go_annotate_existing, {
     removeModal()
     
     selected_cells <- as.numeric(event_data('plotly_selected', source = 'main')$key)
     
     if (input$var_annotation_type == "Merge"){
       
       previous <-  overall_vars$metadata[[input$existing_variable_name]][selected_cells]
       
       merged <- paste(previous, input$replace_annotation, sep ="-")
       
       overall_vars$metadata[[input$existing_variable_name]][selected_cells] <- merged 
     }
     
     else{
       
       overall_vars$metadata[[input$existing_variable_name]][selected_cells] <- input$replace_annotation
     }
     
     save_session(overall_vars$session_token, overall_vars, "metadata") 
     
     showNotification("New annotation available.")
     
   }) #close observe Annotate new    

}#close server
  
  
 
 shinyApp(ui = ui, server = server)


