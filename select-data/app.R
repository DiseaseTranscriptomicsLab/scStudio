# INFORMATION--------------------------------------------------

# WORK DIR ----------------------------------------------------

#LOAD LIBRARIES ---------------------------------------------
library(data.table)
library(dplyr)
library(DT)
library(GEOquery)
library(Matrix)
library(Seurat)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinythemes)
library(stringi)
library(stringr)
library(tools)
library(vroom)


print("Sucessfully loaded libraries.")

##-------------------------------------##

#LOAD TABS-------------------------------

source("setup.R")

source("tab_SELECTDATA/UI.R")
source("tab_SELECTDATA/server.R")

print("Successfully loaded tabs.")

##-------------------------------------##

#USER INFERFACE ----------------------------

ui <- fluidPage(
  ## theme ----------------------------------------
  
  useShinyjs(),
  theme = shinytheme("flatly"),
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  
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
  
  tags$script(src="https://cdn.plot.ly/plotly-latest.min.js"),
  

## navbarPage --------------------------------------- 
  navbarPage("scStudio",
           tabPanel("Upload data",
            tab_SELECTDATA       
           ) #close SELECTDATA
      )#close navbarPage 
  
) #close UI 

##-------------------------------------##

#SERVER ------------------------------------
server <- function(input, output, session) {

##anchors ---------------------------------
  try({
    
  output$count_matrix_preview <- renderDataTable(data.frame())
  
  output$metadata_preview <- renderDataTable(data.frame())
  
  }, silent = TRUE)
  


## OVERALL SESSION VARIABLES -----------------------------------
overall_vars <- reactiveValues(tempCountMatrices = list(),
                               metadata = data.frame(), 
                               countMatrices = list(),
                               mat_names= c(),
                               session_token = get_random_token(),
                               md5 = list()) #initiate when uploading or adding samples
  
##md5sum check ---------------------------------------------
  autoInvalidate <- reactiveTimer(30000)
  
  
##page configuration -----------------------------------  
output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)})  
  

#UPLOAD SESSION ------------------------------------
##main modal -----------------------------------  
   uploadModal <- function(failed = FALSE) {
     modalDialog(
       h4("Option 1"),
       textInput("upload_token", "Insert token",
                 placeholder = 'E.g.: exzy8s'),
       h4("Option 2"),
       fileInput(inputId = 'upload_session_zip', label = 'Upload token.zip', multiple = FALSE),
       span(""),
       if (failed)
         div(tags$b("Invalid token or file", style = "color: red;")),
       
       footer = tagList(
         actionButton("ok_upload_session", "Upload session"),
         modalButton("Cancel")
         )
       )} #close dataModal
   
   observeEvent(input$upload_session,{
     showModal(uploadModal())
   }) # close upload session

##ok_upload_session --------------------------------   
   observeEvent(input$ok_upload_session, {
     
###zip modality -------------------------------------    
     if (!is.null(input$upload_session_zip)) {
       removeModal()
       print(input$upload_session_zip)
       
       token <- gsub("\\.zip","",input$upload_session_zip$name)
       unzip(zipfile = input$upload_session_zip$datapath,
             exdir = paste0("./tokens/", "temp"))
       old <- paste0("./tokens/temp","/tokens/",token)
       print(old)
       new <- paste0("./tokens/",token)
       print(new)
       file.move(old, new)
       unlink("./tokens/temp", recursive = TRUE)
       
       overall_vars$session_token <- token
       
       withProgress(message = 'Uploading data', value = 0, {
         incProgress(0.1, detail = "Retrieving count matrices...")
         overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
         
         overall_vars$md5$countMatrices <-  md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
         
         incProgress(0.5, detail = "Retrieving metadata...")
         overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
         
         overall_vars$md5$metadata <-  md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
         
         incProgress(0.6, detail = "Finishing...")
  
         incProgress(1, detail = "Done.")
       }) #close progress
     }
     
##token modality -------------------------------------      
     # Check if token exists
     else if (!is.null(input$upload_token)) {
       print(input$upload_token)
       removeModal()
       print(list.dirs(path = "./tokens/", full.names = FALSE)[-1])
       
       upload_token_in_dir <- input$upload_token %in% list.dirs(path = "./tokens/", full.names = FALSE)[-1]
       
       print(upload_token_in_dir)
       
       if (!upload_token_in_dir ) {showModal(uploadModal(failed = TRUE))} #close second if
       else{
         removeModal()
         overall_vars$session_token <- input$upload_token
         
         withProgress(message = 'Uploading data', value = 0, {
         incProgress(0.1, detail = "Retrieving count matrices...")
         overall_vars$countMatrices <- upload_session(input$upload_token, "countMatrices")
         
         overall_vars$md5$countMatrices <-  md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
         
         
         incProgress(0.5, detail = "Retrieving metadata...")
         overall_vars$metadata <- upload_session(input$upload_token, "metadata")
         
         overall_vars$md5$metadata <-  md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
         
         incProgress(0.6, detail = "Finishing...")
        
         
         
         incProgress(1, detail = "Done.")
         
         }) #close progress bar
       } #close else
     } #close first if 
    }) #close upload session
   
   
#SAVE SESSION -------------------------------------------
##main modal --------------------------------------------   
     saveModal <- function(failed = FALSE) {
       modalDialog(
         
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
   
##ok_download ------------------------------------------
   output$ok_download <- downloadHandler(
     filename = function() {
       paste0(overall_vars$session_token, ".zip")
     },
     content = function(file) {
       withProgress(message = 'Downloading latest saved session...', value = 0, {
         incProgress(0.4, detail = "")
       on.exit(removeModal())
       
       saveRDS(rownames(overall_vars$countMatrices$rawCountMatrix),
               paste0(getwd(), "/tokens/", overall_vars$session_token, "/genes.rds"))
       
       # Zip the selected files
       zip(file, paste0("tokens/", overall_vars$session_token))
       }) # close progress bar
       }
     )

##ok_create_session -------------------------------------      
   observeEvent(input$ok_create_session, {
     
     new_token <- get_random_token()
     
     removeModal()
     
     showModal(modalDialog(
       title = "New token generated", #if token is null, we create a new one
       paste("Use the following token to retrieve your session:",new_token)
     )) #close showModal
     
     withProgress(message = 'Saving data', value = 0, {
     incProgress(0.1, detail = "Saving count matrices...")
       
     save_session(new_token,overall_vars, c("countMatrices","metadata", "mat_names"))
     saveRDS(rownames(overall_vars$countMatrices$rawCountMatrix),
             paste0(getwd(), "/tokens/", overall_vars$session_token, "/genes.rds"))
     
     
     incProgress(0.9, detail = "Finishing...")
    
     }) #close progress bar
   }) # close create new session

##ok_current_session --------------------------------------   
   observeEvent(input$ok_current_session, {
    
     removeModal()
     
     withProgress(message = 'Saving data', value = 0, {
       
     token <- overall_vars$session_token
     
       incProgress(0.1, detail = "Saving count matrices...")
       
       save_session(token, overall_vars, c("countMatrices","metadata", "mat_names"))
       saveRDS(rownames(overall_vars$countMatrices$rawCountMatrix),
               paste0(getwd(), "/tokens/", overall_vars$session_token, "/genes.rds"))
       
       incProgress(0.9, detail = "Finishing...")
       
     }) #close progress bar
   }) # close create new session
  
   
# UPDATE SELECT OPTIONS/PLOTS ---------------------------------------------------
  observe({
 
## md5sum checks --------------------------------------------------------------
    try({
    # check if countMatrices file changed
    check_md5sum_countMatrices <- md5sum(
      paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
    
    print("Here")
    
    print(check_md5sum_countMatrices)
    
    print(overall_vars$md5$countMatrices)
    
    
    # update the list of countMatrices if new data is available from other scStudio tools
    if (check_md5sum_countMatrices != overall_vars$md5$countMatrices){
      
      
      overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
      overall_vars$mat_names <- names(overall_vars$countMatrices)
    } # close if
    
    # check if metadata file changed
    check_md5sum_metadata <- md5sum(
      paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
    
    # update metadata if new data is available from other scStudio tools
    
    if (check_md5sum_metadata != overall_vars$md5$metadata){
      
      overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
    } # close if
  }) #close try
    
## input options -------------------------------------------------  
    
    if (input$select_file_type == "table (.csv/.tsv/.txt)"){
      
      output$seurat_upload <- NULL
      output$sce_upload <- NULL
      output$matrix_10X <- NULL
      output$barcodes_10X <- NULL
      output$features_10X <- NULL
      output$h5_10X <- NULL
      
      output$matrix_upload <- renderUI(fileInput(inputId = 'matrix_upload', label = 'Gene expression matrix', multiple = FALSE))
      
      output$metadata_upload <- renderUI(fileInput(inputId = 'metadata_upload', label = 'Metadata', multiple = FALSE))
      
    } #close if 
    
    else if (input$select_file_type == "Seurat (.rds)" ){
      
      output$seurat_upload <- renderUI(fileInput(inputId = 'seurat_upload', label = 'Upload seurat.rds', multiple = FALSE))
      
      output$sce_upload <- NULL
      output$matrix_upload <- NULL
      output$metadata_upload <- NULL
      output$example_matrix <- NULL
      output$example_metadata <- NULL
      output$matrix_10X <- NULL
      output$barcodes_10X <- NULL
      output$features_10X <- NULL
      output$h5_10X <- NULL
      
    } #close else if
    
    else if (input$select_file_type == "SingleCellExperiment (.rds)" ){
      
      output$sce_upload <- renderUI(fileInput(inputId = 'sce_upload', label = 'Upload sce.rds', multiple = FALSE))
      
      output$seurat_upload <- NULL
      output$matrix_upload <- NULL
      output$metadata_upload <- NULL
      output$example_matrix <- NULL
      output$example_metadata <- NULL
      output$matrix_10X <- NULL
      output$barcodes_10X <- NULL
      output$features_10X <- NULL
      output$h5_10X <- NULL
      
    } #close else if
    
    else if (input$select_file_type == "10X Genomics (triplet format)" ){
      
      output$matrix_10X <- renderUI(
        fileInput(inputId = 'matrix_10X', label = 'Upload 10X matrix.mtx.gz', 
                  multiple = FALSE))
      
      output$barcodes_10X <- renderUI(
        fileInput(inputId = 'barcodes_10X', label = 'Upload 10X barcodes.tsv.gz', 
                  multiple = FALSE))
      
      output$features_10X <- renderUI(
        fileInput(inputId = 'features_10X', label = 'Upload 10X features.tsv.gz', 
                  multiple = FALSE))
      
      output$metadata_upload <- renderUI(fileInput(inputId = 'metadata_upload', label = 'Metadata', multiple = FALSE))
      
      output$seurat_upload <- NULL
      output$matrix_upload <- NULL
      output$example_matrix <- NULL
      output$example_metadata <- NULL
      output$sce_upload<- NULL
      output$h5_10X <- NULL
      
    } #close else if
    
    else if (input$select_file_type == "10X Genomics (.h5)" ){
      
      output$h5_10X <- renderUI(
        fileInput(inputId = 'h5_10X', label = 'Upload 10X matrix.h5', 
                  multiple = FALSE))
      
      output$metadata_upload <- renderUI(fileInput(inputId = 'metadata_upload', label = 'Metadata', multiple = FALSE))
      
      output$seurat_upload <- NULL
      output$matrix_upload <- NULL
      output$example_matrix <- NULL
      output$example_metadata <- NULL
      output$sce_upload<- NULL
      output$matrix_10X <- NULL
      output$barcodes_10X <- NULL
      output$features_10X <- NULL
      
    } # close else if 
   
   updateSelectInput(session, "availableCountMatrices",choices = names(overall_vars$countMatrices),
                     selected =  names(overall_vars$countMatrices))   
    
   updateSelectInput(session, "selectCountMatrices",choices = names(overall_vars$tempCountMatrices), 
                       selected = c(names(overall_vars$tempCountMatrices)[1]))   
  
   output$countMatricesProperties = renderDataTable(describe_matrices(overall_vars$countMatrices))
  
   output$sampleProperties = renderDataTable(describe_matrices(overall_vars$tempCountMatrices))
   
   output$metadata_preview = renderDataTable({datatable(overall_vars$metadata, 
                                                        options = list(
                                                          searching = FALSE,
                                                          pageLength = 5,
                                                          lengthMenu = c(5, 10, 20),
                                                          scrollX = T))})
   
   updateSelectInput(session, "subsetMetadata",choices = colnames(overall_vars$metadata), 
                     selected = c(colnames(overall_vars$metadata)))
   }) #close observe
   

# RUN UPLOAD pre-processed public data --------------------------------------------------   
   
    observeEvent(input$view_preprocessed,{ 
      
      existing_tokens <- list.files(paste0(getwd(),"/tokens/"))
      
      if (!(overall_vars$session_token %in% existing_tokens)){
        dir.create(paste0(getwd(),"/tokens/", overall_vars$session_token))
        
        withProgress(message = 'Uploading dataset', value = 0, {
          
          incProgress(0.4, detail = "...")
          
          overall_vars$countMatrices <- upload_preprocessed(input$selectData, "countMatrices")
          
          names(overall_vars$countMatrices) <- names(overall_vars$countMatrices)
          
          print("countMatrices uploaded")
          
          incProgress(0.4, detail = "...")
          
          overall_vars$metadata <- upload_preprocessed(input$selectData, "metadata")
          
          print("metadata uploaded")
          
          incProgress(0.1, detail = "Preparing session...")
          
          files <- list.files(paste0(getwd(), "/public_datasets/", input$selectData))
          for (file in files){
            file.copy(paste0(getwd(), "/public_datasets/", input$selectData,"/",file), paste0(getwd(),"/tokens/", overall_vars$session_token),
                      overwrite = TRUE)}
        }) #close incProgress
        
        try({
          overall_vars$md5$countMatrices <-  md5sum(
            paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
          
          overall_vars$md5$metadata <-  md5sum(
            paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
        }) #close try
  
      }
      else {
        showModal(modalDialog(
          title = "Would you like to overwrite the current session?",
          span(""),
          footer = tagList(
           
            actionButton("ok_overwrite_withPublic", "Overwrite"),
            actionButton("ok_newSession_withPublic", "Create new session"),
            modalButton("Cancel")
          ))) 
      }
    }) #close observeEvent
   
   observeEvent(input$ok_overwrite_withPublic ,{
     removeModal()
     withProgress(message = 'Uploading dataset', value = 0, {
       incProgress(0.4, detail = "...")
       overall_vars$countMatrices <- upload_preprocessed(input$selectData, "countMatrices")
       print("countMatrices uploaded")
       overall_vars$tempCountMatrices <- overall_vars$countMatrices
       names(overall_vars$countMatrices) <- names(overall_vars$countMatrices)
       
       incProgress(0.4, detail = "...")
       overall_vars$metadata <- upload_preprocessed(input$selectData, "metadata")
       print("metadata uploaded")
      
       incProgress(0.1, detail = "Preparing session...")
       files <- list.files(paste0(getwd(), "/public_datasets/", input$selectData))
       for (file in files){
       file.copy(paste0(getwd(), "/public_datasets/", input$selectData,"/",file), paste0(getwd(),"/tokens/", overall_vars$session_token),
                 overwrite = TRUE)}
     }) #close progress
     
   }) #close ObserveEvent
     
     observeEvent(input$ok_newSession_withPublic ,{
       removeModal()
       print("Creating new session for use with public dataset")
       
      overall_vars$session_token <- get_random_token() 
      dir.create(paste0(getwd(), "/tokens/", overall_vars$session_token))
       
      withProgress(message = 'Uploading dataset', value = 0, {
        incProgress(0.4, detail = "...")
        
      overall_vars$countMatrices <- upload_preprocessed(input$selectData, "countMatrices")
      
      overall_vars$mat_names <- names(overall_vars$countMatrices)
      
      print("countMatrices uploaded")
      
      overall_vars$tempCountMatrices <- overall_vars$countMatrices
      
      names(overall_vars$countMatrices) <- names(overall_vars$countMatrices)
      
      incProgress(0.4, detail = "...")
      
      overall_vars$metadata <- upload_preprocessed(input$selectData, "metadata")
      
      print("metadata uploaded")
    
      incProgress(0.1, detail = "Preparing session...")
      
      files <- list.files(paste0(getwd(), "/public_datasets/", input$selectData))
      
      for (file in files){
      
        file.copy(paste0(getwd(), "/public_datasets/", input$selectData,"/",file), paste0(getwd(),"/tokens/", overall_vars$session_token),
                overwrite = TRUE)}
      
      }) #close progress
    
   }) # close ObserveEvent


# RUN publish dataset -------------------------------------------------- 
   
   observeEvent(input$publish_analysis, {
     showModal(modalDialog(
       
       span(""),
       footer = tagList(
          textInput("public_dataset_name", "Name of the project:",
                    placeholder = "FernandesAM_2023"),
         actionButton("ok_publish", "Publish"),
        
         modalButton("Cancel"),
         p("This action will automatically save your session before publishing.")
       ))) 
    
   })
   
   observeEvent(input$ok_publish, {
     removeModal()
     names <- list.files(paste0(getwd(),"/public_datasets"))
     
     if (input$public_dataset_name %in% names){
       showModal(modalDialog(
         title = "Warning: name already in use."
       )) #close showModal
       
     }
     else{
     
     token_folder <- paste0("./tokens/", overall_vars$session_token)
     new_folder <- paste0("./public_datasets/", input$public_dataset_name)
     dir.create(new_folder)
     list_of_files <- list.files(token_folder, ".rds$") 
     # ".py$" is the type of file you want to copy. Remove if copying all types of files. 
     file.copy(file.path(token_folder,list_of_files), new_folder)
     # update the selection 
     updateSelectInput(session, "selectData",choices = get_dataset_names())
     }
   }) 
   
# RUN upload GEO -------------------------------------------------- 
    
   observeEvent(input$view_GEO,{ 
     
    withProgress(message = 'Downloading data', value = 0, {
      
      incProgress(0.5, detail = "Retrieving count matrices")
      
      try_geo <- try({
        
    temp <- get_geo(input$GEO_id)
    overall_vars$tempCountMatrices <- c(overall_vars$tempCountMatrices ,temp[[1]])
    overall_vars$metadata$percentage_mt_genes <- NULL
    overall_vars$metadata$library_size <- NULL
    overall_vars$metadata$total_features <- NULL
    overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp[[2]])
    
    },
    silent = FALSE
      ) #close try
    }) #close progress bar
     print(class(try_geo))
     if(class(try_geo) == "try-error"){showNotification("Download failed.")}
  }) #close view_geo 
   

    
# RUN user upload -------------------------------------------------- 
  
  observeEvent(input$view_upload,{

    # mandatory - count matrix needs to have 1st column for genes
    try({
    withProgress(message = 'Uploading data', value = 0, {
      
      if (input$upload_id %in% names(overall_vars$tempCountMatrices)){
      showModal(modalDialog(
        title = "Dataset name is duplicated.",
        p("Please choose a different name."),
        span(""),
        footer = tagList(
          modalButton("OK")
        ))) # close modal
      } # close if 
      
      else{
      if (input$select_file_type == "Seurat (.rds)"){
        
        temp <- load_srt(input$seurat_upload, input$upload_id)
        overall_vars$tempCountMatrices[[input$upload_id]] <- temp[[1]]
        overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp[[2]])
        rm(temp)
        
      } # close if
      
      else if (input$select_file_type == "SingleCellExperiment (.rds)"){
        temp <- load_sce(input$sce_upload, input$upload_id)
        overall_vars$tempCountMatrices[[input$upload_id]] <- temp[[1]]
        overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp[[2]])
        rm(temp)
        
      } # close else i
        
      else if (input$select_file_type == "10X Genomics (triplet format)"){
        temp <- load_10X_triplet(barcodes_path = input$barcodes_10X, 
                                 features_path = input$features_10X, 
                                 matrix_path = input$matrix_10X,
                                 name = input$upload_id) 
        
        overall_vars$tempCountMatrices[[input$upload_id]] <- temp[[1]]
        overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp[[2]])
        rm(temp)
        
        infile2 <- input$metadata_upload
        temp <- as.data.frame(vroom(infile2$datapath))
        
        print("Uploaded metadata:")
        orig.ident <- rep(input$upload_id, dim(temp)[1])
        temp$orig.ident <- orig.ident
        temp$dataset <- orig.ident
        
        cell_ids <- define_cell_ids(ncol(overall_vars$tempCountMatrices[[input$upload_id]]))
        #for (cell in 1:ncol(overall_vars$tempCountMatrices[[input$upload_id]])){
        #  cell_ids <- c(cell_ids, paste0(key, stri_rand_strings(1, 3, pattern = "[a-z0-9]")))
        #}
        temp$cell_id <- cell_ids
        
        overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp)
        
      } # close else if 
        
        else if (input$select_file_type == "10X Genomics (.h5)"){
          temp <- load_10X_h5(h5_path = input$h5_10X, 
                                   name = input$upload_id) 
          
          overall_vars$tempCountMatrices[[input$upload_id]] <- temp[[1]]
          overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp[[2]])
          rm(temp)
          
          infile2 <- input$metadata_upload
          temp <- as.data.frame(vroom(infile2$datapath))
          
          print("Uploaded metadata:")
          orig.ident <- rep(input$upload_id, dim(temp)[1])
          temp$orig.ident <- orig.ident
          temp$dataset <- orig.ident
          
          cell_ids <- define_cell_ids(ncol(overall_vars$tempCountMatrices[[input$upload_id]]))
          #for (cell in 1:ncol(overall_vars$tempCountMatrices[[input$upload_id]])){
          #  cell_ids <- c(cell_ids, paste0(key, stri_rand_strings(1, 3, pattern = "[a-z0-9]")))
          #}
          temp$cell_id <- cell_ids
          
          overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp)
          
        } # close else if     
        
      else if (input$select_file_type == "table (.csv/.tsv/.txt)"){
      infile1 <- input$matrix_upload
      infile2 <- input$metadata_upload
      
      incProgress(0.4, detail = "Uploading count matrix")
      
      overall_vars$tempCountMatrices[[input$upload_id]] <- upload_inhouse(infile1$datapath)
      
      print("Uploaded matrix:")
     
      incProgress(0.4, detail = "Uploading metadata")
      
      temp <- as.data.frame(vroom(infile2$datapath))
      
      print("Uploaded metadata:")
      orig.ident <- rep(input$upload_id, dim(temp)[1])
      temp$orig.ident <- orig.ident
      temp$dataset <- orig.ident
      
      cell_ids <- define_cell_ids(ncol(overall_vars$tempCountMatrices[[input$upload_id]]))
      #for (cell in 1:ncol(overall_vars$tempCountMatrices[[input$upload_id]])){
      #  cell_ids <- c(cell_ids, paste0(key, stri_rand_strings(1, 3, pattern = "[a-z0-9]")))
      #}
      temp$cell_id <- cell_ids
      
      overall_vars$metadata <- merge_metadata(overall_vars$metadata, temp)
      } # close else
      incProgress(0.4, detail = "Finishing")
      showNotification("Successfully uploaded data.")
      } #close else
      
    }) # close progress bar
    
    
  }) #close try 
    
  }) #close view_upload 
  
# RUN Add samples -----------------------------------------------------
  # when make.unique colnames, update barcodes column in metadata
  observeEvent(input$kickoff,{
  try({ 
  # if we already have count data in the analysis, verify if user wants to delete previous analyses  
    if ("rawCountMatrix" %in% names(overall_vars$countMatrices)){
      
      showModal(modalDialog(
        title = "Would you like to overwrite the current session?",
        p("Merging datasets will erase previous analyses."),
        span(""),
        footer = tagList(
          
          actionButton("ok_overwrite_addSamples", "Overwrite"),
          actionButton("ok_newSession_addSamples", "Create new session"),
          modalButton("Cancel")
        ))) 
    }
  ### session without previous analyses -------------------------------------    
    else{
      withProgress(message = 'Preparing the dataset for analysis', value = 0, {
        incProgress(0.2, detail = "Merging data...")  
        
        run_merging <- merge_matrices(overall_vars$tempCountMatrices[input$selectCountMatrices])
        raw_mat <- run_merging[[1]]
        orig_ids <- run_merging[[2]]
        
        print("New raw count matrix ready")
        print(dim(raw_mat))
        
      
        overall_vars$metadata <- overall_vars$metadata[(overall_vars$metadata$cell_id %in% orig_ids), #& 
                                                       #(overall_vars$metadata$orig.ident %in% input$selectCountMatrices),
                                                        input$subsetMetadata]
        
        incProgress(0.7, detail = "Saving session...")
        
        overall_vars$countMatrices <- list("rawCountMatrix" = raw_mat)
        
        overall_vars$mat_names <- names(overall_vars$countMatrices)
        
        
      
      # run save session
        save_session(overall_vars$session_token, overall_vars, c("countMatrices", "metadata", "mat_names"))
        
        overall_vars$tempCountMatrices <- list()
        
        updateSelectInput(session, "selectCountMatrices",choices = c("")) 
        
        overall_vars$md5$countMatrices <-  md5sum(
          paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
        
        overall_vars$md5$metadata <-  md5sum(
          paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
        
      })# close withProgress
      showNotification("Dataset ready for analysis.", duration = 10)
    } # close else
    
  }) #close try   
}) #close add samples
  #***************************************************************************************  
  #*
  ### session overwrite -------------------------------------  
   observeEvent(input$ok_overwrite_addSamples, {
     removeModal()
     try({
     withProgress(message = 'Preparing the dataset for analysis', value = 0, {
       incProgress(0.2, detail = "Merging data...")  
       
       cm_list <- c(list(overall_vars$countMatrices[["rawCountMatrix"]]),overall_vars$tempCountMatrices[input$selectCountMatrices])
      
       run_merging <- merge_matrices(cm_list)
       
       raw_mat <- run_merging[[1]]
       orig_ids <- run_merging[[2]]
       
       print("New raw count matrix ready")
       print(dim(raw_mat))
       
       overall_vars$metadata <- overall_vars$metadata[(overall_vars$metadata$cell_id %in% orig_ids), #& 
                                                        #(overall_vars$metadata$orig.ident %in% input$selectCountMatrices),
                                                      input$subsetMetadata]
       incProgress(0.7, detail = "Saving session...")  
       
       overall_vars$countMatrices <- list("rawCountMatrix" = raw_mat)
       
       overall_vars$mat_names <- names(overall_vars$countMatrices)
       
       # run save session
       save_session(overall_vars$session_token, overall_vars, c("countMatrices", "metadata", "mat_names"))
       remove_objects(overall_vars$session_token)
       
       overall_vars$md5$countMatrices <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
       
       overall_vars$md5$metadata <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       
       overall_vars$tempCountMatrices <- list()
       
       updateSelectInput(session, "selectCountMatrices",choices = c(""))
       
     })# close withProgress
     showNotification("Dataset ready for analysis.", duration = 10)
     }) #close try  
   }) #close ok_overwrite_addSamples
#***************************************************************************************
#*
### new session -------------------------------------   
   observeEvent(input$ok_newSession_addSamples, {
     removeModal()
     try({
     withProgress(message = 'Preparing the dataset for analysis', value = 0, {
       incProgress(0.2, detail = "Merging data...")  
       
       overall_vars$session_token <- get_random_token()
       
       run_merging <- merge_matrices(
         c(list(overall_vars$countMatrices[["rawCountMatrix"]]),overall_vars$tempCountMatrices[input$selectCountMatrices]))
       
       raw_mat <- run_merging[[1]]
       orig_ids <- run_merging[[2]]
       
       print("New raw count matrix ready")
       print(dim(raw_mat))
       
       overall_vars$metadata <- overall_vars$metadata[(overall_vars$metadata$cell_od %in% orig_ids), #& 
                                                        #(overall_vars$metadata$orig.ident %in% input$selectCountMatrices),
                                                      input$subsetMetadata]
       incProgress(0.7, detail = "Saving new session...")  
       
       overall_vars$countMatrices <- list("rawCountMatrix" = raw_mat)
       overall_vars$mat_names <- names(overall_vars$countMatrices)
       
       # run save session
       save_session(overall_vars$session_token, overall_vars, c("countMatrices", "metadata", "mat_names"))
       
       overall_vars$md5$countMatrices <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
       
       overall_vars$md5$metadata <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       
       overall_vars$tempCountMatrices <- list()
       
       updateSelectInput(session, "selectCountMatrices",choices = c(""))
       
     })# close withProgress
     showNotification("Dataset ready for analysis.", duration = 10)
     }) #close try 
   }) #close ok_overwrite_addSamples
   
   
  # RUN Delete matrix --------------------------------------------------------------
   observeEvent(input$delete_countMatrix,{
     if("rawCountMatrix" %in% input$availableCountMatrices){
       showModal(modalDialog(
         title = "Cannot delete raw counts. Please delete session instead.",
         span(""),
         footer = tagList(
           modalButton("OK")
         ))) 
     } #close if 
     
     else{
       showModal(modalDialog(
         title = paste0("Confirm deletion of the the following data: ", input$availableCountMatrices ),
         p("You will not be able to recover the data."),
         span(""),
         footer = tagList(
           actionButton("ok_delete_matrix", "Proceed"),
           modalButton("Cancel")
         ))) 
     } #close else
  }) # close delete matrix

### ok delete matrix -----------------------------------------------    
     observeEvent(input$ok_delete_matrix, {
       removeModal()
       for (mat in input$availableCountMatrices){
         overall_vars$countMatrices[[mat]] <- NULL
       } #close for
       overall_vars$mat_names <- names(overall_vars$countMatrices)
       
     }) #close ok_delete_matrix
     
     
# DELETE SESSION  
   observeEvent(input$delete_session,{
   showModal(modalDialog(
     title = paste0("Confirm deletion of the the following session: ", overall_vars$session_token ),
     p("You will not be able to recover the data."),
     span(""),
     footer = tagList(
       actionButton("ok_delete_session", "Proceed"),
       modalButton("Cancel")
     ))) 
   }) #close observe delete session   
   
## ok delete session ---------------------------------------------------------
   observeEvent(input$ok_delete_session, {
     removeModal()
     unlink(paste0(getwd(), "/tokens/", overall_vars$session_token), recursive = TRUE)
     overall_vars <- reactiveValues(tempCountMatrices = list(),
                                    metadata = data.frame(), 
                                    countMatrices = list(),
                                    mat_names= c(),
                                    session_token = get_random_token(),
                                    md5 = list())
     
     updateSelectInput(session, "availableCountMatrices",choices = "")   
     
     updateSelectInput(session, "selectCountMatrices",choices = "")   
     
     output$countMatricesProperties = renderDataTable(data.frame())
     
     output$sampleProperties = renderDataTable(data.frame())
     
     output$metadata_preview = renderDataTable(data.frame())
     
     updateSelectInput(session, "subsetMetadata",choices = "")
     
     output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)}) 
     
    showNotification("Session deleted sucessfully.") 
   })
  
# To keep the app from running without the user selecting the initial dataset 
# all functions are inside an observeEvent, except the Select Data tab functions   

}#close server
  
  
 #shiny::runApp(launch.browser = FALSE)
 shinyApp(ui = ui, server = server)


