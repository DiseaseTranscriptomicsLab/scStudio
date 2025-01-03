#INFORMATION-----------------------------

#LOAD LIBRARIES ------------------------
library(data.table)
library(dplyr)
library(DT)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(gplots)
library(limma)
library(plotly)
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

print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_VAREXPLAINED/UI.R")
source("tab_VAREXPLAINED/server.R")

source("tab_DEPTH/UI.R")
source("tab_DEPTH/server.R")

source("tab_BEC/UI.R")
source("tab_BEC/server.R")

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
  navbarPage("scStudio: Normalization",
          tabPanel("Variance Explained",
                     tab_VAREXPLAINED      
             ), 
           tabPanel("Depth Normalization",
                    tab_DEPTH     
             ),
          tabPanel("Batch Effect Correction",
                   tab_BEC   
          )
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$boxplot_raw_depth <- renderDataTable(data.frame())
  output$boxplot_norm_depth <- renderDataTable(data.frame())
  output$var_exp <- renderDataTable(data.frame())
  output$boxplot_be <- renderDataTable(data.frame())


  }, silent = TRUE)
  
# OVERALL SESSION VARIABLES -----------------------------------
  overall_vars <- reactiveValues(metadata = data.frame())
  
# timer  ---------------------------------------------
autoInvalidate <- reactiveTimer(30000)  

##page configuration -----------------------------------  
output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)})

# UPLOAD SESSION ------------------------------------
## main modal -----------------------------------  
   uploadModal <- function(failed = FALSE) {
     modalDialog(
       h1("Welcome to scStudio - Normalization"),
       p("In this tab, you will be able to perform normalization of your dataset.
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
             
             incProgress(0.5, detail = "Retrieving metadata...")
             overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
             
             overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
             
             overall_vars$md5$metadata <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
             
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
             
             updateSelectInput(session, "features_sample",choices = names(identify_discrete(overall_vars$metadata)))   
             
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
        
           updateSelectInput(session, "features_sample",choices = names(identify_discrete(overall_vars$metadata)),
                             selected = "orig.ident") 
           
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

    ## md5sum checks --------------------------------------------------------------
    
    try({
      # check if files changed
      check_md5sum_metadata <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
      
      check_md5sum_mat_names <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
      } # close if
      
      
    }) #close try

 
# update inputs ----------------------------------------------------   
    print(overall_vars$mat_names)
    updateSelectInput(session, 
                      "select_matrix_varexp",
                      choices = overall_vars$mat_names)
                      
    updateSelectInput(session, 
                      "variance_vars",
                      choices = c(colnames(overall_vars$metadata)), 
                      selected = c("library_size", "total_features"))
    
    updateSelectInput(session, 
                      "select_matrix_depth",
                      choices = overall_vars$mat_names[-1])
    
    updateSelectInput(session, 
                      "var_depth_boxplot",
                      choices = colnames(identify_discrete(overall_vars$metadata)))
  
    ##bec -----------------------------------------------   
    
    updateSelectInput(session, "select_matrix_bec",choices = overall_vars$mat_names[-1])
    
    updateSelectInput(session, "select_matrix_plot_bec",choices = overall_vars$mat_names[-1])
    
    updateSelectInput(session, "var_be_boxplot",
                      choices =  names(identify_discrete(overall_vars$metadata)), 
                      selected = c("orig.ident"))
    
    updateSelectInput(session, "batch_label1",
                      choices =  names(identify_discrete(overall_vars$metadata)), 
                      selected = c("orig.ident"))
    
    updateSelectInput(session, "batch_label2",
                      choices =  c(names(identify_discrete(overall_vars$metadata)), "NULL"),
                      selected = c("NULL"))
 
    updateSelectInput(session, "adjust_covariate",
                      choices =  c(names(identify_continuous(overall_vars$metadata)), "NULL"),
                      selected = c("NULL"))
    
    updateSelectInput(session, "keep_biological",
                      choices =  c(names(identify_discrete(overall_vars$metadata)), "NULL"),
                      selected = c("NULL"))
}) #close observe
 
   
# RUN VARIANCE EXPLAINED -----------------------------------------------  

 observeEvent(input$run_variance, {
   output$var_exp <- renderPlotly(
     { 
       input$run_variance
         withProgress(message = 'Calculating variance...', value = 0, {
           incProgress(0.4, detail = "")
          
            countMatrices <- upload_session(overall_vars$session_token, "countMatrices") 
            isolate(p <- plot_variance_explained(
                         countMatrices[[input$select_matrix_varexp]],
                         input$variance_vars,
                         overall_vars$metadata)
                    
                    
            ) # close isolate
            
            toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                      layout(dragmode = "lasso", height = 450))    
         }) #close progress 
       }) # close var_exp
   }) # close observe run_variance
   
# RUN DEPTH NORMALIZATION -----------------------------------------------
   
   observeEvent(input$run_normalization, {
     withProgress(message = 'Normalizing data...', value = 0, {
       incProgress(0.4, detail = "")
       
       # get count matrices
       
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
       
       overall_vars$countMatrices[[input$norm_mat_name]] <- choose_norm_method(
                                                                  input$norm_method, 
                                                                  as.matrix(overall_vars$countMatrices[["rawCountMatrix"]]))
       
       overall_vars$mat_names <- names(overall_vars$countMatrices)
       
       overall_vars$report$normalization[[input$norm_mat_name]]$method <- input$norm_method
       
       save_session(overall_vars$session_token, overall_vars, "countMatrices")
       
       save_session(overall_vars$session_token, overall_vars, "mat_names")
       overall_vars$md5$mat_names <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
       
       save_session(overall_vars$session_token, overall_vars, "report")
       overall_vars$md5$report <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/report.rds"))
       
       norm_mat <- (2**overall_vars$countMatrices[[input$norm_mat_name]]) - 1
       
       overall_vars$metadata[paste0("library_size_", input$norm_mat_name)] <- 
         colSums(overall_vars$countMatrices[[input$norm_mat_name]])
       
       save_session(overall_vars$session_token, overall_vars, "metadata")
       overall_vars$md5$metadata <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       
       overall_vars$countMatrices <- NULL
       
     }) # close progress bar
   })# close observe run_normalization   
   
## RUN depth normalization boxplots    
   
   observeEvent(input$generate_norm_boxplots, {

     cols <- c(as.vector(unlist(overall_vars$colours[[input$var_depth_boxplot]])), extra_cols)
     
     output$boxplot_raw_depth <- renderPlotly(
       
       { isolate({p <- plot_boxplot(overall_vars$metadata, 
                                    "library_size", 
                                    input$var_depth_boxplot, 
                                    "Before normalization",
                                    "log10(library size)",
                                    input$select_boxplot_range1,
                                    cols =  cols)
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 600, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))
       }) #close isolate
       }) #close boxplot_raw_libsize
     
     output$boxplot_norm_depth <- renderPlotly(
       
       { isolate({p <- plot_boxplot(overall_vars$metadata, 
                                    paste0("library_size_", input$select_matrix_depth), 
                                    input$var_depth_boxplot, 
                                    "After normalization",
                                    "log10(library size)",
                                    input$select_boxplot_range1,
                                    cols =  cols)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 600, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))
       
       }) #close isolate
       }) #close boxplot_norm_libsize
     
   })# close observe generate_norm_boxplots 
   
# RUN BATCH EFFECT CORRECTION -----------------------------------------------  
   observeEvent(input$correct_be, {
     print("Pushed button to correct BE.")
     withProgress(message = 'Batch effect correction', value = 0, {
       incProgress(0.4, detail = "Correcting data...")
       
       try({
         overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
         overall_vars$countMatrices[[input$mat_name]] <- correct_be(
           overall_vars$countMatrices[[input$select_matrix_bec]], 
                                     input$be_method, 
                                     overall_vars$metadata[[input$batch_label1]], 
                                     verify_null(input$batch_label2, overall_vars$metadata), 
                                     verify_null(input$keep_biological, overall_vars$metadata),
           verify_null(input$adjust_covariate, overall_vars$metadata))
        
         overall_vars$report$bec[[input$mat_name]]$method <- input$be_method
         overall_vars$report$bec[[input$mat_name]]$batch_label1 <- input$batch_label1
         overall_vars$report$bec[[input$mat_name]]$batch_label2 <- input$batch_label2
         overall_vars$report$bec[[input$mat_name]]$biological_condition <- input$keep_biological
         
         save_session(overall_vars$session_token,overall_vars, "report")
         
         corrected_mat <- (2**overall_vars$countMatrices[[input$mat_name]]) - 1
         
         save_session(overall_vars$session_token,overall_vars, "countMatrices")
         
         overall_vars$mat_names <- names(overall_vars$countMatrices)
         
         overall_vars$countMatrices <- NULL
         
         save_session(overall_vars$session_token,overall_vars, "mat_names")
         overall_vars$md5$mat_names <- md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
         
         overall_vars$metadata[paste0("library_size_", input$mat_name)] <- colSums(corrected_mat)
         save_session(overall_vars$session_token, overall_vars, "metadata")
         overall_vars$md5$metadata <- md5sum(
           paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
         
         rm(corrected_mat)
         print(names(overall_vars$mat_names))
         
         
       })
     }) # close progress bar
   }) # close observe correct_be

## RUN BATCH EFFECT CORRECTION boxplots -----------------------------------------------        
   observeEvent(input$generate_be_boxplots, {
     output$boxplot_be <- renderPlotly({
       cols <- c(as.vector(unlist(overall_vars$colours[[input$var_be_boxplot]])), extra_cols)
       
       isolate({p <- plot_boxplot(  df = overall_vars$metadata, 
                                    var = paste0("library_size_", input$select_matrix_plot_bec), 
                                    sample = input$var_be_boxplot, 
                                    title = "After batch effect correction",
                                    ytitle = "log10(library size)",
                                    yrange = input$select_boxplot_range3,
                                    cols = cols)})
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 600, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))
     })
   }) #close generate plots
   

}#close server
  
  
 
 shinyApp(ui = ui, server = server)


