#INFORMATION-----------------------------

#LOAD LIBRARIES ------------------------
library(data.table)
library(DT)
library(dplyr)
library(ff)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(gridExtra)
library(Matrix)
library(matrixStats)
library(magrittr)
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

source("tab_DIMRED/UI.R")
source("tab_DIMRED/server.R")

source("tab_VIOLINS/UI.R")
source("tab_VIOLINS/server.R")

source("tab_DOTPLOTS/UI.R")
source("tab_DOTPLOTS/server.R")

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
  navbarPage("scStudio: Feature Plots",
          tabPanel("Dimensionality Reduction",
                     tab_DIMRED),
          tabPanel("Violin Plots",
                   tab_VIOLINS),
          tabPanel("Dot Plots",
                   tab_DOTPLOTS)
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$dotplot <- renderDataTable(data.frame())
  output$grid_violinplots <- renderDataTable(data.frame())
  output$grid_featplots <- renderDataTable(data.frame())



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
       h1("Welcome to scStudio - Feature Plots"),
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
             
             incProgress(0.5, detail = "Retrieving metadata...")
             
             overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
             overall_vars$md5$metadata <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
             
             overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
             overall_vars$md5$mat_names <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
             
             overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
             overall_vars$md5$countMatrices <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
             
             overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
             overall_vars$md5$dimred <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
            
             incProgress(0.6, detail = "Finishing...")
             
             if ("colours.rds" %in% files){
               overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
             } #close if 
             
             else{overall_vars$colours <- list()}
             
             if ("report.rds" %in% files){
               overall_vars$report<- upload_session(overall_vars$session_token, "report")
             } #close if
             else {overall_vars$report <- list()}
             
             if ("scores.rds" %in% files){
               overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
             } #close if 
             
             else{overall_vars$scores <- list()}
             
             
  
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
           
           overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
           overall_vars$md5$countMatrices <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/countMatrices.rds"))
           
           overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
           
           overall_vars$md5$genes <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
           
           overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
           overall_vars$md5$dimred <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
           
           incProgress(0.6, detail = "Finishing...")
           
           if ("colours.rds" %in% files){
             overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
           } #close if 
          
            else{overall_vars$colours <- list()}
           
           if ("report.rds" %in% files){
             overall_vars$report <- upload_session(overall_vars$session_token, "report")
           } #close if 
          
            else{overall_vars$report <- list()}
           
           if ("scores.rds" %in% files){
             overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
           } #close if 
           
           else{overall_vars$scores <- list()}
           
           overall_vars$md5$scores <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/scores.rds"))
           
        
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
      
      check_md5sum_dimred <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
      
      check_md5sum_scores <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/scores.rds"))
      
      check_md5sum_genes <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
        overall_vars$md5$metadata <- check_md5sum_metadata
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
        overall_vars$md5$mat_names <- check_md5sum_mat_names
      } # close if
      
      if (check_md5sum_countMatrices != overall_vars$md5$countMatrices){
        overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
        overall_vars$md5$countMatrices <- check_md5sum_countMatrices
      } # close if
      
      if (check_md5sum_dimred != overall_vars$md5$dimred){
        overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
        overall_vars$md5$dimred <- check_md5sum_dimred
      } # close if
      
      if (check_md5sum_scores != overall_vars$md5$scores){
        overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
        overall_vars$md5$scores <- check_md5sum_scores
      } # close if
      
      if (check_md5sum_genes != overall_vars$md5$genes){
        overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
        overall_vars$md5$genes <- check_md5sum_genes
      } # close if
      
    }) #close try
  
}) #close observe
  

   
## update inputs ----------------------------------------------------   
observe({

  updateSelectInput(session, "select_matrix_dimred",
                    choices = names(overall_vars$countMatrices), 
                    selected = c(names(overall_vars$countMatrices)[2]))
  
  updateSelectInput(session, "select_matrix_violins",
                    choices = names(overall_vars$countMatrices), 
                    selected = c(names(overall_vars$countMatrices)[2]))
  
  updateSelectInput(session, "select_matrix_dotplot",
                    choices = names(overall_vars$countMatrices)[-1], 
                    selected = c(names(overall_vars$countMatrices)[2]))
  
  updateSelectInput(session, "featplot_gsets",
                    choices = names(overall_vars$scores), 
                    selected = "")
  
  updateSelectInput(session, "select_featplot_dimred",
                    choices = names(overall_vars$dimred[[input$select_featplot_type_dimred]]))
  
  updateSelectInput(session, "select_featplot_violin_var",
                    choices = names(identify_discrete(overall_vars$metadata)), 
                    selected = "orig.ident")
  
  updateSelectInput(session, "featplot_violins_gsets",
                    choices = names(overall_vars$scores), 
                    selected = "")
  
  updateSelectInput(session, "select_featplot_violin_subvar",
                    choices = c("none", names(identify_discrete(overall_vars$metadata))), 
                    selected = "none")
  
  updateSelectInput(session, "select_featplot_dotplot_var",
                    choices = names(identify_discrete(overall_vars$metadata)), 
                    selected = "orig.ident")
  
  updateSelectInput(session, "select_featplot_dotplot_subvar",
                    choices = c("none", names(identify_discrete(overall_vars$metadata))), 
                    selected = "none")
  
  observeEvent(input$select_featplot_dimred, {
    if(input$select_featplot_type_dimred != "pca"){
      try({
        
        updateSelectInput(session, "select_featplot_xx", 
                          choices = 1:dim(overall_vars$dimred[[input$select_featplot_type_dimred]][[input$select_featplot_dimred]])[2])
        updateSelectInput(session, "select_featplot_yy", 
                          choices = 1:dim(overall_vars$dimred[[input$select_featplot_type_dimred]][[input$select_featplot_dimred]])[2],
                          selected = 2)
      })
    }
    else {
      try({
        updateSelectInput(session, "select_featplot_xx", 
                          choices = 1:dim(overall_vars$dimred[[input$select_featplot_type_dimred]][[input$select_featplot_dimred]][["x"]])[2])
        updateSelectInput(session, "select_featplot_yy", 
                          choices = 1:dim(overall_vars$dimred[[input$select_featplot_type_dimred]][[input$select_featplot_dimred]][["x"]])[2],
                          selected = 2)
      }, silent = TRUE)
    }
    
  }) # close observeEvent  
  
}) #close observe
 
   
# RUN DIMRED -----------------------------------------------  
   
observeEvent(input$select_matrix_dimred, 
             updateSelectizeInput(session, 
                                  "featplot_genes", 
                                   choices =  rownames(overall_vars$countMatrices[[input$select_matrix_dimred]]), 
                                   server = TRUE))

observeEvent(input$run_featplot, {
  
     new_mat <- overall_vars$countMatrices[[input$select_matrix_dimred]]
     
     if(!is.null(input$featplot_gsets)){
       
       for (gset in input$featplot_gsets){
         print(dim(new_mat))
         print(gset)
         new_mat <- rbind(new_mat, overall_vars$scores[[gset]])
         print(dim(new_mat))
       }
       rownames(new_mat) <- c(
         rownames(overall_vars$countMatrices[[input$select_matrix_dimred]]), unlist(input$featplot_gsets))}
     
     output$grid_featplots <- bindEvent(renderPlot({
       plot_featplot(mat = new_mat, 
                     genes = c(input$featplot_genes, input$featplot_gsets), 
                     dimred = overall_vars$dimred[[input$select_featplot_type_dimred]][[input$select_featplot_dimred]], 
                     type = input$select_featplot_type_dimred, 
                     x = as.numeric(input$select_featplot_xx), 
                     y = as.numeric(input$select_featplot_yy),
                     density_lines = input$featplot_density)
     }, height = input$featplot_dimred_height, width = input$featplot_dimred_width) , input$run_featplot) #close bindEvent

   }) 

# RUN VIOLINS ------------------------------------------------------------

observeEvent(input$select_matrix_violins, 
             updateSelectizeInput(session, 
                                  "featplot_violin_genes", 
                                  choices =  rownames(overall_vars$countMatrices[[input$select_matrix_violins]]), 
                                  server = TRUE))

observeEvent(input$run_featplot_violin, {
  print("Pushed violinplot button")
  
  new_mat <- overall_vars$countMatrices[[input$select_matrix_violins]]
  
  if(!is.null(input$featplot_violins_gsets)){
    
    for (gset in input$featplot_violins_gsets){
      print(dim(new_mat))
      print(gset)
      new_mat <- rbind(new_mat, overall_vars$scores[[gset]])
      print(dim(new_mat))
    }
    
    print(unlist(input$featplot_violins_gsets))
    
    rownames(new_mat) <- c(
      rownames(overall_vars$countMatrices[[input$select_matrix_violins]]), 
      unlist(input$featplot_violins_gsets))
    }
  
  
  output$grid_violinplots <- bindEvent(renderPlot({
    
    plot_grid_violins(
      mat = new_mat, 
      genes = c(input$featplot_violin_genes, input$featplot_violins_gsets), 
      metadata = overall_vars$metadata,
      var = input$select_featplot_violin_var,
      subvar = input$select_featplot_violin_subvar,
      cols = c(as.vector(unlist(overall_vars$colours[[input$select_featplot_violin_var]])), extra_cols),
      pt.size = input$featplot_violin_ptsize,
      gene_list = overall_vars$genes
    )
  }, height = input$featplot_violin_height, width = input$featplot_violin_width) , input$run_featplot_violin) #close bindEvent
 
})

# RUN DOTPLOT ----------------------------------------------------------
observeEvent(input$select_matrix_dotplot, 
             updateSelectizeInput(session, 
                                  "featplot_dotplot_genes", 
                                  choices =  rownames(overall_vars$countMatrices[[input$select_matrix_dotplot]]), 
                                  server = TRUE))

observeEvent(input$run_featplot_dotplot, {
  print("Pushed dotplot button")
  
  output$dotplot <- bindEvent(renderPlot({
    
    plot_dotplot(mat = overall_vars$countMatrices[[input$select_matrix_dotplot]], 
                 genes = input$featplot_dotplot_genes, 
                 var = overall_vars$metadata[,input$select_featplot_dotplot_var])
      
      
    
    
    #plot_dotplot(mat = overall_vars$countMatrices[[input$select_matrix_dotplot]], 
     #            genes = input$featplot_dotplot_genes, 
      #           gsets = input$featplot_dotplot_gsets,
       #          scores = overall_vars$scores, 
        #         metadata = overall_vars$metadata, 
         #        var = input$select_featplot_dotplot_var, 
          #       subvar = input$select_featplot_dotplot_subvar, 
           #      scale_TF = as.logical(input$featplot_dotplot_scale), 
            #     cols = c("blue", "green"), 
             #    flip = input$featplot_dotplot_reverse)
    
  }, height = input$featplot_dotplot_height, width = input$featplot_dotplot_width), input$run_featplot_dotplot)

})


}#close server
  
  
 
 shinyApp(ui = ui, server = server)


