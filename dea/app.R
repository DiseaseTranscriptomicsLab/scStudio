#INFORMATION-----------------------------

# updated on the appserver 

#LOAD LIBRARIES ------------------------
library(data.table)
library(DT)
library(dplyr)
library(ff)
library(ggheatmap) #install dev github version
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

#library(markdown)
#library(scico)
#library(dqshiny)
#library(scater)
#library(scDblFinder)
#library(scran)
#library(limma)

print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_DEA/UI.R")
source("tab_DEA/server.R")

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
  navbarPage("scStudio: DEA",
          tabPanel("Gene Ranking",
                     tab_DEA)
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {
  
# anchors ---------------------------------
  try({
  
  output$volcano_dea <- renderDataTable(data.frame())
  output$dea_plot <- renderDataTable(data.frame())
  output$heatmap_dea <- renderDataTable(data.frame())
  output$table_dea <- renderDataTable(data.frame())




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
       h1("Welcome to scStudio - Differential Expression Analysis"),
       p("In this tab, you will rank genes based on their power to discern conditions.
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
             
             vars <- colnames(overall_vars$metadata) ###IMPROVE
             options <- c()
             for (var in vars){
               if (!is.numeric(overall_vars$metadata[[var]])) {
                 options <- c(options, var)} #close if
             } #close for
             
             updateSelectInput(session, "dea_condition",choices = options)
             
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
           
           if ("dea.rds" %in% files){
             overall_vars$dea <- upload_session(overall_vars$session_token, "dea")
           } #close if 
           else {overall_vars$dea <- list()
           save_session(overall_vars$session_token ,overall_vars, "dea")}
           
           overall_vars$md5$dea <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
           
           vars <- colnames(overall_vars$metadata) ###IMPROVE
           options <- c()
           for (var in vars){
             if (!is.numeric(overall_vars$metadata[[var]])) {
               options <- c(options, var)} #close if
           } #close for
           
           updateSelectInput(session, "dea_condition",choices = options)
           
           
           
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
    
  md5_hvgs <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/hvgs.rds"))
  
  try({
    
    if (md5_hvgs != overall_vars$md5$hvgs){
      
      overall_vars$hvgs <- readRDS(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/hvgs.rds"))
      overall_vars$md5$hvgs <- md5_hvgs
      
    #
      #
      #
      
    }
    
  }, silent = FALSE)
  
  
  try({
    
    md5_dea <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
    if (md5_dea != overall_vars$md5$dea){
      print("Updating DEA results...")
      
      overall_vars$dea <- readRDS(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/dea.rds"))
      
      overall_vars$md5$dea <- md5_dea
      
      updateSelectInput(session, "select_dea",choices = names(overall_vars$dea))
   
      showNotification("New DEA available.")
    }
  }, silent = TRUE)
  
}) #close observe
  

   
## update inputs ----------------------------------------------------   
observe({
  
  updateSelectInput(session, "select_matrix_dea",choices = names(overall_vars$countMatrices)[-1])
  
  updateSelectInput(session, "select_dea",choices = names(overall_vars$dea))
  
  
  try({
    
    options <- c()
    
    for (condition in input$dea_condition){
      options <- c(options, as.character(overall_vars$metadata[,condition]))
    }
    
    updateSelectInput(session, "dea_group1",choices = options, selected = "")
    
    updateSelectInput(session, "dea_group2",choices = options, selected = "")
  })
  
  isolate(output$table_dea <- DT::renderDataTable({ 
    
    table <- overall_vars$dea[[input$select_dea]]
    
    table <- table[order(table[["Log2FC (mean)"]], decreasing = TRUE),]
    
    DT::datatable({ table },
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
  )#close isolate

 
}) #close observe
 
   
# RUN DEA -----------------------------------------------  

observeEvent(input$run_dea, {

  withProgress(message = 'Preparing data...', value = 0, {
    incProgress(0.4, detail = "")
     
     group1 <- get_groups_dea(input$dea_group1, input$dea_grouping1, overall_vars$metadata)
     group2 <- get_groups_dea(input$dea_group2, input$dea_grouping2, overall_vars$metadata)
     
     print("Number of cells on group 1:")
     print(length(group1))
     print("Number of cells on group 2:")
     print(length(group2))
    
     overall_vars$jobs[[input$dea_id]] <- callr::r_bg(get_dea_allMethods,
                                                      args = list(
                                                      input$dea_id, 
                                                      overall_vars$countMatrices[[input$select_matrix_dea]], 
                                                      group1, 
                                                      group2, 
                                                      input$select_method_dea,
                                                      token = overall_vars$session_token
                                                      ), supervise = FALSE)
     }) # close progress bar
     
  showNotification("Running DEA in the background.")
     
  overall_vars$report$DEA[[input$dea_id]][["matrix"]] <- input$select_matrix_dea
  overall_vars$report$DEA[[input$dea_id]][["method"]] <- input$select_method_dea
  overall_vars$report$DEA[[input$dea_id]][["group1"]] <- input$dea_group1
  overall_vars$report$DEA[[input$dea_id]][["group1_type"]] <- input$dea_grouping1
  overall_vars$report$DEA[[input$dea_id]][["group1_index"]] <- group1
  overall_vars$report$DEA[[input$dea_id]][["group2"]] <- input$dea_group2
  overall_vars$report$DEA[[input$dea_id]][["group2_type"]] <- input$dea_grouping2
  overall_vars$report$DEA[[input$dea_id]][["group2_index"]] <- group2
  
  save_session(overall_vars$session_token, overall_vars, "report")
  print("Report saved successfully.")
     
   }) #close run_dea
   
observeEvent(input$get_volcano, {
  try({
  
  ranking <- overall_vars$dea[[input$select_dea]]
  ranking[["Adjusted p-value"]] <- -log10(ranking[["Adjusted p-value"]] +10**-100)
  
  gene_label <- rep("cell", dim(ranking)[1])
  
  
  gene_label[ranking[["Adjusted p-value"]] > input$volcano_pvalue &
               ranking[["Log2FC (mean)"]] >= input$volcano_fc] <- "UP"
  
  gene_label[ranking[["Adjusted p-value"]] > input$volcano_pvalue &
               ranking[["Log2FC (mean)"]] <= -input$volcano_fc] <- "DOWN"
  
  output$volcano_dea <- renderPlotly({
    
    req(input$select_dea)
    
    isolate(p <- make_volcano(overall_vars$dea[[input$select_dea]], 
                              input$select_dea, 
                              gene_label, 
                              as.numeric(input$volcano_pvalue),
                              logFC = input$volcano_fc)
    )#close isolate
    
    toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
              layout(dragmode = "lasso", height = 400, 
                     legend = list(orientation = "h", x = 0, y = -0.2)))
  })
}, silent = FALSE) #close try
  
}) #close get_volcano   
   

observeEvent(input$get_heatmap,{
  print("Generating heatmap")
  output$heatmap_dea <- renderPlot({
    print("Getting top genes")
    genes <- get_top_genes(overall_vars$dea[[input$select_dea]], 
                           as.numeric(input$heatmap_topGenes))
    print("Getting ID")
    
    id <- remove_method_from_id(input$select_dea)
    
    print("plotting")
    isolate(h <- make_heatmap(mat = overall_vars$countMatrices[[overall_vars$report$DEA[[id]][["matrix"]]]],
                      genes = genes, 
                      group1 = overall_vars$report$DEA[[id]][["group1_index"]], 
                      group2 = overall_vars$report$DEA[[id]][["group2_index"]],
                      cluster = input$cluster_heatmap,
                      scale = input$dea_scale_heatmap)
  )# close isolate
    
       h
       
       }, 
    height = 650, width = 1250)
  
})#close observeEvent
   
   
observeEvent(input$get_dea_plot,{
  
  output$dea_plot <- renderPlotly({
       
    isolate(p <- make_dea_scatter(overall_vars$dea[[input$select_dea]], 
                                  input$dea_cutoff)
            )#close isolate
       
    toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 600, width = 500, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))
    })
  })# close observeEvent
}#close server
  
  
 
 shinyApp(ui = ui, server = server)


