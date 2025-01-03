#INFORMATION-----------------------------

#LOAD LIBRARIES ------------------------
library(data.table)
library(dplyr)
library(DT)
library(ff)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(gridExtra)
library(Matrix)
library(matrixStats)
library(plotly)
library(scran)
library(scuttle)
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
#library(Seurat)
#library(scater)
#library(scDblFinder)
#library(limma)

print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_FEATURE_SELECTION/UI.R")
source("tab_FEATURE_SELECTION/server.R")

source("tab_PCA/UI.R")
source("tab_PCA/server.R")

source("tab_TSNE/UI.R")
source("tab_TSNE/server.R")

source("tab_UMAP/UI.R")
source("tab_UMAP/server.R")


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
  navbarPage("scStudio: Dimensionality Reduction",
          tabPanel("Feature Selection",
                     tab_FEATURE_SELECTION),
          tabPanel("PCA",
                   tab_PCA),
          tabPanel("t-SNE",
                   tab_TSNE),
          tabPanel("UMAP",
                   tab_UMAP)
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$pca_plot <- renderDataTable(data.frame())
  output$scree_plot <- renderDataTable(data.frame())
  output$tsne_plot1 <- renderDataTable(data.frame())
  output$tsne_plot2 <- renderDataTable(data.frame())
  output$umap_plot1 <- renderDataTable(data.frame())
  output$umap_plot2 <- renderDataTable(data.frame())
  output$varvsmean <- renderDataTable(data.frame())
  output$pvalvsbio <- renderDataTable(data.frame())



  }, silent = TRUE)
  
# OVERALL SESSION VARIABLES -----------------------------------
  overall_vars <- reactiveValues(metadata = data.frame())
  
# timer  ---------------------------------------------
autoInvalidate <- reactiveTimer(15000)  

## page configuration -----------------------------------  
output$session_id <- renderText({ paste0("Session token: ", overall_vars$session_token)})

# UPLOAD SESSION ------------------------------------
## main modal -----------------------------------  
   uploadModal <- function(failed = FALSE) {
     modalDialog(
       h1("Welcome to scStudio - Dimensionality Reduction"),
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
             
             if ("genes.rds" %in% files){
               overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
             } #close if 
             else {
               overall_vars$genes <- list()
               save_session(overall_vars$session_token ,overall_vars, "genes")
             }
             
             overall_vars$md5$genes <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
             
           
             
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
           
           if ("genes.rds" %in% files){
             overall_vars$genes <- upload_session(overall_vars$session_token, "genes")
           } #close if 
           else {
             overall_vars$genes <- list()
             save_session(overall_vars$session_token ,overall_vars, "genes")
           }
           
           overall_vars$md5$genes <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
           
           
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
      
      check_md5sum_genes <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/genes.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
      } # close if
      
      if (check_md5sum_genes != overall_vars$md5$genes){
        overall_vars$genes<- upload_session(overall_vars$session_token, "genes")
        
      } # close if
    }) #close try
    
    
  md5_dimred <-  md5sum(paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
  
    try({
      
      if (md5_dimred != overall_vars$md5$dimred){
        
        overall_vars$dimred <- readRDS(
          paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
        overall_vars$md5$dimred <- md5_dimred
        

        updateSelectInput(session, "choose_pca",
                          choices = names(overall_vars$dimred$pca))
        
        updateSelectInput(session, "select_pcs_tsne",choices = names(overall_vars$dimred$pca))
        updateSelectInput(session, "select_umap_pcs",choices = names(overall_vars$dimred$pca))
      
        updateSelectInput(session, "choose_tsne1",choices = names(overall_vars$dimred$tsne))
        updateSelectInput(session, "choose_tsne2",choices = names(overall_vars$dimred$tsne))
        updateSelectInput(session, "choose_umap1",choices = names(overall_vars$dimred$umap))
        updateSelectInput(session, "choose_umap2",choices = names(overall_vars$dimred$umap))
        
        showNotification("New dimensionality reduction plot available.")
      }
    
        
    }, silent = TRUE)

  }) #close observe
## update inputs ----------------------------------------------------   
   
observe({
  updateSelectInput(session, "select_matrix_fselect",
                    choices = overall_vars$mat_names[-1])
  
  updateSelectizeInput(session, 
                       "gene_highlight", 
                       choices =  overall_vars$genes, 
                       selected = "",
                       server = TRUE)
### PCA ----------------------------------------------------------------------  
  
  
  updateSelectInput(session, "select_matrix_PCA",
                    choices = overall_vars$mat_names)
  
  updateSelectInput(session, "select_features",
                    choices = c("All",names(overall_vars$hvgs)))
  
  updateSelectInput(session, "choose_pca",
                    choices = names(overall_vars$dimred$pca))
  
    
    try({
      updateSelectInput(session, 
                        "pca_x",
                        choices = c(1:dim(overall_vars$dimred$pca[[input$choose_pca]][["x"]])[2]), 
                        selected = 1)
      
      updateSelectInput(session, 
                        "pca_y",
                        choices = c(1:dim(overall_vars$dimred$pca[[input$choose_pca]][["x"]])[2]), 
                        selected = 2)

            #updateSelectInput(session, "pca_z",choices = c(1:dim(overall_vars$dimred$pca[[input$choose_pca]][["x"]])[2]), 
      #                 selected = 3)
   }, silent = TRUE)
    
  updateSelectInput(session, "pca_color",
                      choices = c(colnames(identify_discrete(overall_vars$metadata))), 
                      selected = "orig.ident")  
    
  #updateSelectInput(session, "choose_pca",
  #                    choices = names(overall_vars$dimred$pca))
### TSNE ----------------------------------------------------------------------   
  if (input$genesvspcs == "Genes"){
    output$select_features_tsne <- renderUI(
      selectInput(inputId = "select_features_tsne", 
                  label = "Select features", 
                  choices = c("All",names(overall_vars$hvgs)), 
                  multiple = FALSE, 
                  width = "100%"))
    
    output$select_matrix_tsne <- renderUI(  
      selectInput(inputId = "select_matrix_tsne", 
                  label = "Select count matrix", 
                  choices = overall_vars$mat_names, 
                  multiple = FALSE, 
                  width = "100%"),)
    
    output$select_pcs_tsne <- NULL
    
    output$ncomponents_tsne <- NULL
    
  } #close if 
  else {
    output$select_pcs_tsne <- renderUI(
      selectInput(inputId = "select_pcs_tsne", 
                  label = "Select PCA", 
                  choices = names(overall_vars$dimred$pca), 
                  multiple = FALSE, 
                  width = "100%"))
    
    output$ncomponents_tsne <- renderUI(
      textInput(inputId = "ncomponents_tsne", 
                label = "Number of PCs:",
                value = 25))
    
    output$select_features_tsne <- NULL
    
    output$select_matrix_tsne <- NULL
  } #close else
  
  updateSelectInput(session, "choose_tsne1",
                    choices = names(overall_vars$dimred$tsne))
  
  updateSelectInput(session, "choose_tsne2",
                    choices = names(overall_vars$dimred$tsne))
  
  updateSelectInput(session, "tsne_color1",
                    choices = c(colnames(identify_discrete(overall_vars$metadata))), 
                    selected = "orig.ident")
  
  updateSelectInput(session, "tsne_color2",
                    choices = c(colnames(identify_discrete(overall_vars$metadata))), 
                    selected = "orig.ident")
  
### UMAP ----------------------------------------------------------------------  
  
  updateSelectInput(session, "choose_umap1",choices = names(overall_vars$dimred$umap))
  updateSelectInput(session, "choose_umap2",choices = names(overall_vars$dimred$umap))
  
  if (input$genesvspcs_umap == "Genes"){
    output$select_features_umap <- renderUI(
      selectInput(inputId = "select_features_umap", 
                  label = "Select features", 
                  choices = c("All",names(overall_vars$hvgs)), 
                  multiple = FALSE, width = "100%"))
    
    output$select_matrix_umap <- renderUI(
      selectInput(inputId = "select_matrix_umap", 
                  label = "Select count matrix", 
                  choices = overall_vars$mat_names, 
                  multiple = FALSE, width = "100%"))
    
    output$select_umap_pcs <- NULL
    
    output$ncomponents_umap <- NULL
    
  }
  else {
    output$select_umap_pcs <- renderUI(
      selectInput(inputId = "select_umap_pcs", 
                  label = "Select PCA", 
                  choices = names(overall_vars$dimred$pca), 
                  multiple = FALSE, 
                  width = "100%"))
    
    output$ncomponents_umap <- renderUI(
      textInput(inputId = "ncomponents_umap", 
                label = "Number of PCs:",
                value = 25))
    
    output$select_features_umap <- NULL
    
    output$select_matrix_umap <- NULL
  }
  
  
  updateSelectInput(session, "umap_color1",choices = c(colnames(identify_discrete(overall_vars$metadata))), 
                    selected = "orig.ident")
  
  updateSelectInput(session, "umap_color2",choices = c(colnames(identify_discrete(overall_vars$metadata))), 
                    selected = "orig.ident")
    
}) #close observe
 
   
# RUN FEATURE SELECTION -----------------------------------------------  
  
   observeEvent(input$run_fs, {
     withProgress(message = 'Calculating variance', value = 0, {
       incProgress(0.4, detail = "Estimating biological variability and technical noise.")
       
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
       
       hvgs <- calculate_hvg(overall_vars$countMatrices[[input$select_matrix_fselect]])
       
       overall_vars$countMatrices <- NULL
       
       overall_vars$report$hvgs[[input$hvgs_id]]$table <- hvgs
       overall_vars$report$hvgs[[input$hvgs_id]]$matrix <- input$select_matrix_fselect
       overall_vars$report$hvgs[[input$hvgs_id]]$top<- input$top_hvg
       overall_vars$report$hvgs[[input$hvgs_id]]$pval<- input$max_pval
       
       save_session(overall_vars$session_token, overall_vars, "report")
       
       output$varvsmean<- renderPlotly({
         p <- plot_varvsmean(hvgs, 
                             as.numeric(input$top_hvg), 
                             as.numeric(input$max_pval),
                             input$gene_highlight)  
         
         toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "Gene") %>% 
                   layout(dragmode = "lasso", height = 500, 
                          legend = list(orientation = "h", x = 0, y = -0.2)))
       })
       
       output$pvalvsbio <- renderPlotly({
          p <- plot_pvalvsbio(hvgs,
                              as.numeric(input$top_hvg), 
                              as.numeric(input$max_pval),
                              input$gene_highlight)
          
          toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "Gene") %>% 
                     layout(dragmode = "lasso", height = 500, 
                            legend = list(orientation = "h", x = 0, y = -0.2)))
          })
       
       #output$hvgs_table <- DT::renderDataTable({ 
      #     
      #     DT::datatable(
      #       { hvgs},
      #       extensions = 'Buttons',
      #       rownames= FALSE,
      #       options = list(
      #         paging = TRUE,
      #         searching = TRUE,
      #         fixedColumns = FALSE,
      #         autoWidth = FALSE,
      #         ordering = TRUE,
      #         dom = 'Bfrtip',
      #         lengthMenu=c(10,20,50),
      #         buttons = c('copy', 'csv', 'excel')
      #       ),
      #       class = "display"
      #     )
      # }) #close table hvgs
       
     }) # close progress bar
   }) # close observe feature_selection
   
   observeEvent(input$subset_features, {

     hvgs <- overall_vars$report$hvgs[[input$hvgs_id]]$table

     include_genes <- c(rep(TRUE, input$top_hvg), rep(FALSE, dim(hvgs)[1] - as.numeric(input$top_hvg))) &
       hvgs$p.value <= as.numeric(input$max_pval)
     
     include_genes <- rownames(hvgs)[include_genes]
     
     
     
     overall_vars$hvgs[[input$hvgs_id]] <- include_genes

     
     updateSelectInput(session, "select_features",choices = c("All",names(overall_vars$hvgs)))
     save_session(overall_vars$session_token, overall_vars, "hvgs")

     showNotification("New subset of features available.")
     
     
   }) # close observe subset_features
   
# RUN PCA -------------------------------------------------------------------
   observeEvent(input$run_pca, {
     
     withProgress(message = 'Preparating PCA calculation.', value = 0, {
       incProgress(0.4, detail = "Accessing count data.")
       
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
       
     if (input$select_features == "All"){
       subset_row <- rownames(overall_vars$countMatrices[[input$select_matrix_PCA]])} 
     else {subset_row <- overall_vars$hvgs[[input$select_features]]}
     
     
     }) # close progress bar
     
     showNotification("PCA is running in the background. You will be notified when the job finishes.")
     
     print("Running PCA...")
     overall_vars$jobs[[input$pca_id]] <- callr::r_bg(getPCA, 
                                                      args = list(mat = overall_vars$countMatrices[[input$select_matrix_PCA]], 
                                                                          ncomponents = as.numeric(input$ncomponents), 
                                                                          subset_row = subset_row, 
                                                                          do_scaling = input$do_scaling, 
                                                                          token =  overall_vars$session_token, 
                                                                          session_obj = overall_vars$dimred, 
                                                                          ID = input$pca_id))
     overall_vars$countMatrices <- NULL
     
     overall_vars$report$pca[[input$pca_id]]$matrix <- input$select_matrix_PCA
     overall_vars$report$pca[[input$pca_id]]$ncomponents <- input$ncomponents
     overall_vars$report$pca[[input$pca_id]]$features <- input$select_features
     overall_vars$report$pca[[input$pca_id]]$scaling <- input$do_scaling
     
    save_session(overall_vars$session_token, overall_vars, "report")
   }) # close observe run_pca
   
   observeEvent(input$plot_pca, {
     try({
       
       cols <- c(get_ordered_colors(overall_vars$colours[[input$pca_color]]), extra_cols)
       
       output$pca_plot <- renderPlotly({
         isolate({p <- plot_dimRed(coord = overall_vars$dimred$pca[[input$choose_pca]][["x"]], 
                                   x = as.numeric(input$pca_x), 
                                   y = as.numeric(input$pca_y), 
                                   metadata = overall_vars$metadata, 
                                   sel_var = input$pca_color, 
                                   title = input$pca_title,
                                   pointSize = input$pca_pointSize,
                                   cols = cols,
                                   density_lines = input$pca_density)
         
         toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                   layout(dragmode = "lasso", height = 800, 
                          legend = list(orientation = "h", x = 0, y = -0.2)))}) 
         
       }) #close pca_plot
     }, silent = FALSE) #close try
     
     output$scree_plot <- renderPlotly({
       isolate(p <- plot_scree(overall_vars$dimred$pca[[input$choose_pca]]))
     })
   }) # close observe plot_pca
   
# RUN TSNE -------------------------------------------------------------------
   
   observeEvent(input$run_tsne, {
     
     withProgress(message = 'Preparating t-SNE calculation.', value = 0, {
       incProgress(0.4, detail = "Accessing count data.")
       
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, 
                                                    "countMatrices") 
       
     }) # close progress bar
     
     showNotification("TSNE is running in background. You will be notified when the job finishes.")

     if (input$genesvspcs == "PCs"){
       print("Running TSNE with PCs.")
       overall_vars$jobs[[input$tsne_id]] <- callr::r_bg(getTSNE,
                                                         args = list(
                                                           mat = t(overall_vars$dimred$pca[[input$select_pcs_tsne]][["x"]]),
                                                           subset_row = c(1:as.numeric(input$ncomponents_tsne)),
                                                           perplexity = input$perplexity,
                                                           token =  overall_vars$session_token, 
                                                           session_obj = overall_vars$dimred, 
                                                           ID = input$tsne_id))
       print("TSNE job submitted successfully.")
       overall_vars$report$tsne[[input$tsne_id]]$method <- "PCA"
       overall_vars$report$tsne[[input$tsne_id]]$matrix <- input$select_pcs_tsne
       overall_vars$report$tsne[[input$tsne_id]]$features <- input$ncomponents_tsne
       overall_vars$report$tsne[[input$tsne_id]]$perplexity <- input$perplexity
       print("Report for new tsne was updated.")
       
     }
     
     else {
       if (input$select_features_tsne == "All"){

         subset_row <- rownames(overall_vars$countMatrices[[input$select_matrix_tsne]])} 

       else {subset_row <- overall_vars$hvgs[[input$select_features_tsne]]}
       
       overall_vars$jobs[[input$tsne_id]] <- callr::r_bg(getTSNE,
                                                         args = list(
                                                           mat = overall_vars$countMatrices[[input$select_matrix_tsne]],
                                                           subset_row = overall_vars$hvgs[[input$select_features_tsne]],
                                                           perplexity = input$perplexity,
                                                           token =  overall_vars$session_token, 
                                                           session_obj = overall_vars$dimred, 
                                                           ID = input$tsne_id))
       
       overall_vars$report$tsne[[input$tsne_id]]$method <- "Genes"
       overall_vars$report$tsne[[input$tsne_id]]$matrix <- input$select_matrix_tsne
       overall_vars$report$tsne[[input$tsne_id]]$features <- input$select_features_tsne
       overall_vars$report$tsne[[input$tsne_id]]$perplexity <- input$perplexity
     }
   }) # close run_tsne
   
   observeEvent(input$plot_tsne1, {
     cols <- c(get_ordered_colors(overall_vars$colours[[input$tsne_color1]]), extra_cols)
     
     output$tsne_plot1 <- renderPlotly({
       isolate({ p <- plot_dimRed(coord = overall_vars$dimred$tsne[[input$choose_tsne1]], 
                                  x = 1, y = 2, 
                                  metadata = overall_vars$metadata, 
                                  sel_var = input$tsne_color1, 
                                  title = input$tsne_title1,
                                  pointSize = input$tsne1_pointSize,
                                  cols = cols,
                                  density_lines = input$tsne1_density)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 800, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))}) 
     }) #close plot_tsne1
   }) #close observe plot_tsne1 
   
   observeEvent(input$plot_tsne2, {
     cols <- c(get_ordered_colors(overall_vars$colours[[input$tsne_color2]]), extra_cols)
     
     output$tsne_plot2 <- renderPlotly({
       isolate({ p <- plot_dimRed(coord = overall_vars$dimred$tsne[[input$choose_tsne2]], 
                                  x = 1, y = 2, 
                                  metadata = overall_vars$metadata, 
                                  sel_var = input$tsne_color2, 
                                  title = input$tsne_title2,
                                  pointSize = input$tsne2_pointSize,
                                  cols = cols,
                                  density_lines = input$tsne2_density)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 800, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))}) 
     }) #close plot_tsne2
   }) #close observe plot_tsne2    
   
# RUN UMAP ------------------------------------------------------------------- 
   observeEvent(input$run_umap, {
     
     withProgress(message = 'Preparating UMAP calculation.', value = 0, {
       incProgress(0.4, detail = "Accessing count data.")
       
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, 
                                                    "countMatrices") 
       
     }) # close progress bar
     
     showNotification("UMAP is running in background - you will be notified when the job finishes.")
     
     if (input$genesvspcs_umap == "PCs"){
       print("Running UMAP with PCs.")
       overall_vars$jobs[[input$umap_id]] <- callr::r_bg(getUMAP,
                                                         args = list(
                                                           mat = t(overall_vars$dimred$pca[[input$select_umap_pcs]][["x"]]),
                                                           subset_row = c(1:as.numeric(input$ncomponents_umap)),
                                                           min_dist = as.numeric(input$min_dist_umap),
                                                           nneigh = as.numeric(input$n_neighbors_umap),
                                                           token =  overall_vars$session_token, 
                                                           session_obj = overall_vars$dimred, 
                                                           ID = input$umap_id))
       
       overall_vars$report$umap[[input$umap_id]]$method <- "PCA"
       overall_vars$report$umap[[input$umap_id]]$matrix <- input$select_umap_pcs
       overall_vars$report$umap[[input$umap_id]]$features <- input$ncomponents_umap
       overall_vars$report$umap[[input$umap_id]]$min_dist <- as.numeric(input$min_dist_umap)
       overall_vars$report$umap[[input$umap_id]]$nneigh <- as.numeric(input$n_neighbors_umap)
       
       save_session(overall_vars$session_token, overall_vars, "report")
     }
     
     else {
       if (input$select_features_umap == "All"){
         subset_row <- rownames(overall_vars$countMatrices[[input$select_matrix_umap]])} 
       else {subset_row <- overall_vars$hvgs[[input$select_features_umap]]}
       
       overall_vars$jobs[[input$umap_id]] <- callr::r_bg(getUMAP,
                                                         args = list(
                                                           mat = overall_vars$countMatrices[[input$select_matrix_umap]],
                                                           subset_row = overall_vars$hvgs[[input$select_features_umap]],
                                                           min_dist = as.numeric(input$min_dist_umap),
                                                           nneigh = as.numeric(input$n_neighbors_umap),
                                                           token =  overall_vars$session_token, 
                                                           session_obj = overall_vars$dimred, 
                                                           ID = input$umap_id))
       
       overall_vars$report$umap[[input$umap_id]]$method <- "Genes"
       overall_vars$report$umap[[input$umap_id]]$matrix <- input$select_matrix_umap
       overall_vars$report$umap[[input$umap_id]]$features <- input$select_features_umap
       overall_vars$report$umap[[input$umap_id]]$min_dist <- as.numeric(input$min_dist_umap)
       overall_vars$report$umap[[input$umap_id]]$nneigh <- as.numeric(input$n_neighbors_umap)
       
       save_session(overall_vars$session_token, overall_vars, "report")
     }
   }) # close run_umap
   
   observeEvent(input$plot_umap1, {
     cols <- c(get_ordered_colors(overall_vars$colours[[input$umap_color1]]), extra_cols)
     output$umap_plot1 <- renderPlotly({
       isolate({ p <- plot_dimRed(coord = overall_vars$dimred$umap[[input$choose_umap1]], 
                                  x = 1, y = 2, 
                                  metadata = overall_vars$metadata, 
                                  sel_var = input$umap_color1, 
                                  title = input$umap_title1,
                                  pointSize = input$umap1_pointSize,
                                  cols = cols,
                                  density_lines = input$umap1_density)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 800, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))}) 
     }) #close umap_plot1
   }) #close observe umap_plot1
   
   observeEvent(input$plot_umap2, {
     cols <- c(as.vector(unlist(overall_vars$colours[[input$umap_color2]])), extra_cols)
     output$umap_plot2 <- renderPlotly({
       isolate({ p <- plot_dimRed(coord = overall_vars$dimred$umap[[input$choose_umap2]], 
                                  x = 1, y = 2, 
                                  metadata = overall_vars$metadata, 
                                  sel_var = input$umap_color2, 
                                  title = input$umap_title2,
                                  pointSize = input$umap2_pointSize,
                                  cols = cols,
                                  density_lines = input$umap2_density)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 800, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))}) 
     }) #close umap_plot2
   })  #close observe umap_plot2
   

   
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


