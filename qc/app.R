#INFORMATION-----------------------------

# add doublet methods
# allow to store doublet analyses 
# add doublet analysis to report 
# add cut off line to violins. 

#LOAD LIBRARIES ------------------------
library(bslib)
library(data.table)
library(dplyr)
library(DT)
library(ff)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggridges)
library(gplots)
library(gridExtra)
library(Matrix)
library(matrixStats)
library(limma)
library(plotly)
library(scater)
library(scDblFinder)
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

#library(markdown)
#library(scico)
#library(dqshiny)
#library(reshape2)


print("Sucessfully loaded libraries.")

#LOAD TABS-------------------------------

source("setup.R")

source("tab_QC_CELLS/UI.R")
source("tab_QC_CELLS/server.R")

source("tab_QC_GENES/UI.R")
source("tab_QC_GENES/server.R")

source("tab_QC_BARPLOTS/UI.R")
source("tab_QC_BARPLOTS/server.R")

source("tab_QC_DOUBLETS/UI.R")
source("tab_QC_DOUBLETS/server.R")

print("Successfully loaded tabs.")

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
  navbarPage("scStudio: Quality Control",
          tabPanel("Cells",
                     tab_QC_CELLS      
             ), #close QC_CELLS
           tabPanel("Genes",
                    tab_QC_GENES      
           ), #close QC_GENES
           tabPanel("Doublets",
                    tab_QC_DOUBLETS     
           ), #close QC_GENES
          tabPanel("Barplots",
                   tab_QC_BARPLOTS     
          ) #close QC_GENES
      )#close navbarPage 
  
) #close UI 


# SERVER ------------------------------------
server <- function(input, output, session) {

# anchors ---------------------------------
  try({
  
  output$summary_table <- renderDataTable(data.frame())
  output$summary_barplot <- renderPlotly(data.frame())
  output$discarded_features <- renderPlot(data.frame())
  output$doublet_plots<- renderPlotly(data.frame())
  output$boxplot_libraries <- renderPlotly(data.frame())
  output$boxplot_features <- renderPlotly(data.frame())
  output$boxplot_mtgenes <- renderPlotly(data.frame())

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
       h1("Welcome to scStudio - QC"),
       p("In this tab, you will be able to perform quality control of your dataset.
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
             
             incProgress(0.6, detail = "Finishing...")
             
             if ("gene_qc.rds" %in% files){
               overall_vars$gene_qc <- upload_session(overall_vars$session_token, "gene_qc")
             } #close if 
             
             else{overall_vars$gene_qc <- list()}
             
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
           overall_vars$md5$metadata <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
           
           incProgress(0.6, detail = "Finishing...")
          
            if ("gene_qc.rds" %in% files){
             overall_vars$gene_qc <- upload_session(overall_vars$session_token, "gene_qc")
           } #close if 
           
           else{overall_vars$gene_qc <- list() } #close else
           
           if ("doublets.rds" %in% files){
             overall_vars$doublets <- upload_session(overall_vars$session_token, "doublets")
           } #close if 
           
           else{overall_vars$doublets <- list()}
           
           if ("colours.rds" %in% files){
             overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
           } #close if 
          
            else{overall_vars$colours <- data.frame()}
           
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
      # check if metadata file changed
      check_md5sum_metadata <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
      
      # update metadata if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
      } # close if
    }) #close try
    
## CELL QC  ----------------------------------------------------------   
   try({
     req(input$ok_upload_session)
     if (is.null(overall_vars$metadata$library_size)){
       
### calculate QC (CELL AND GENE) metrics ---------------------------------------------       
       print("Calculating QC metrics.")
       withProgress(message = 'Calculating QC metrics', value = 0, {
         
       incProgress(0.4, detail = "Preparing data")
         
       # we calculate everything after we upload the count matrix to save time
       mat <- upload_session(overall_vars$session_token, "countMatrices") 
       mat <- mat[["rawCountMatrix"]]  
       
         incProgress(0.4, detail = "Finishing calculations")
         
         try({
           QCmetrics <- calculateQCmetrics(mat)
           overall_vars$metadata <- cbind(overall_vars$metadata, QCmetrics)
        
           overall_vars$gene_qc[["metrics"]] <- calculate_gene_qc(mat)
  
           cols <- c(as.vector(unlist(overall_vars$colours[[input$features_sample]])), extra_cols)
          
            overall_vars$gene_qc[["top_genes"]] <- 
              get_top_expressed_genes(mat, overall_vars$metadata) 
            
            
         }) #close try
         
         incProgress(0.1, detail = "Preparing plots")
         
       }) # close progress bar
     } #close if 

### update sliders ----------------------------------------------------        
   scale <-  input$selectScale
   
   updateSliderInput(session, "selectLibSizes",
                     value = c(round(select_vline(min(overall_vars$metadata$library_size), 
                                                  scale),2), 
                               round(select_vline(max(overall_vars$metadata$library_size),
                                                  scale),2)),
                      min = round(select_vline(min(overall_vars$metadata$library_size), 
                                              scale),2), 
                      max = round(select_vline(max(overall_vars$metadata$library_size),
                                              scale),2), 
                      step = if(scale == "Linear"){5} else{0.1})
#***** 
   updateSliderInput(session, "selectNrFeatures",
                     value = c(round(select_vline(min(overall_vars$metadata$total_features), 
                                                  scale),2), 
                              round(select_vline(max(overall_vars$metadata$total_features),
                                                  scale),2)),
                     min = round(select_vline(min(overall_vars$metadata$total_features), 
                                              scale),2), 
                     max = round(select_vline(max(overall_vars$metadata$total_features),
                                              scale),2), 
                     step = if(scale == "Linear"){5} else{0.1})
#****   
   updateSliderInput(session, "selectMtCutoff",
                     value = round(as.vector(quantile(overall_vars$metadata$percentage_mt_genes, 0.99)),2),
                     min = round(min(overall_vars$metadata$percentage_mt_genes),2), 
                     max = round(max(overall_vars$metadata$percentage_mt_genes),2), 
                     step = 0.1)
   
   updateSelectInput(session, "features_subsample",
                     choices = unique(overall_vars$metadata[[input$features_sample]]))
   
## histogram - library size ---------------------------------------------------
   output$histLibSize <- renderPlotly({ 
     # make plots wait for initialized values
     if (input$selectLibSizes[1] != -1){
       
       if (!is.null(input$features_subsample)){
         metadata_df <- overall_vars$metadata[overall_vars$metadata[[input$features_sample]] %in% input$features_subsample,]
       } 
       else{metadata_df <- overall_vars$metadata}
       
       if (input$selectScale == "Log2"){
         x_legend = "Log2(Total counts)"
       }
       else if (input$selectScale == "Log10"){
         x_legend = "Log10(Total counts)"
       }
       else(x_legend = "Total counts")
       
     p <- plot_histogram(metadata_df, 
                         "library_size", 
                         input$selectLibBins, 
                         title = "Library size",
                         xlab = x_legend,
                         min = input$selectLibSizes[1],
                         max = input$selectLibSizes[2],
                         input$selectScale)
       
     toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 400, 
                        legend = list(orientation = "h", x = 0, y = -0.2)))
        } #close if 
       }) #close plot_histLibSize  
    
## histogram - number of features ---------------------------------------------------
   output$histFeatures <- renderPlotly({
     if (input$selectLibSizes[1] != -1){
       
       if (!is.null(input$features_subsample)){
         metadata_df <- overall_vars$metadata[overall_vars$metadata[[input$features_sample]] %in% input$features_subsample,]
       } 
       else{metadata_df <- overall_vars$metadata}
       
       if (input$selectScale == "Log2"){
         x_legend = "Log2(Total features)"
       }
       else if (input$selectScale == "Log10"){
         x_legend = "Log10(Total features)"
       }
       else(x_legend = "Total features")
       
      p <- plot_histogram(metadata_df,"total_features", 
                          input$selectFeatBins, 
                          title = "Unique features",
                          xlab = x_legend,
                          min = input$selectNrFeatures[1],
                          max = input$selectNrFeatures[2],
                          input$selectScale)
       
       toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 400, 
                        legend = list(orientation = "h", x = 0, y = -0.2))) 
       } #close if 
      }) #close histFeatures
   
## scatter plot - library size vs number of features -------------------------------------------
   output$scatterLibFeat <- renderPlotly({
     if (input$selectLibSizes[1] != -1){
 
       if (!is.null(input$features_subsample)){
         metadata_df <- overall_vars$metadata[overall_vars$metadata[[input$features_sample]] %in% input$features_subsample,]
       } 
       else{metadata_df <- overall_vars$metadata}

      p <- plot_scatterLibFeat(metadata_df, 
                                input$selectScale,
                                input$selectLibSizes[1],
                                input$selectLibSizes[2],
                                input$selectNrFeatures[1],
                                input$selectNrFeatures[2],
                                leg = "% MT-genes",
                                midpt = input$selectMtCutoff)
     
     toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
               layout(dragmode = "lasso", height = 400, 
                      legend = list(orientation = "h", x = 0, y = -0.2))) 
       } #close if 
      }) #close scatterLibFeat
   
## scatter plot - library size vs % mt-genes -------------------------------------------
   output$scatterLibMt <- renderPlotly({ 
     if (input$selectLibSizes[1] != -1){
       
       if (!is.null(input$features_subsample)){
         metadata_df <- overall_vars$metadata[overall_vars$metadata[[input$features_sample]] %in% input$features_subsample,]
       } 
       else{metadata_df <- overall_vars$metadata}
       
     p <- plot_scatterLibMt(metadata_df, 
                            input$selectScale,
                            input$selectLibSizes[1],
                            input$selectLibSizes[2],
                            leg = "% MT-genes",
                            midpt = input$selectMtCutoff)
     
     toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
               layout(dragmode = "lasso", height = 400, 
                      legend = list(orientation = "h", x = 0, y = -0.2))) 
       } #close if 
      }) #close scatterLibMt
   
  }) #close try

# GENE QC  ----------------------------------------------------------       
    
    #updateSelectInput(session, "top_genes_sample",
    #                   choices = names(identify_discrete(overall_vars$metadata)),
    #                   selected = "orig.ident")
   
   try({
     filter <- overall_vars$gene_qc$metrics$total_exp > as.numeric(input$select_nr_counts) & 
       overall_vars$gene_qc$metrics$expr_cells > as.numeric(input$select_nr_cells)
     
   

## scatter plot - number of expressing cells vs total counts -------------------------------------
     output$discarded_features <- renderPlotly({
       p <- plot_scatter_features(!filter, 
                                  filterx = as.numeric(input$select_nr_counts), 
                                  filtery = as.numeric(input$select_nr_cells), 
                                  gene_qc = overall_vars$gene_qc$metrics)
       
       toWebGL(ggplotly(p, source = "gene_qc", key = "pointNumber", tooltip = "text") %>% 
                 layout(dragmode = "lasso", height = 500, width = 500)) 
     })
 
## scatter plot - top expressed genes -------------------------------------         
   # observeEvent( input$top_genes_n ,{
      output$topGenes <- renderPlot({
       #p <- overall_vars$gene_qc$top_genes
       
       # need to reorder the colors according to the order of plotting the points
       
       #p$data$colour_by <- rep(overall_vars$metadata[[input$top_genes_sample]], each = 30)
  
       #n <- as.numeric(input$top_genes_n)

       #p + ylim(rev(p$plot_env$sub_names[1:n]))
        
      plot_top_genes(df = overall_vars$gene_qc$top_genes, 
                     log_scale = input$top_genes_scale,
                     nr_genes = input$top_genes_n)
      
    
      
     # toWebGL(ggplotly(p, source = "gene_qc", key = "pointNumber", tooltip = "text") %>% 
       #         layout(dragmode = "lasso", height = 800, width = 500)) 
     
     }) 
    #}) #close observe top_genes_n    
       
   }) # close try
   
# DOUBLETS -----------------------------------------------------
  
   updateSelectInput(session, "doublet_channel",
                     choices = names(identify_discrete(overall_vars$metadata)),
                     selected = "orig.ident")
    
    output$doublet_plots <- renderPlotly({ 
      req(overall_vars$metadata$doublets_class)
      g <- arrange_doublet_plots(df_qc = overall_vars$metadata , 
                                   scale = input$selectScale,
                                   sample = input$doublet_channel) 
      
      toWebGL(ggplotly(g, source = "main", key = "pointNumber", tooltip = "text") %>% 
                layout(dragmode = "lasso", height = 600, 
                       legend = list(orientation = "h", x = 0, y = -0.2))) 
      }) #close doublet_plots 
    
    output$doublet_metrics <- renderDataTable({
      req(overall_vars$metadata$doublets_class)
      get_doublets_table(input$doublet_channel, overall_vars$metadata)
     }) #close doublet_metrics table
    
    try({
 
      overall_vars$metadata$doublets_class[overall_vars$metadata$doublets_score >= input$selectDoubletCutoff] <- "doublet"
      overall_vars$metadata$doublets_class[overall_vars$metadata$doublets_score < input$selectDoubletCutoff] <- "singlet"

    }, silent = FALSE) #close try
    
# BARPLOTS -----------------------------------------------------    
    updateSelectInput(session, "summary_var",
                      choices = names(identify_discrete(overall_vars$metadata)),
                      selected = "orig.ident")    

}) #close observe
 
   
# RUN CELL QC -----------------------------------------------  

## boxplots -------------------------------------------------
   # define colors for boxplots
   observeEvent(input$update_violins, {
     
   cols <- c(as.vector(unlist(overall_vars$colours[[input$features_sample]])), extra_cols)
   try({
   ### violinplot for library sizes ----------------------------------------------------------------
     
     if (!is.null(input$features_subsample)){
       metadata_df <- overall_vars$metadata[overall_vars$metadata[[input$features_sample]] %in% input$features_subsample,]
     } 
     else{metadata_df <- overall_vars$metadata}
     
  output$boxplot_libraries <- renderPlotly({ isolate(p <- plot_boxplot(df = metadata_df, 
                                                               var = "library_size", 
                                                               sample = input$features_sample, 
                                                               title = "",
                                                               ytitle = "Library size",
                                                               cols = cols,
                                                               scale = input$selectScale,
                                                               showCells = input$violin_qc_showDots))
   
   toWebGL(ggplotly(p, source = "gene_qc", key = "pointNumber", tooltip = "text") %>% 
             layout(dragmode = "lasso", height = 500, width = 500,
                    legend = list(orientation = "h", x = 0, y = -0.2))) 
   }) #close violinplot_features
   
   ### violinplot for number of unique features --------------------------------------------------------
   
   output$boxplot_features <- renderPlotly({ p <- isolate(plot_boxplot(df = metadata_df, 
                                                               var = "total_features", 
                                                               sample = input$features_sample, 
                                                               title = "",
                                                               ytitle = "Unique features",
                                                               cols = cols,
                                                               scale = input$selectScale,
                                                               showCells = input$violin_qc_showDots))
   
   toWebGL(ggplotly(p, source = "gene_qc", key = "pointNumber", tooltip = "text") %>% 
             layout(dragmode = "lasso", height = 500, width = 500,
                    legend = list(orientation = "h", x = 0, y = -0.2))) 
   }) #close violinplot_features
   
   ###boxplot for %mt-genes ----------------------------------------------------------------
   output$boxplot_mtgenes <- renderPlotly({ p <- isolate(plot_boxplot(df = metadata_df, 
                                                               var = "percentage_mt_genes", 
                                                               sample = input$features_sample, 
                                                               title = "",
                                                               ytitle = "% mt-genes",
                                                               cols = cols,
                                                               scale = input$selectScale,
                                                               showCells = input$violin_qc_showDots))
   
   toWebGL(ggplotly(p, source = "gene_qc", key = "pointNumber", tooltip = "text") %>% 
             layout(dragmode = "lasso", height = 500, width = 500,
                    legend = list(orientation = "h", x = 0, y = -0.2))) 
   }) #close boxplot_features
  })  #close try 
}) #close observeEvent
  
## remove cells modal -----------------------------------------------------
   observeEvent(input$remove_cells,{
     showModal(modalDialog(
       h1("Warning"),
       p("This action will permanently remove selected cells from the session . Are you sure
         you want to continue?"),

       footer = tagList(
         modalButton("Cancel"),
         actionButton("ok_remove_cells", "Remove cells"))
     )) #close modalDialog
   }) #close observe remove_cells
   
## remove poor quality cells ---------------------------------------
  observeEvent(input$ok_remove_cells,{
    removeModal()
    withProgress(message = 'Removing cells', value = 0, {
     
    try({
      
    if (input$selectScale == "Linear"){
      library_size <- overall_vars$metadata$library_size
      total_features <- overall_vars$metadata$total_features}
    
    if (input$selectScale == "Log2"){
      library_size <- log2(overall_vars$metadata$library_size)
      total_features <- log2(overall_vars$metadata$total_features)}
     
    if (input$selectScale == "Log10"){
      library_size <- log10(overall_vars$metadata$library_size)
      total_features <- log10(overall_vars$metadata$total_features)}
      
    keep_cells <- library_size >= input$selectLibSizes[1] &
                  library_size <= input$selectLibSizes[2] &
                  total_features >= input$selectNrFeatures[1] &
                  total_features <= input$selectNrFeatures[2] &
                  overall_vars$metadata$percentage_mt_genes <= input$selectMtCutoff
    
    print("Filtering cells:")
    print(table(keep_cells))
    
    # save parameters in report
    overall_vars$report$QC$filter_cells <- keep_cells
    overall_vars$report$QC$library_size_min <- input$selectLibSizes[1]
    overall_vars$report$QC$library_size_max <- input$selectLibSizes[2]
    overall_vars$report$QC$total_features_size_min <- input$selectNrFeatures[1]
    overall_vars$report$QC$total_features_size_max <- input$selectNrFeatures[2]
    overall_vars$report$QC$percentage_mt_genes <- input$selectMtCutoff
    overall_vars$report$QC$scale <- input$selectScale
    
    save_session(overall_vars$session_token,overall_vars, "report")
    
    incProgress(0.3, detail = "")
    
    # update count matrices
    
    overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
    
    for (mat in names(overall_vars$countMatrices)){
      overall_vars$countMatrices[[mat]] <- overall_vars$countMatrices[[mat]][,keep_cells]
      } # close for 
    
    save_session(overall_vars$session_token,overall_vars, "countMatrices")
    
    incProgress(0.3, detail = "")
    
    # update metadata
 
    overall_vars$metadata <- overall_vars$metadata[keep_cells,]

    save_session(overall_vars$session_token, overall_vars, "metadata")
    overall_vars$md5$metadata <- md5sum(
      paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
    
    incProgress(0.4, detail = "")
    
    # re-calculate gene QC metrics
    
    overall_vars$gene_qc[["metrics"]] <- calculate_gene_qc(overall_vars$countMatrices[["rawCountMatrix"]])
    
    cols <- c(as.vector(unlist(overall_vars$colours[[input$features_sample]])), extra_cols)
    
    overall_vars$gene_qc[["top_genes"]] <- 
      get_top_expressed_genes(overall_vars$countMatrices[["rawCountMatrix"]], overall_vars$metadata) 
    
    overall_vars$countMatrices <- NULL
    
    save_session(overall_vars$session_token, overall_vars, "gene_qc")
    
    # update dimensionality reduction - verify first if object exists
    files <- list.files(paste0(getwd(), "/tokens/", overall_vars$session_token))
    
    if ("dimred.rds" %in% files){
    
    overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
    
    if (!is.null(overall_vars$dimred$pca)){
      for(pca in names(overall_vars$dimred$pca)){
          overall_vars$dimred$pca[[pca]]$x <- overall_vars$dimred$pca[[pca]][keep_cells,]
      } #close for
    } # close if 
    
    if (!is.null(overall_vars$dimred$tsne)){
      for(tsne in names(overall_vars$dimred$tsne)){
        overall_vars$dimred$tsne[[tsne]] <- overall_vars$dimred$tsne[[tsne]][keep_cells,]
      } #close for
    } # close if 
    
    if (!is.null(overall_vars$dimred$umap)){
      for(umap in names(overall_vars$dimred$umap)){
        overall_vars$dimred$umap[[umap]] <- overall_vars$dimred$umap[[umap]][keep_cells,]
      } #close for
    } # close if 
    
    save_session(overall_vars$session_token, overall_vars, "dimred")
  } # close if 
    
    # update clusters 
    # verify first if object exists
    
    if ("clusters.rds" %in% files){
      overall_vars$clusters <- upload_session(overall_vars$session_token, "clusters")
      for (cluster_analysis in names(overall_vars$clusters)){
        overall_vars$clusters[[cluster_analysis]] <- overall_vars$clusters[[cluster_analysis]][keep_cells,]
      } #close for 
      save_session(overall_vars$session_token, overall_vars, "clusters")
    } #close if 
    
    # update scores 
    # verify first if object exists
    
    if ("scores.rds" %in% files){
      overall_vars$scores <- upload_session(overall_vars$session_token, "scores")
      for (score in names(overall_vars$scores)){
        overall_vars$scores[[score]] <- overall_vars$scores[[score]][keep_cells]
      } #close for 
      save_session(overall_vars$session_token, overall_vars, "scores")
    } #close if 
    
    #overall_vars$gene_qc <- calculate_gene_qc(overall_vars$countMatrices[["rawCountMatrix"]])
    #save_session(overall_vars$session_token, overall_vars, "gene_qc")
    
  showNotification(paste("Successfully removed",table(keep_cells)[1] ,"cells", sep =" "))
     }, silent = FALSE) #close Try
   }) # close progress bar
      
}) #close ok_remove_cells
  
## flag cells ---------------------------------------------
  observeEvent(input$flag_cells,{
    try({
      if (input$selectScale == "Linear"){
        library_size <- overall_vars$metadata$library_size
        total_features <- overall_vars$metadata$total_features}
      
      if (input$selectScale == "Log2"){
        library_size <- log2(overall_vars$metadata$library_size)
        total_features <- log2(overall_vars$metadata$total_features)}
      
      if (input$selectScale == "Log10"){
        library_size <- log10(overall_vars$metadata$library_size)
        total_features <- log10(overall_vars$metadata$total_features)}
      
      keep_cells <- library_size >= input$selectLibSizes[1] &
        library_size <= input$selectLibSizes[2] &
        total_features >= input$selectNrFeatures[1] &
        total_features <= input$selectNrFeatures[2] &
        overall_vars$metadata$percentage_mt_genes <= input$selectMtCutoff
      
      overall_vars$metadata$pass_qc <- keep_cells
      
    }) #close try
    
  }) #close observe flag_cells
  
##gene qc -----------------------------------------------  
### remove genes modal -----------------------------------------------------
   observeEvent(input$remove_features,{
     showModal(modalDialog(
       h1("Warning"),
       p("This action will permanently remove selected genes from the session . Are you sure
         you want to continue?"),
       
       footer = tagList(
         modalButton("Cancel"),
         actionButton("ok_remove_genes", "Remove genes"))
     )) #close modalDialog
   }) #close observe genes_cells
   
### remove lowly expressed genes ---------------------------------------
  observeEvent(input$ok_remove_genes,{
    removeModal()  
    withProgress(message = 'Removing features...', value = 0, {
      
      incProgress(0.6, detail = "")
      
      keep_genes <- overall_vars$gene_qc$metrics$total_exp > as.numeric(input$select_nr_counts) & 
                       overall_vars$gene_qc$metrics$expr_cells > as.numeric(input$select_nr_cells)
      
      print("Gene filter:")
      print(table(keep_genes)) 
      
      # annotate and save the analysis parameters in the report
      overall_vars$report$QC$filter_genes <- keep_genes
      
      overall_vars$report$QC$total_gene_expression <- as.numeric(input$select_exp_counts)
      
      overall_vars$report$QC$number_expressing_cells <- as.numeric(input$select_number_cells)
      
      save_session(overall_vars$session_token, overall_vars, "report")
    
      
      try({ 

      # remove the genes from the count matrices available in the session 
      overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")  
        
      for (mat in names(overall_vars$countMatrices)){
        overall_vars$countMatrices[[mat]] <- overall_vars$countMatrices[[mat]][keep_genes,]
        print(dim(overall_vars$countMatrices[[mat]]))
        
      } #close for 
        
      save_session(overall_vars$session_token, overall_vars, "countMatrices")
      saveRDS(rownames(overall_vars$countMatrices$rawCountMatrix),
              paste0(getwd(), "/tokens/", overall_vars$session_token, "/genes.rds"))
      
      incProgress(0.2, detail = "")
      
       # re-calculate the cell QC metrics 
       overall_vars$metadata$library_size <- NULL
       overall_vars$metadata$total_features <- NULL 
       overall_vars$metadata$percentage_mt_genes <- NULL
       
       mat <- overall_vars$countMatrices[["rawCountMatrix"]]  
     
       QCmetrics <- calculateQCmetrics(mat)
       
       print(QCmetrics[1:5,])
      
       overall_vars$metadata <- cbind(overall_vars$metadata, QCmetrics)
       
       save_session(overall_vars$session_token, overall_vars, "metadata")
       overall_vars$md5$metadata <- md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
        
       # subset the gene_qc object 
       overall_vars$gene_qc$metrics <- overall_vars$gene_qc$metrics[keep_genes,]
       
       save_session(overall_vars$session_token, overall_vars, "gene_qc")
       
       overall_vars$countMatrices <- list()

      }) #close try
      
      incProgress(0.2, detail = "")
      
   }) # close progress bar
    
   showNotification(paste("Successfully removed",table(keep_genes)[1] ,"genes", sep =" "))
   
  })

  
# CALCULATE DOUBLETS -----------------------------------------------
  
  observeEvent(input$calculate_doublets,{
    withProgress(message = 'Identifying doublets...', value = 0, {
      incProgress(0.4, detail = "")
      
      mat <- upload_session(overall_vars$session_token, "countMatrices") 
      mat <- as.matrix(mat[["rawCountMatrix"]])  
      
      doublets_result <-  calculate_doublets(mat = mat, 
                                             sample = overall_vars$metadata[,input$doublet_channel])
      
     #overall_vars$doublets <- doublets_result[[2]]  
     #save_session(overall_vars$session_token,overall_vars, c("doublets"))
     
     overall_vars$metadata$doublets_class <- "class"
     overall_vars$metadata$doublets_class[doublets_result$score >= input$selectDoubletCutoff] <- "doublet"
     overall_vars$metadata$doublets_class[doublets_result$score < input$selectDoubletCutoff] <- "singlet"
     overall_vars$metadata$doublets_score <- doublets_result$score
     
     rm(mat)
    
     }) # close progress bar
    showNotification("Calculated doublets.")
  }) #close calculate_doublets
   
# REMOVE DOUBLETS -----------------------------------------------
     
  observeEvent(input$remove_doublets,{
    withProgress(message = 'Removing doublets...', value = 0, {
      incProgress(0.4, detail = "")
      
      is_singlet <- overall_vars$metadata$doublets_class  
      is_singlet <- is_singlet == "singlet"
      
      print("Initial size of count matrix:")
      print(dim(overall_vars$countMatrices[["rawCountMatrix"]]))
      print(table(is_singlet))
      
      overall_vars$report$doublets$score <- overall_vars$metadata$doublets_score
      
      overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
      
      for (mat in names(overall_vars$countMatrices)){
        overall_vars$countMatrices[[mat]] <- overall_vars$countMatrices[[mat]][,is_singlet]
      }
      
      save_session(overall_vars$session_token, overall_vars, "countMatrices")
      
      overall_vars$countMatrices <- list()
      
      overall_vars$metadata <- overall_vars$metadata[is_singlet,]
      
      save_session(overall_vars$session_token, overall_vars, "metadata")
      
      # save parameters in report
      
      overall_vars$report$doublets$is_singlet <- is_singlet
      
      save_session(overall_vars$session_token, overall_vars, "report")
      
    }) # close progress bar
    
    showNotification("Successfully removed doublets.")
    
}) #close remove_doublets
  
  
# BARPLOTS -----------------------------------------------     
     observeEvent(input$go_summary,{
     try({
       new_var <- overall_vars$metadata[[input$summary_var[1]]]
       
       new_var <- str_sub(new_var, 1, 15)
       
       for (var in input$summary_var[-1]){
         new_var <- paste(new_var, overall_vars$metadata[[var]], sep ="_")
       }
       
       output$summary_table <- DT::renderDataTable({ 
         
         table <- as.data.frame(table(new_var)) 
         colnames(table) <- c("Condition", "#cells")
         
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
       
     }, silent = TRUE) #close try
  
       try({
       output$summary_barplot <- renderPlotly(
         { isolate(p <- get_summary_barplot(var1 = overall_vars$metadata[[input$summary_var[1]]], 
                                    var2 = overall_vars$metadata[[input$summary_var[2]]], 
                                    colors = c(overall_vars$colours[[input$summary_var[2]]], extra_cols), 
                                    label = input$summary_var[2]))
         
           toWebGL(ggplotly(p, source = "bp", key = "pointNumber", tooltip = "text") %>% 
                     layout(dragmode = "lasso", height = 800, 
                            legend = list(orientation = "h", x = 0, y = -0.8)))
         }) #close summary_barplot  
       }, silent = TRUE)#close try
       
     }) #close ObserveEvent
  
}#close server
  
  
 #shiny::runApp(launch.browser = FALSE)
 shinyApp(ui = ui, server = server)


