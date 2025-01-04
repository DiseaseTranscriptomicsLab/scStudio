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
library(Matrix)
library(matrixStats)
library(magrittr)
library(plotly)
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

source("tab_COLOURS/UI.R")
source("tab_COLOURS/server.R")

source("tab_METADATA/UI.R")
source("tab_METADATA/server.R")

source("tab_SUBSET/UI.R")
source("tab_SUBSET/server.R")

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
  navbarPage("scStudio: Settings",
          tabPanel("Colours",
                     tab_COLOURS),
          tabPanel("Metadata", 
                   tab_METADATA),
          tabPanel("Subset",
                   tab_SUBSET)

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
       h1("Welcome to scStudio - Settings"),
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
             
             updateSelectInput(session, "select_discrete_col_var",
                               choices = names(identify_discrete(overall_vars$metadata)))
             
             updateSelectInput(session, "select_var_annotation",
                               choices = names(identify_discrete(overall_vars$metadata)),
                               selected = "")
             
             overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
             overall_vars$md5$mat_names <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
             
             
             overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
             overall_vars$md5$dimred <- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
            
             incProgress(0.6, detail = "Finishing...")
             
             if ("colours.rds" %in% files){
               overall_vars$colours <- upload_session(overall_vars$session_token, "colours")
             } #close if 
             
             else{overall_vars$colours <- list()}
             
             if ("clusters.rds" %in% files){
               overall_vars$clusters<- upload_session(overall_vars$session_token, "clusters")
               
               updateSelectInput(session, "select_clustering_analysis",
                                 choices = names(overall_vars$clusters))
             } #close if 
             
             else{overall_vars$clusters <- list()}
             
             overall_vars$md5$clusters<- md5sum(
               paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
             
             if ("report.rds" %in% files){
               overall_vars$report<- upload_session(overall_vars$session_token, "report")
             } #close if
             else {overall_vars$report <- list()}
             
             
  
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
           
           updateSelectInput(session, "select_var_annotation",
                             choices = names(identify_discrete(overall_vars$metadata)),
                             selected = "")
           
           updateSelectInput(session, "select_discrete_col_var",
                             choices = names(identify_discrete(overall_vars$metadata)))
           
           overall_vars$md5$mat_names <- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/mat_names.rds"))
           
           
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
           
           if ("clusters.rds" %in% files){
             overall_vars$clusters<- upload_session(overall_vars$session_token, "clusters")
             
             updateSelectInput(session, "select_clustering_analysis",choices = names(overall_vars$clusters))
           } #close if 
           
           else{overall_vars$clusters <- list()}
           
           overall_vars$md5$clusters<- md5sum(
             paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
           
           
           if ("report.rds" %in% files){
             overall_vars$report <- upload_session(overall_vars$session_token, "report")
           } #close if 
          
            else{overall_vars$report <- list()}
           
        
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
      
      check_md5sum_dimred <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
      
      check_md5sum_clusters <- md5sum(
        paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
      
      # update objs if new data is available from other scStudio tools
      
      if (check_md5sum_metadata != overall_vars$md5$metadata){
        overall_vars$metadata <- upload_session(overall_vars$session_token, "metadata")
        
        overall_vars$md5$metadata <- check_md5sum_metadata
        
        updateSelectInput(session, "select_discrete_col_var",
                          choices = names(identify_discrete(overall_vars$metadata)))
        
        updateSelectInput(session, "select_var_annotation",
                          choices = names(identify_discrete(overall_vars$metadata)),
                          selected = "")
        
       } # close if
      
      if (check_md5sum_mat_names != overall_vars$md5$mat_names){
        overall_vars$mat_names <- upload_session(overall_vars$session_token, "mat_names")
        overall_vars$md5$mat_names <- check_md5sum_mat_names
      } # close if
      
      if (check_md5sum_dimred != overall_vars$md5$dimred){
        overall_vars$dimred <- upload_session(overall_vars$session_token, "dimred")
        overall_vars$md5$dimred <- check_md5sum_dimred
      } # close if
      
      if (check_md5sum_clusters != overall_vars$md5$clusters){
        overall_vars$clusters<- upload_session(overall_vars$session_token, "clusters")
        overall_vars$md5$clusters <- check_md5sum_clusters
        
        updateSelectInput(session, "select_clustering_analysis",choices = names(overall_vars$clusters))
        
      } # close if
      
 
      
    }) #close try
  
}) #close observe
  

   
## update inputs ----------------------------------------------------   
observe({
 try({ 
   plots <- c(names(overall_vars$dimred$tsne),
              names(overall_vars$dimred$umap))
   
   updateSelectInput(session, "choose_plot_subset",
                     choices = plots,
                     selected = plots[1])
   
  choices <- unique(overall_vars$metadata[[input$select_discrete_col_var]])
  updateSelectInput(session, "select_discrete_col_group",choices = choices)
  
  
  choices <- unique(overall_vars$metadata[[input$select_var_annotation]])
  updateSelectInput(session, "select_group_annotation",choices = choices)
  
  choices <- names(overall_vars$clusters[[input$select_clustering_analysis]])
  updateSelectInput(session, "select_clustering_annotation",choices = choices)
  
  choices <- colnames(identify_discrete(overall_vars$metadata))
  updateSelectInput(session, "subset_var",choices = choices, selected = choices[1])

  df <- data.frame(colors = unlist(overall_vars$colours[[input$select_discrete_col_var]]),
                   group = names(overall_vars$colours[[input$select_discrete_col_var]]),
                   add = rep(1, length(names(overall_vars$colours[[input$select_discrete_col_var]]))))

  
  try(df <- df[order(df$group),])

  output$colors_plot <- renderPlotly({
    p <- make_show_color_plot(df)
    
    toWebGL(ggplotly(p, source = "colors", key = "pointNumber", tooltip = "text") %>% 
              layout(dragmode = "lasso", height = 300, 
                     legend = list(orientation = "h", x = 0, y = -0.2))) 
  })
  
 }, silent = FALSE) #close try 
}) #close observe
 
   
# RUN COLOURS -----------------------------------------------  
observeEvent(input$save_discrete_col, {
     
     overall_vars$colours[[input$select_discrete_col_var]][[input$select_discrete_col_group]] <- input$select_discrete_col 
     print(overall_vars$colours) 
     save_session(overall_vars$session_token, overall_vars, "colours")
     showNotification("Successfully saved colour.")
     
   }) #close save_discrete_col   

# RUN METADATA -----------------------------------------------
   
observeEvent(input$save_new_annotation, {
  if (input$change_var_name != ""){
    
    if (! (input$change_var_name %in% names(overall_vars$metadata))){
     
      names(overall_vars$metadata) <- gsub(input$select_var_annotation, input$change_var_name, names(overall_vars$metadata))
      
      names(overall_vars$colours) <- replace_label(
        names(overall_vars$colours), 
        input$select_var_annotation, input$change_var_name)
      
      save_session(overall_vars$session_token, overall_vars, "colours")   
      
      if (input$select_new_annotation != ""){

           new_var <- replace_label(var = overall_vars$metadata[[input$change_var_name]],
                                    old_label = input$select_group_annotation,
                                    new_label = input$select_new_annotation)
   
           overall_vars$metadata[[input$change_var_name]] <- new_var
           
         }
      
      updateSelectInput(session, "select_discrete_col_var",
                        choices = names(identify_discrete(overall_vars$metadata)),
                        selected = NULL)
      
      updateSelectInput(session, "select_var_annotation",
                        choices = names(identify_discrete(overall_vars$metadata)),
                        selected = input$change_var_name)
      
      updateTextInput(session, "change_var_name", value  = "")
      
      save_session(overall_vars$session_token, overall_vars, "metadata")
      
      
      
      
         
         }
       
       
       else{showNotification("Cannot have duplicated names in metadata.")}
  }
  
  else{

     if (input$select_new_annotation != ""){

       new_var <- replace_label(var = overall_vars$metadata[[input$select_var_annotation]],
                                old_label = input$select_group_annotation,
                                new_label = input$select_new_annotation)

       overall_vars$metadata[[input$select_var_annotation]] <- new_var
       

       updateSelectInput(session, "select_var_annotation",
                         choices = names(identify_discrete(overall_vars$metadata)),
                         selected = input$change_var_name)
       
       save_session(overall_vars$session_token, overall_vars, "metadata")
       
       
       names(overall_vars$colours[[input$select_var_annotation]]) <- replace_label(
         names(overall_vars$colours[[input$select_var_annotation]]), 
         input$select_group_annotation, input$select_new_annotation)
       
       save_session(overall_vars$session_token, overall_vars, "colours")
     }
  } #close else
   })
   
#*****
   observeEvent(input$add_clustering_metadata, {
     var <- overall_vars$clusters[[input$select_clustering_analysis]][[input$select_clustering_annotation]]
     overall_vars$metadata[[paste(input$select_clustering_analysis, 
                                  input$select_clustering_annotation, sep = "_")]] <- var
     
     showNotification("Clustering analysis added to metadata.")
     
     #update select options    
     updateSelectInput(session, "select_var_annotation",choices = names(identify_discrete(overall_vars$metadata)),
                       selected = "")
     updateSelectInput(session, "dea_condition",choices = names(identify_discrete(overall_vars$metadata)))
   })
   #*****  
   observeEvent(input$subset_var, {
     choices <- overall_vars$metadata[,input$subset_var]
     updateSelectInput(session, "subset_group",choices = choices)
   })   
   
   
   
# RUN subset dataset -------------------------------------------------------------------    
   observeEvent(input$subset_group, {
     isolate({
       output$subset_tsne <- renderPlotly({ 
         
         print("Making tsne with subset group:")
         print(input$subset_group)
         
         coordinates <- get_groups_dea(input$subset_group, input$subset_grouping, overall_vars$metadata)
         subset_by <- rep(FALSE, nrow(overall_vars$metadata))
         subset_by[coordinates] <- TRUE
         
         print(table(subset_by))
         
         
         if(input$choose_plot_subset %in%  names(overall_vars$dimred$tsne)){
         p <- plot_dimRed(coord = overall_vars$dimred$tsne[[input$choose_plot_subset]], 
                          x = 1, 
                          y = 2, 
                          metadata = data.frame(subset_by = subset_by), 
                          sel_var = "subset_by", 
                          title = paste0("Number of selected cells: ", as.numeric(table(subset_by)[2])),
                          pointSize = 0.5,
                          alpha = 1, 
                          cols = c("grey", "green"))
         
         
         toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                   layout(dragmode = "lasso", height = 800, 
                          legend = list(orientation = "h", x = 0, y = -0.2)))
         } #close if
         
         else{
           p <- plot_dimRed(coord = overall_vars$dimred$umap[[input$choose_plot_subset]], 
                            x = 1, y = 2, 
                            metadata = data.frame(subset_by = subset_by), 
                            sel_var = "subset_by", 
                            title = paste0("Number of selected cells: ", as.numeric(table(subset_by)[2])),
                            pointSize = 0.5,
                            alpha = 1, 
                            cols = c("grey", "green"))
           
           toWebGL(ggplotly(p, source = "main", key = "pointNumber", tooltip = "text") %>% 
                     layout(dragmode = "lasso", height = 800, 
                            legend = list(orientation = "h", x = 0, y = -0.2)))
         } #close else 
       })
     })
   })
   
  observeEvent(input$subset_dataset, {
     
     showModal(modalDialog(
       title = "Important message",
       tags$p("Subsetting your dataset will save your current session and open a new one. 
        Are you sure you want to proceed?"),
       footer = tagList(
         modalButton("Cancel"),
         actionButton("ok_subset", "OK")
       )#close footer
     )#close modalDialog
     ) #close showModal
   })
   
   observeEvent(input$ok_subset, {
     removeModal()
     withProgress(message = 'Subsetting in progress', value = 0, detail="0%", {  
    
       # Initiate new session with subset 
       overall_vars$countMatrices <- upload_session(overall_vars$session_token, "countMatrices")
       
       overall_vars$session_token <- get_random_token()
       print("Defined new token session.")
       print(overall_vars$session_token)
       output$session_id <- renderText({ paste0("Session ID: ", overall_vars$session_token)})
       
       coordinates <- get_groups_dea(input$subset_group, input$subset_grouping, overall_vars$metadata)
       subset_by <- rep(FALSE, nrow(overall_vars$metadata))
       subset_by[coordinates] <- TRUE
       
       
       
       for (mat in names(overall_vars$countMatrices)){
         overall_vars$countMatrices[[mat]] <- overall_vars$countMatrices[[mat]][,subset_by]}
       print("Subsetting matrices done.")
       
       incProgress(0.3, detail = "30%")
       
       for (reduction in names(overall_vars$dimred$tsne)){
         overall_vars$dimred$tsne[[reduction]] <- overall_vars$dimred$tsne[[reduction]][subset_by,]}
       print("Subsetting tsnes done.")
       
       for (reduction in names(overall_vars$dimred$umap)){
         overall_vars$dimred$umap[[reduction]] <- overall_vars$dimred$umap[[reduction]][subset_by,]}
       print("Subsetting umaps done.")
       
       for (reduction in names(overall_vars$dimred$pca)){
         overall_vars$dimred$pca[[reduction]]$x <- overall_vars$dimred$pca[[reduction]]$x[subset_by,]}
       print("Subsetting pcas done.")
       
       for (cluster in names(overall_vars$clusters)){
         overall_vars$clusters[[cluster]] <- overall_vars$clusters[[cluster]][subset_by,]}
       print("Subsetting clusters done.")
       
       overall_vars$metadata <- overall_vars$metadata[subset_by,]
       print("Subsetting metadata done.")
       
       incProgress(0.4, detail = "40%")
       
       overall_vars$gene_qc <- data.frame()
       overall_vars$dea <- list()
       overall_vars$gsea <- list()
       overall_vars$scores <- list()
       print("Resetting remaining objects done.")
       
       # Save new session objects
       save_session(overall_vars$session_token, overall_vars, "countMatrices")
       
       incProgress(0.5, detail = "50%")
       
       save_session(overall_vars$session_token, overall_vars, "metadata")
       save_session(overall_vars$session_token, overall_vars, "mat_names")
       save_session(overall_vars$session_token, overall_vars, "gene_qc")
       save_session(overall_vars$session_token, overall_vars, "hvgs")
       save_session(overall_vars$session_token, overall_vars, "dimred")
       save_session(overall_vars$session_token, overall_vars, "clusters")
       save_session(overall_vars$session_token, overall_vars, "mkgs")
       save_session(overall_vars$session_token, overall_vars, "dea")
       save_session(overall_vars$session_token, overall_vars, "gsea")
       save_session(overall_vars$session_token, overall_vars, "colours")
       save_session(overall_vars$session_token, overall_vars, "report")
       save_session(overall_vars$session_token, overall_vars, "scores")
       save_session(overall_vars$session_token, overall_vars, "genes")
       
       incProgress(0.6, detail = "60%")
       
       # Verify md5sum
       print("Verifying md5sum")
       overall_vars$md5$dimred <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/dimred.rds"))
       print("dimred done")
       overall_vars$md5$clusters <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/clusters.rds"))
       print("clusters done")
       overall_vars$md5$metadata <-  md5sum(
         paste0(getwd(),"/tokens/", overall_vars$session_token, "/metadata.rds"))
       print("metadata done")
       
       observeEvent(input$subset_var, {
         choices <- overall_vars$metadata[,input$subset_var]
         updateSelectInput(session, "subset_group",choices = choices)
       })
       
      
       showModal(modalDialog(
         title = "Important message",
         tags$p(paste0("Subsetting done. New session token: ", overall_vars$session_token)),
         footer = tagList(
           modalButton("OK"),

         )#close footer
       )#close modalDialog
       ) #close showModal

       
     })  #%>%  tagAppendAttributes(class = 'test') # close progress bar
   })
   
   
   

}#close server
  
  
 
 shinyApp(ui = ui, server = server)


