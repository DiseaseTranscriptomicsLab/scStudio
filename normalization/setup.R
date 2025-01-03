##-------------------------------------##
#####          SETUP OPTIONS        #####
##-------------------------------------##
options(warn = -1)

set.seed(08071993)

options(shiny.maxRequestSize=900000*1024^2) 

options(spinner.color="#E7F5F6", 
        spinner.color.background="#ffffff", 
        spinner.size=0.5)
#https://projects.lukehaas.me/css-loaders/

extra_cols <- c("#6fccc1", "#cf1fa7", "#37ee7c", "#da0322", "#ffd94a", "#207ede", "#fea111", "#b09de9",
                "#b86fbb", "#150de2", "#ab3f2f", "#d49b6c", "#773aab", "#759474", "#e4b2b8", "#f88a40",
                "#74b5cf", "#e6208e", "#0a3535", "#5b37c3", "#cc4025", "#d2dc7c", "#344d76", "#ba8b07",
                "#2e10b6", "#4d0cb0", "#e1668a", "#47a58b", "#d7734a", "#ff58c7", "#edceef", "#a21d59",
                "#8c8da5", "#ff3e5e", "#688d4b", "#214fca", "#48a8f5", "#752758", "#3b4a4a", "#11674b",
                "#9a1bd9", "#16146d", "#7d7277", "#2d8cce", "#2460b8", "#0f9a8f", "#b11e90", "#54ed4f",
                "#987217", "#980119", "#cf0606", "#dd182c", "#88a1e4", "#3b89ba", "#12defd", "#feba3e",
                "#e31fa5", "#7b3537", "#098a95", "#1f1716", "#df6c94", "#9965ee", "#b438f9", "#0140f5",
                "#e6c4a8", "#d94f7c", "#0dae56", "#b86d1d", "#0577ad", "#464551", "#6b0959", "#ccf705",
                "#12271d")

upload_session <- function(token, type){
  dir_path <- paste0("./tokens/", token, "/", type, ".rds")
  print(dir_path)
  obj <- readRDS(dir_path)
  return (obj)
}

save_session <- function(token, overall_vars, objects){
  available_sessions <- list.files(paste0(getwd(), "/tokens/"))
  dir_path <- paste0("./tokens/", token)
  if (!(token %in% available_sessions )){
    
    dir.create(path = dir_path)
  } # close if
  
  for (name in objects){
    saveRDS(overall_vars[[name]], paste0(dir_path,"/",name, ".rds"))}
} # close save sesion

identify_discrete <- function(df){
  is_discrete <- vector()
  
  for (col in colnames(df)){
    if (length(unique(df[[col]])) > 75){
      is_discrete <- c(is_discrete, FALSE)
    }
    else {is_discrete <- c(is_discrete, TRUE)}
  }
  return(df[,is_discrete])
}

identify_continuous <- function(df){
  is_continuous <- c()
  
  for (col in colnames(df)){

    
    if (is.character(df[[col]])){

      is_continuous <- c(is_continuous, FALSE)
    }
    else {
      is_continuous<- c(is_continuous, TRUE)
  }
  }

  return(df[,is_continuous])
}


verify_null <- function(string, metadata){
  if (string == "NULL"){ return(NULL)}
  else return(metadata[[string]])
}
