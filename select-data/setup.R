##-------------------------------------##
#####          SETUP OPTIONS        #####
##-------------------------------------##

options(timeout = max(1000, getOption("timeout")))

options(download.file.method.GEOquery = "auto")

options(warn = -1)

set.seed(08071993)

options(shiny.maxRequestSize=900000*1024^2) 

options(spinner.color="#E7F5F6", 
        spinner.color.background="#ffffff", 
        spinner.size=0.5)
#https://projects.lukehaas.me/css-loaders/

get_random_token <- function(){
  set.seed(as.numeric(Sys.time()))
  token <- stri_rand_strings(1, 10, pattern = "[a-z0-9]")
  print(token)
  return(token)
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

remove_objects <- function(token){
  files <- list.files(paste0("tokens/",token))
  for (file in files){
    if (!(file %in% c("countMatrices.rds", "metadata.rds"))){
      unlink(paste0(getwd(),"/tokens/", token,"/", file))
    } #close for
  } #close if 
} #close remove_objects

upload_session <- function(token, type){
  dir_path <- paste0("./tokens/", token, "/", type, ".rds")
  print(dir_path)
  obj <- readRDS(dir_path)
  return (obj)
}


