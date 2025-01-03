##-------------------------------------##
##             TSNE TAB                ##
##-------------------------------------##

 getTSNE <- function(mat, subset_row, perplexity, token, session_obj, ID){
 
  tsne <- scater::calculateTSNE(x = mat,
                  ncomponents = 3, 
                  subset_row = subset_row,
                  perplexity)
  
  session_obj$tsne[[ID]] <- as.data.frame(tsne)
  saveRDS(session_obj, 
          paste0(getwd(),"/tokens/", token, "/dimred.rds"))
  

 }

 
 
 


