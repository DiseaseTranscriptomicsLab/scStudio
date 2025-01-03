##-------------------------------------##
##              UMAP TAB               ##
##-------------------------------------##

getUMAP <- function(mat, subset_row, min_dist, nneigh, token, session_obj, ID){
  
  umap <- scater::calculateUMAP(x = mat,
                                ncomponents = 10,
                                subset_row = subset_row,
                                min_dist = min_dist,
                                n_neighbors = nneigh)
  
  session_obj$umap[[ID]] <- as.data.frame(umap)
  saveRDS(session_obj, paste0(getwd(),"/tokens/", token, "/dimred.rds"))
  
} 
 
 
 


