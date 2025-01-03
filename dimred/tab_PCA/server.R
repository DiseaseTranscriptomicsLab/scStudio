##-------------------------------------##
##              PCA TAB                ##
##-------------------------------------##

getPCA <- function(mat, ncomponents, subset_row, do_scaling, token, session_obj, ID){
  mat <- as.matrix(mat)
  if(subset_row == "All"){subset_row <- rownames(mat)}
  
      pca <- BiocSingular::runPCA(x = t(mat[subset_row,]), 
                                  rank = ncomponents,
                                  center = TRUE,
                                  scale = as.logical(do_scaling))
      pca$x <- as.data.frame(pca$x)
      session_obj$pca[[ID]] <- pca
      saveRDS(session_obj, 
              paste0(getwd(),"/tokens/", token, "/dimred.rds")
              )
      
  
}



plot_dimRed <- function(coord, x, y, metadata, sel_var, title, 
                        pointSize, alpha = 1, cols,
                        density_lines)
  {
    
  text <- paste("Condition: ", metadata[,sel_var],
                paste0("\nDim",x,": "), coord[,x],
                paste0("\nDim",y,": "), coord[,y])
  
  p <- ggplot(data = coord, 
              aes(x = coord[,x], 
                  y = coord[,y],
                  color = metadata[,sel_var],
                  text = text,
                  key = 1:nrow(coord))) +  
    geom_point(size = pointSize, alpha = alpha) +   
    scale_color_manual(values = cols) + 
    ggtitle(title) + 
    labs(color = sel_var) + 
    xlab(paste("Dim",x, sep ="_")) +
    ylab(paste("Dim",y, sep ="_")) +
    
    theme(text = element_text(size=15), 
          legend.position="bottom",
          legend.title = element_text(size=14), 
          legend.text = element_text(size=10),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  if(density_lines){p <- p + geom_density_2d(colour = "black")}
  
  return(p)
}


plot_scree <- function(pca){
  
  var_explained <- pca$sdev**2
  total_variance_explained <- sum(pca$sdev**2)
  proportion_var_explained <- (var_explained/total_variance_explained)*100 
  
  df <- data.frame(proportion_var_explained = proportion_var_explained, 
                   components = 1:length(proportion_var_explained))
  
  p <- ggplot(data = df, 
              aes(x = components, 
                  y = proportion_var_explained)) +  
    geom_point(size = 2.5) +   
    geom_line() +
    ggtitle("Scree plot: PCA") + 
    xlab("PC") +
    ylab("Variance explained (%)") +
    theme(text = element_text(size=15), 
          legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
 

 


