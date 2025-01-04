##-------------------------------------##
##           FEATPLOTS TAB             ##
##-------------------------------------##


plot_featplot <- function(mat, genes, dimred, type, x, y, density_lines)
{ 
  plots <- list()
  if (type == "pca"){dimred <- dimred$x}
  
  dimred <- dimred[,c(x,y)]
  names(dimred) <- c("V1", "V2")
  list_of_dimreds <- list()
  for (gene in genes){
    print(gene)
  mat2 <- mat[gene,]
  dimred$var_col <- mat2
  list_of_dimreds[[gene]] <- dimred[order(dimred$var_col, decreasing = FALSE),]
  print(names(list_of_dimreds))
  

  p <- ggplot(data = list_of_dimreds[[gene]], 
              aes(x = V1, 
                  y = V2,
                  color = var_col,
                  key = 1:nrow(list_of_dimreds[[gene]]))) +  
    geom_point() +   
    scale_colour_gradient(low = "gray88", high = "red" ) + 
    labs(color = gene) +
    xlab(paste("Dim",x, sep ="_")) +
    ylab(paste("Dim",y, sep ="_")) +
    theme(text = element_text(size=10), 
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
    
  if (density_lines){p <- p + geom_density_2d()}
  
  plots[[gene]] <- p

  } 
 if (length(genes) %in% c(1,2,3,4)){number_cols <- 2}
 else if (length(genes) == 1 ){number_cols <- 1}
 else if (length(genes) %in% c(5,6,7,8,9) ){number_cols <- 3}
 else if (length(genes) > 9 ){number_cols <- 5}
 #return(plots)
 #return(grid.arrange(grobs = plots, ncol = number_cols , top=""))
  return(ggarrange(plotlist = plots))
 
}


