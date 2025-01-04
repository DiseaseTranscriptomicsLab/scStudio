##-------------------------------------##
##            VIOLINS TAB              ##
##-------------------------------------##

plot_grid_violins <- function(mat, genes, metadata, var, subvar, cols, pt.size,
                              gene_list){
  #https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
  set.seed(1234)
  #text <- paste("Sample: ", df[["orig.ident"]], 
  #              "Dataset: ", df[["dataset"]])
  plots <- list()
  

  for (gene in genes){

    if (gene %in% gene_list){
 
      ylabel <- "Log2(counts + 1)"
    } else{

      ylabel <- "Score"}
    

    
    plots[[gene]] <- local({
      gene <- gene
  
  if (subvar != "none"){        
  df <- data.frame(plot_gene = mat[gene,], 
                   condition = metadata[[var]],
                   subcondition = metadata[[subvar]])
  p <- ggplot(df, aes(x = df[["condition"]], 
                      y = df[["plot_gene"]], 
                      color = df[["condition"]],
                      fill = df[["subcondition"]])) + 
    geom_violin(size = 0.2, colour = "black", scale = "width")  +
    scale_fill_manual(values = cols) + 
    ggtitle(gene) +
    xlab("") + 
    labs(fill = subvar) + 
    ylab(ylabel) +
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle=45, hjust=1), 
          text = element_text(size=17), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  if (pt.size == 0){print(p)}
  else {
    p <- p +  geom_point(position = position_jitterdodge(), size = pt.size,
                         alpha = 0.5, dodge.width = 0.5, 
                         seed = 1234,
                         colour = "black")}
      #geom_jitter(shape= 20, position=position_jitter(0.1),  
              #            alpha = 0.5, 
              #            size = pt.size , colour = "black")}
  print(p)
  }
      
  else {
    df <- data.frame(plot_gene = mat[gene,], 
                     condition = metadata[[var]])
    
    p <- ggplot(df, aes(x = df[["condition"]], 
                        y = df[["plot_gene"]], 
                        color = df[["condition"]],
                        fill = df[["condition"]])) + 
      geom_violin(size = 0.2, colour = "black", scale = "width")  +
      stat_summary(fun.y=median, geom="point", size=2, color="red") + 
      stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
      scale_fill_manual(values = cols) + 
      ggtitle(gene) +
      xlab("") + 
      ylab(ylabel) +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=45, hjust=1), 
            text = element_text(size=17), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    if (pt.size == 0){print(p)}
    else {
      p <- p +  geom_jitter(shape= 20, position=position_jitter(0.3),  
                            alpha = 0.5, 
                            size = pt.size , colour = "black")}
    print(p)
    }    
  
  
    }) #close local
  } #close loop
  
  if (length(genes) %in% c(2,3,4)){number_cols <- 2}
  else if (length(genes) == 1 ){number_cols <- 1}
  else if (length(genes) %in% c(5,6,7,8,9) ){number_cols <- 3}
  else if (length(genes) > 9 ){number_cols <- 5}
  
  return(grid.arrange(grobs = plots, ncol = number_cols , top=""))
  #return(plots)
}


