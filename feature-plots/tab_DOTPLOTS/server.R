##-------------------------------------##
##            DOTPLOTS TAB             ##
##-------------------------------------##

get_average_scores <- function(score, clusters, meta){
  average_score <- c()
  for (cluster in clusters){
    average_score <- c(average_score, mean(score[meta == cluster]))
  }
  return(average_score)
}

get_percent_cells <- function(score, clusters, meta){
  total_over <- c()
  for (cluster in clusters){
    score_of_cluster <-  score[meta == cluster]
    total_cells <- length(score_of_cluster )
    over <- sum(score_of_cluster > 0) 
    over_percentage <- (over/total_cells)*100
    total_over <- c(total_over, over_percentage)
  }
  return(total_over)
}

plot_dotplot <- function(mat, genes, var){
  total_average_exp <- c()
  total_percent_cells <- c()
  groups <- unique(var)
  mat <- mat[genes,]
  mat <- (2**mat) - 1
  
  for (gene in unique(genes)){
    
    gene_average_exp <- get_average_scores(score = mat[gene,],
                                           clusters = groups,
                                           meta = var)
    
    total_average_exp <- c(total_average_exp, gene_average_exp)
    
    gene_percent_cells <- get_percent_cells(score = mat[gene,], 
                                            clusters = groups, 
                                            meta = var)
    
    total_percent_cells <- c(total_percent_cells, gene_percent_cells)
  } # close for
  
  df <- data.frame(group = rep(groups, length(unique(genes))),
                   total_average_exp = total_average_exp,
                   total_percent_cells = total_percent_cells,
                   gene = rep(genes, each = length(groups))
                   )

colnames(df) <- c("Group", "Average expression (counts)", "% cells with expression > 0", "Gene") 

ggplot(df)+
  geom_point(mapping = aes(x = `Group`, 
                           y = `Gene`, 
                           color = `Average expression (counts)`, 
                           size = `% cells with expression > 0`)) +  
  scale_color_gradientn(colors = c("blue", "#F9F487", "red", "#A00018")) +

  
  theme_classic() + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +

  theme(text = element_text(size = 15),
        axis.text.x = element_text(size = 17,colour="black"),
        axis.text.y = element_text(size = 16,colour="black")) + ylab("") + xlab("") +
  ggtitle("")

}






#plot_dotplot <- function(mat, genes, gsets, scores, metadata, var, subvar, scale_TF, cols, flip){
#  print("Building dotplot...")
#  # AverageExpression() assumes "data" in log(counts+1)
#  norm_mat <- (2**mat)-1
#  ln_mat <- log(norm_mat+1)
#  genes <- gsub("_", "-", genes)
#  srt <- Seurat::CreateSeuratObject(counts = norm_mat, min.cells = 0, min.features = 0)
#  print("Created seurat object successfully")
#  srt$group <- metadata[[var]]
#  for (gs in gsets){
#    srt[[gs]] <- scores[gs]
#  }
#
#  if(subvar == "none"){
#    
#    p <- Seurat::DotPlot(object = srt,
#                         features = rev(c(genes, gsets)),
#                         cols = cols, 
#                         scale = scale_TF,
#                         group.by = "group") + 
#      xlab("") +ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#  
#    }
#  
#  else{
#    group <- paste(metadata[[var]], metadata[[subvar]], sep ="_")
#    srt$group <- group
#    print("Here2")
#    p <- Seurat::DotPlot(object = srt,
#                         features = rev(c(genes, gsets)),
#                         scale = scale_TF,
#                         cols = cols,
#                         group.by = "group",
##    ) +  xlab("") +ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) }
#    
# 
#  print("Finished")  
#  if(flip){return(p)} 
#  
#  else{return(p + coord_flip())}
#    #
#
#}




