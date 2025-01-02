##-------------------------------------##
##               QC - GENES            ##
##-------------------------------------##

calculate_gene_qc <- function(countMatrix){
  print("Calculating GENE QC")
  average_exp <- rowMeans(countMatrix) 
  total_exp <- rowSums(countMatrix)
  
  expr_cells <- countMatrix > 0 
  expr_cells <- rowSums(expr_cells)
  return(data.frame(average_exp = average_exp, expr_cells = expr_cells, total_exp = total_exp))
}


plot_scatter_features <- function(filter, filterx, filtery, gene_qc){
  
  text <- paste0("Gene: ",rownames(gene_qc), "\n",
                 "#cells: ",gene_qc$expr_cells)
  
  ggplot(gene_qc, aes(x = log10(total_exp + 1), y = log10(expr_cells + 1))) +
    #stat_density2d(aes(alpha = as.factor(..density..)), geom = "raster", contour = FALSE, 
    #               show.legend = FALSE) +  
    
  
    geom_point(aes(text = text), size = 0.2, alpha = 0.5 ) +
        theme(text = element_text(size=15), 
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_density_2d()+ 
    ylab("log10(number of expressing cells + 1)") + xlab("log10(total counts + 1)") +
    
    geom_vline(xintercept = log10(filterx+1), linetype="dotted",
             color = "red", size=0.5)+
    
    geom_hline(yintercept = log10(filtery+1), linetype="dotted",
               color = "red", size=0.5) +
    
    ggtitle(paste("Features to be removed:", sum(filter), sep = " ")) + 
    
    labs(color = "log10(average counts + 1)")
}


plot_tufte_features <- function(mat, filter, gene_qc, vars, var){
  mat <- mat[filter,]
  gene_qc <- gene_qc[filter,]
  df <- as.data.frame(mat)
  df <- df[order(gene_qc$average_exp),]
  
  df$genes <- factor(rownames(df), levels = rownames(df))
  melted <- reshape2::melt(df, id.vars = c("genes"))
  indx <- match(melted$variable, vars$barcodes)
  melted$var_color <- vars[[var]][indx]
  
  print(melted[1:10,])
  
  p <- ggplot(data = melted, aes(x = genes, y = value)) +

    geom_tufteboxplot() + 
    xlab("Features") +
    ylab("Counts") +
    ggtitle(paste("Features to be removed:", dim(df)[1], sep = " ")) +
    theme(text = element_text(size=15), 
       
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none") + 
    facet_grid(~var_color) +
    stat_summary(fun.y=median, geom="point", size=0.1, color="red")
  
  print("Finished tufte plot")
  
  return(p)
}

#plotTopGenes <- function(mat, metadata, sample, cols){
#  print("Plotting top most expressed genes...")
#  
#  sce <- SingleCellExperiment(assays=list(counts=as.matrix(mat)))
#  sce$sample <- metadata[[sample]]
#  
#  p <- scater::plotHighestExprs(sce, n = 30, colour_cells_by = "sample") + 
#       geom_point(size = 5, shape = 124, alpha = 1) +
#       scale_color_manual(values = cols) +
#    
#       stat_summary(fun.y=median, geom="point", size=2, color="red") + 
#       stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
#
#       theme(text = element_text(size=20),
#          legend.position="bottom") + 
#    
#       ggtitle("Top most expressed genes") +
#    
#       ylab("") + 
#       xlab("% of counts in each cell")  + 
#       labs(colour = "") +
#    
#       guides(color = guide_legend(override.aes = list(
#         linetype = 0, size = 3, shape = 16)))
#  
#    return(p)
#}

get_top_expressed_genes <- function(mat, metadata){
  total_feats <- rowSums(mat)
  total_feats <- total_feats[order(total_feats, decreasing = TRUE)][1:50]
  names_feats <- names(total_feats)
  mat <- mat[names_feats,]
  ratio <- sweep(mat, 2, metadata$library_size, "/") * 100
  # Convert the matrix to a data frame with gene names as a column
  percentage_df <- as.data.frame(ratio)
  percentage_df$Gene <- factor(rownames(ratio), levels = rev(rownames(ratio)))  # Add gene names as a column
  percentage_long <- melt(
    percentage_df, id.vars = c("Gene"), variable.name = "Cell", value.name = "Percentage")
  
  percentage_long$Percentage <- round(percentage_long$Percentage,2)
  
  
  
  return(percentage_long)
}

plot_top_genes <- function(df, log_scale, nr_genes){
  x_legend <- "% of counts in each cell"
  if (log_scale == "Log2"){
    df$Percentage <- log2(df$Percentage + 0.001)
    x_legend <- "Log2(% counts in each cell)"
  }
  
  
  
  gene_names <- df$Gene[1:nr_genes]
  df <- df[df$Gene %in% gene_names,]
  
  ggplot(df, aes(x = Percentage, y = Gene, fill = Gene)) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position='none') +
    ggtitle("Top most expressed genes") +
    xlab(x_legend) +
    ylab("")

}














