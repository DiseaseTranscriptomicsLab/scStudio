##-------------------------------------##
##              DOUBLETS               ##
##-------------------------------------##

#sce <- SingleCellExperiment(list(counts=mat))

#doublets <- cxds(sce,retRes = TRUE)

#doublets$cxds_score



calculate_doublets <- function(mat, sample){
  set.seed(2024)
  sample <- as.character(sample)
  
 doublets <- scDblFinder(sce = mat, samples = sample)
 
 df1 <- data.frame(
                   doublets = doublets$scDblFinder.class,
                   score = doublets$scDblFinder.score)
 
 #uniq_sample <- unique(sample)
 #total_cells <- as.vector(table(sample))
 #sample_doublets <- as.data.frame(table(sample, df1$doublets))
 #sample_doublets <- sample_doublets[sample_doublets$Var2 == "doublet","Freq"]
 #percentage_doublets <- sample_doublets/total_cells*100
# 
# df2 <- data.frame(sample = uniq_sample,
#                   total_cells = total_cells,
#                   sample_doublets = sample_doublets,
#                   percentage_doublets = round(percentage_doublets, 2))
# colnames(df2) <- c("Sample", "Total cells", "#doublets", "%doublets")
 
 return(df1)
}

get_doublets_table <- function(sample, metadata){
  
  uniq_sample <- unique(metadata[[sample]])
  total_cells <- as.vector(table(metadata[[sample]]))
  sample_doublets <- as.data.frame(table(metadata[[sample]], metadata[["doublets_class"]]))
  sample_doublets <- sample_doublets[sample_doublets$Var2 == "doublet","Freq"]
  percentage_doublets <- sample_doublets/total_cells*100
  
  df2 <- data.frame(sample = uniq_sample,
                    total_cells = total_cells,
                    sample_doublets = sample_doublets,
                    percentage_doublets = round(percentage_doublets, 2))
  
  colnames(df2) <- c("Sample", "Total cells", "#doublets", "%doublets")
  
  return(df2)
  
}


arrange_doublet_plots <- function(df_qc, scale, sample){
  
  
  
  df_qc$library_size <- select_scale(df_qc, "library_size", scale)
  df_qc$total_features <- select_scale(df_qc, "total_features", scale)
  
  df_qc <- df_qc[order(df_qc$doublets_score, decreasing = FALSE),]
  print(colnames(df_qc))
  df_qc[[sample]] <- str_sub(df_qc[[sample]], 1, 20)
  
  text <- paste("Class: ", df_qc[["doublets_class"]],
                paste0("\nsample: ", df_qc[["orig.ident"]]),
                paste0("\nScore: ",round(df_qc[["doublets_score"]], 2)))
  
  p <- ggplot(data = df_qc, aes(x = library_size, 
                                y = total_features,
                                text = text,
                                color = doublets_score)) + 
    
    geom_point( size = 0.3, alpha = 1) + 
    
    scale_color_gradient2(midpoint = 0.5, low = "darkgreen", mid = "yellow", high="red") + 
    
    xlab("Library size") +
    
    ylab("Total features") +
    
    labs(color = "Doublet score", caption = sample) +
    
    theme(text = element_text(size=20), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
   facet_wrap(~ df_qc[[sample]])
  
  return(p)
}


