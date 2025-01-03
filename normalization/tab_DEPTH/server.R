##-------------------------------------##
##             DEPTH TAB               ##
##-------------------------------------##

choose_norm_method <- function(method, mat){
  print(method)
  
  if(method == "scran"){
    print("running scran")
    norm_matrix <- normalize_with_scran(mat)
  }
  
  else if(method == "SCTransform") {
    print("running sctransform")
    norm_matrix <- normalize_with_sctransform(mat)
  }
  
  return(norm_matrix)
}

normalize_with_sctransform <- function(mat){
  
  # remove genes with 0 total counts
  
  gene_counts <- rowSums(mat)
  
  mat <- mat[gene_counts > 0,]
  
  srt <- Seurat::CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
  
  options(future.globals.maxSize = 8000 * 1024^2)
  
  srt <- Seurat::SCTransform(srt)
  
  norm_mat <- log2(as.matrix(GetAssayData(srt, slot = "counts")) + 1)
  
  return(norm_mat)
}

normalize_with_scran <- function(mat){
  
  qclust <- quickCluster(mat)
 
  factors <- pooledSizeFactors(mat, clusters=qclust)
  
  print(summary(factors))
  norm_mat <- normalizeCounts(mat, 
                              size.factors = factors,
                              transform = "log",
                              pseudo.count = 1)
  return(norm_mat)
}

plot_boxplot <- function(df,var, sample, title, ytitle, yrange, cols){
  set.seed(1234)
  df[[sample]] <- str_sub(df[[sample]], 1, 20)
  
  p <- ggplot(df, aes(x = df[[sample]], 
                      y = log10(df[[var]]), 
                      color = df[[sample]],
                      fill = df[[sample]],
                      text = "")) + 
    geom_violin(size = 0.2, colour = "black", scale = "width") +
    geom_jitter(shape= 20, position=position_jitter(0.3),  
                alpha = 0.5, 
                size = 0.05, colour = "black") +
    stat_summary(fun.y=median, geom="point", size=2, color="red") + 
    stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
    scale_fill_manual(values = cols) + 
    ggtitle(title) +
    xlab("") + 
    ylim(yrange) +
    ylab(ytitle) +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle=45, hjust=1), 
          text = element_text(size=17), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  return(p)
}

