
##-------------------------------------##
##            SEURAT TAB               ##
##-------------------------------------##



run_graph <- function(mat, pca, genesvspcs, features, dims, token, session_obj, ID){
  
  print("Running clustering...")
  
  norm_mat <- (2**as.matrix(mat))-1
  
  ln_mat <- log(norm_mat + 1)
  
  srt <- Seurat::CreateSeuratObject(counts = ln_mat, min.cells = 0, min.features = 0)
  
  srt[["pca"]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(pca), key = "pca_", assay = Seurat::DefaultAssay(srt))
  
  if (genesvspcs == "PCs"){
    srt <- Seurat::FindNeighbors(srt,
                                 reduction = "pca",
                                 dims = dims,
                                 compute.SNN = TRUE)
  } #close if 
  else {
    srt <- Seurat::FindNeighbors(srt,
                                 features = features,
                                 compute.SNN = TRUE)
  } #close else
  
  srt <- Seurat::FindClusters(srt, 
                              resolution = seq(0.1, 2, 0.1), 
                              random.seed = 1234,
                              verbose = FALSE)
  
  clusters <- as.data.frame(srt@meta.data)
  clusters <- clusters[, grep("RNA_snn_res", colnames(clusters))]
  colnames(clusters) <- paste0("res_",seq(0.1, 2, 0.1))
  
  session_obj[[ID]] <- clusters
  saveRDS(session_obj, 
          paste0(getwd(),"/tokens/", token, "/clusters.rds"))
  
  return(clusters)
} #close run_graph



#mat <- countMatrices$muSC_mat_norm
#norm_mat <- (2**mat)-1
#ln_mat <- log(norm_mat + 1)
#srt <- Seurat::CreateSeuratObject(counts = ln_mat, min.cells = 0, min.features = 0)
#ths <- 0.25
#all_clusters <- clusters$clustering_1
#session_obj <- clusters
#token <- "k15nly"
#ID <- "clustering_1"
#methods <- c("MAST")


run_marker_genes <- function(srt, mat, methods, ths, all_clusters, session_obj, token, ID, sres){
  library(spatstat.core)
  library(Seurat)
  for (res in sres){
    
    if(length(unique(all_clusters[[paste0("res_",res)]])) > 1){
    clusters <- all_clusters[[paste0("res_",res)]]
    
    Seurat::Idents(srt) <- clusters
    
    for (method in methods){
      markers <-  Seurat::FindAllMarkers(object = srt,
                                                         logfc.threshold = ths,
                                                         test.use = method,
                                                         random.seed = 1234)
      rownames(mat) <- rownames(srt)
      colnames(mat) <- colnames(srt)
      log2FC_mean_all_genes <- c()
      log2FC_median_all_genes <- c()

      
      for (n in 1:dim(markers)[1]){
  
        gene_name <- markers[n,]$gene
        
        cluster <- markers[n,]$cluster
        
        meanCluster <- log2(mean((2**mat[gene_name, clusters == cluster])-1) + 10^-9)
        meanOthers <- log2(mean((2**mat[gene_name, clusters != cluster])-1) + 10^-9)
        
        medianCluster <- log2(median((2**mat[gene_name,clusters == cluster])-1) + 10^-9)
        medianOthers <- log2(median((2**mat[gene_name,clusters != cluster])-1) + 10^-9)
        
        log2FC_mean <- meanCluster - meanOthers
        log2FC_median <- medianCluster - medianOthers
        
        log2FC_mean_all_genes <- c(log2FC_mean_all_genes, round(log2FC_mean,3))
        log2FC_median_all_genes <- c(log2FC_median_all_genes, round(log2FC_median, 3))
      } #close for
      
      
      if (method == "roc"){
        final_markers <- data.frame(
          gene = markers$gene,
          cluster = markers$cluster,
          AUC = markers$myAUC,
          Power = markers$power,
          pct.1 = (markers$pct.1)*100,
          pct.2 = (markers$pct.2)*100,
          log2FC_mean = log2FC_mean_all_genes,
          log2FC_median = log2FC_median_all_genes
        )
        
        colnames(final_markers) <- c("Gene",
                                     "Cluster",
                                     "AUC",
                                     "Power",
                                     "Percentage expressing cells (Target)",
                                     "Percentage expressing cells (Other)",
                                     "Log2FC (mean)",
                                     "Log2FC (median)")
      } #close if 
      
      else {
        
        final_markers <- data.frame(
          gene = markers$gene,
          cluster = markers$cluster,
          adj_pvalue = markers$p_val_adj,
          pct.1 = (markers$pct.1)*100,
          pct.2 = (markers$pct.2)*100,
          log2FC_mean = log2FC_mean_all_genes,
          log2FC_median = log2FC_median_all_genes
        )
      
       colnames(final_markers) <- c("Gene",
                                     "Cluster",
                                     "Adjusted p-value",
                                     "Percentage expressing cells (Target)",
                                     "Percentage expressing cells (Other)",
                                     "Log2FC (mean)",
                                     "Log2FC (median)")
      } #close else
      
      
      name <- paste(ID, method, res, sep = "_")
      session_obj[[name]] <- final_markers
      saveRDS(session_obj, 
              paste0(getwd(),"/tokens/", token, "/mkgs.rds"))
     
       
      } #close for 
     } #close if
    } #close for 
  
}#close function

plot_violin_clusters <- function(mat, clusters, gene){
  set.seed(1234)
  #text <- paste("Sample: ", df[["orig.ident"]], 
  #              "Dataset: ", df[["dataset"]])
  
  print("Plotting clusters violin-plot.")

  df <- data.frame(gene = mat[gene,], clusters = clusters)
  
  p <- ggplot(df, aes(x = clusters, 
                      y = gene, 
                      color = clusters,
                      fill = clusters)) + 
    geom_violin(size = 0.2, colour = "black", scale = "width") +
    geom_jitter(shape= 20, position=position_jitter(0.3),  
                alpha = 0.5, 
                size = 0.05, colour = "black") +
    stat_summary(fun.y=median, geom="point", size=2, color="red") + 
    stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
    #scale_fill_manual(values = c("red", "blue")) + 
    ggtitle(gene) +
    xlab("") + 
    ylab("log2(counts+1)") +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle=45, hjust=1), 
          text = element_text(size=17), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0,0,0,2), "cm"))
  return(p)
}
  
plot_clustree <- function(clusters, res){
     if (res < 1){
     clustree::clustree(
            clusters[,1:10], 
            prefix = "res_",
            return = "plot"
            ) +
            theme(text = element_text(size=10), 
                  legend.position="right",
                  legend.text=element_text(size=10)) + 
         
         labs(color = "Resolution", size = "Size") +
         
         scale_edge_alpha_continuous(name = "Proportion")  +
         guides(colour = guide_legend(override.aes = list(size=3)))
     }
  else {
    clustree::clustree(
           clusters[,10:20], 
           prefix = "res_",
           return = "plot") +
          theme(text = element_text(size=10), 
                legend.position="right",
                legend.text=element_text(size=15)) + 
      
      labs(color = "Resolution", size = "Size") +
      
      scale_edge_alpha_continuous(name = "Proportion")  +
      guides(colour = guide_legend(override.aes = list(size=3)))
    }
}

plot_dimRed_clusters <- function(coord, x, y, metadata, sel_var, title, size = 0.1, alpha = 1, cols){
  
  
  p <- ggplot(data = coord, 
              aes(x = coord[,x], 
                  y = coord[,y],
                  color = metadata[,sel_var],
                  key = 1:nrow(coord))) +  
    geom_point(size = size, alpha = alpha) +   
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
  
}











