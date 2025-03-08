
##-------------------------------------##
##              DEA TAB                ##
##-------------------------------------##

# OPTIMIZE get_groups_dea with DPLYR 

get_dea_allMethods <- function(ID, mat, group1, group2, methods, token){
    
  library(spatstat.core)
  library(Seurat)
  
  ram <- memuse::Sys.meminfo()
    
  while (ram$freeram@size < 20){
      
    print("Available RAM:")
      
    print(print(ram$freeram@size))
      
    print("Waiting for available memory to run DEA.")
      
    Sys.sleep(30)
      
    ram <- memuse::Sys.meminfo()
  }

    
  get_dea <- function(mat, method, group1, group2){
      
      norm_mat <- (2**as.matrix(mat))-1
      
      keep_genes <- rowSums(norm_mat) > 10
        
      norm_mat <- norm_mat[keep_genes,]
    
      ln_mat <- log(norm_mat + 1)
      
      srt <- Seurat::CreateSeuratObject(counts = ln_mat, min.cells = 0, min.features = 0)
      
      group1_names <- colnames(srt)[group1]
      group2_names <- colnames(srt)[group2]
      markers <- Seurat::FindMarkers(object = srt,
                                     ident.1 = group1_names,
                                     ident.2 = group2_names,
                                     logfc.threshold = -Inf,
                                     test.use = method,
                                     random.seed = 1234
                                     )
      
      log2FC_mean_all_genes <- c()
      log2FC_median_all_genes <- c()
      snr_all_genes <- c()
      print("Organizing final table...")
      for (n in 1:dim(markers)[1]){
        gene_name <- rownames(markers)[n]
        
        meanGroup1 <- log2(mean((2**mat[gene_name, group1])-1) + 10^-9)
        meanGroup2 <- log2(mean((2**mat[gene_name, group2])-1) + 10^-9)
        
        medianGroup1<- log2(median((2**mat[gene_name, group1])-1) + 10^-9)
        medianGroup2 <- log2(median((2**mat[gene_name, group2])-1) + 10^-9)
        
        sigma1 <- sd(2**mat[gene_name, group1])
        sigma2 <- sd(2**mat[gene_name, group2]) 
        
        log2FC_mean <- meanGroup1 - meanGroup2
        #if (is.na(log2FC_mean)){log2FC_mean <- 0} 
        log2FC_median <- medianGroup1 - medianGroup2
        
        snr <- (meanGroup1 - meanGroup2)/(sigma1 + sigma2 + 10^-9)  # Small constant to avoid division by zero
        
        log2FC_mean_all_genes <- c(log2FC_mean_all_genes, round(log2FC_mean,3))
        log2FC_median_all_genes <- c(log2FC_median_all_genes, round(log2FC_median, 3))
        snr_all_genes <- c(snr_all_genes, round(snr, 3))
      }
       
    
      
      if (method == "roc"){
        final_markers <- data.frame(
          gene = rownames(markers),
          AUC = markers$myAUC,
          Power = markers$power,
          pct.1 = round((markers$pct.1)*100, 3),
          pct.2 = round((markers$pct.2)*100,3),
          log2FC_mean = log2FC_mean_all_genes,
          log2FC_median = log2FC_median_all_genes,
          SNR = snr_all_genes 
        )
        
        colnames(final_markers) <- c("Gene",
                                     "AUC",
                                     "Power",
                                     "Percentage expressing cells (Target)",
                                     "Percentage expressing cells (Other)",
                                     "Log2FC (mean)",
                                     "Log2FC (median)",
                                     "SNR")
      }
      
      else {
        final_markers <- data.frame(
          gene = rownames(markers),
          adj_pvalue = markers$p_val_adj,
          pct.1 = round((markers$pct.1)*100,3),
          pct.2 = round((markers$pct.2)*100,3),
          log2FC_mean = log2FC_mean_all_genes,
          log2FC_median = log2FC_median_all_genes,
          SNR = snr_all_genes 
        )
        
        colnames(final_markers) <- c("Gene",
                                     "Adjusted p-value",
                                     "Percentage expressing cells (Target)",
                                     "Percentage expressing cells (Other)",
                                     "Log2FC (mean)",
                                     "Log2FC (median)",
                                     "SNR")
      }
      
      return(final_markers)
    } 
    
    for (method in methods){
      print("Method:")
      print(method)
      final_markers <- get_dea(mat, method, group1, group2)
      name <- paste(ID,method,sep = "_")
      session_obj_dea <- readRDS(paste0(getwd(),"/tokens/", token, "/dea.rds"))
      session_obj_dea[[name]] <- final_markers
      
      saveRDS(session_obj_dea, 
              paste0(getwd(),"/tokens/", token, "/dea.rds"))
    }
   } 
  

 
get_groups_dea <- function(group, grouping, metadata){
  final_group <- c()
  
  for (var in colnames(metadata)){
  if(is.factor(metadata[[var]])){
    print("Detected factor, converting to character.")
    metadata[[var]] <- as.character(metadata[[var]])
  }
  }
  
  if (grouping == "Union"){
  
    for (i in 1:dim(metadata)[1]){ #for each cell/row
        for (cond in group){ #for each condition/col  
       
         if (cond %in% unlist(metadata[i,])){
     
           final_group <- c(final_group, i)
          
         } #close if
        } # close for
       } # close for
      } #close if
  else { # intersection
    for (i in 1:dim(metadata)[1]){ #for each cell/row

      check <- c()
      for (cond in group){ #for each condition/col 
      
        check <- c(check,cond %in% unlist(metadata[i,]))
      } # close for
   
        if (FALSE %in% check){
         #pass 
        } #close if
        else {final_group <- c(final_group, i)}
      
    } # close for
  } #close else
 return (unique(final_group))
}


make_volcano <- function(df, ID, gene_label, pval, logFC) {
  
  df$pval <- -log10(df[["Adjusted p-value"]] + 10**-100)
  df$logfc <- df[["Log2FC (mean)"]]
  
  gene_names <- df$Gene
  
  ggplot(df) +
    geom_point(aes(text= gene_names, x = logfc, y = pval, col = gene_label), alpha = 0.5) + 
    theme_minimal() + 
    ggtitle(ID) +
    xlab("Average Log2FC") + 
    ylab("-log10(Adjusted p-value)") +
    theme(legend.position = "none", text = element_text(size=15))  +
    geom_label_repel(aes(x = logfc, y = pval, 
                         label = ifelse(gene_label == "UP" | gene_label == "DOWN",  
                                                      "" ,"")), size = 5,
                      max.overlaps = 5) + 
    scale_color_manual(values=c("grey45", "red", "green")) +
    
    geom_hline(yintercept = pval, linetype="dashed", 
               color = "blue", size = 0.5) +
    
    geom_vline(xintercept = logFC, linetype="dashed", 
               color = "blue", size = 0.5) +
    
    geom_vline(xintercept = -logFC, linetype="dashed", 
               color = "blue", size = 0.5)
}

#make_heatmap <- function(mat, genes, scale, cluster, token, selected_cols, group1, group2){
#  if (cluster == "none"){
#    cluster <- NULL
#    clust_rows <- FALSE
#    clust_cols <- FALSE
#  }
#  else if (cluster == "row"){
#    clust_rows <- TRUE
#    clust_cols <- FALSE
#  }
#  else if (cluster == "column"){
#    clust_rows <- FALSE
#    clust_cols <- TRUE
#  }
#  else {
#    clust_rows <- TRUE
#    clust_cols <- TRUE
#  }
#  
#  mat <- mat[genes,c(group1,group2)]#
#
#  
#  cols <- rep("grey", ncol(mat))
#  cols[c(1:length(group1))] <- selected_cols[1]
#  cols[c(length(group1)+1:length(group2))] <- selected_cols[2]#
#
#  heatmap.2(x= mat, 
#            dendrogram = cluster,
#            scale = scale,
#            col="bluered",
#            Rowv= clust_rows, 
#            Colv = clust_cols,
#            trace="none",
#            ColSideColors = cols,
#            labRow=rownames(mat),
#            main="",
#            ylab="",
#            key.title = NA,
#            keysize = 1.5,
#            labCol = FALSE,
#            cexRow = 1.5,
#            margins = c(10, 20),
#            xlab="")
#  
#  legend("topright", title = "Group",legend= c("Group 1", "Group 2"), 
#         fill = selected_cols, cex=1, box.lty=0)
#}


get_top_genes <- function(tb, nr){
  tb <- tb[order(tb[["Log2FC (mean)"]], decreasing = TRUE),]
  genes <- tb$Gene[1:nr]
  tb <- tb[order(tb[["Log2FC (mean)"]], decreasing = FALSE),]
  genes <- c(genes, tb$Gene[1:nr])
}

make_dea_scatter <- function(dea, cutoff){
  dea <- dea[dea[["Percentage expressing cells (Target)"]]> as.numeric(cutoff),] 
  gene_names <- dea$Gene
 p <- ggplot(dea, 
             aes(text= gene_names, 
                 x=as.numeric(dea[["Log2FC (mean)" ]]), 
                 y = dea[["Percentage expressing cells (Target)"]])) + 
    geom_point(color = "black", size = 0.5, alpha = 0.5) +
    xlab("Average expression (log(counts+1))") +
    ylab(paste0("Expressing cells (%) - truncated at ",cutoff,"%")) +
    theme_classic() + theme(text = element_text(size = 12)) 
}

remove_method_from_id <- function(name){
  name <- gsub("_wilcox$","", name)
  name <- gsub("_t$","", name)
  name <- gsub("_MAST$","", name)
  name <- gsub("_roc$","", name)
  return(name)
}


make_heatmap <- function(mat, genes, group1, group2, cluster, scale){
  
    if (cluster == "none"){
      cluster <- NULL
      clust_rows <- FALSE
      clust_cols <- FALSE
    }
    else if (cluster == "row"){
      clust_rows <- TRUE
      clust_cols <- FALSE
    }
    else if (cluster == "column"){
      clust_rows <- FALSE
      clust_cols <- TRUE
    }
    else {
      clust_rows <- TRUE
      clust_cols <- TRUE
    }

metadata <- rep("remove", ncol(mat))  

print("HERE1")

mat <- mat[genes,c(group1,group2)]

print("HERE2")
  
Group <- c("green","red")

names(Group) <- c("Group 1","Group 2")

col <- list(Group = Group)

metadata[group1] <- "Group 1"

metadata[group2] <- "Group 2"

metadata <- metadata[metadata != "remove"]

print("HERE3")
  
col_metaData <- data.frame(Group = metadata)

print("HERE4")
rownames(col_metaData) <- colnames(mat)

print("HERE5")

col_metaData <- col_metaData[order(col_metaData$Group),, drop = FALSE]

mat <- mat[, order(col_metaData$Group)]

colnames(mat) <- order(as.character(seq(1, ncol(mat))))

rownames(col_metaData) <- colnames(mat)

if (scale == "none") {legend_name <- "log2(counts + 1)"} 
else if (scale == "row") {
  legend_name <- "Scaled counts"
  mat <- t(scale(t(mat)))
  mat <- scale_matrix_rows(mat)
}
else{ 
  legend_name <- "Scaled counts"
  mat <- scale_matrix_columns(mat)
}

print("HERE6")
print(rownames(mat))

ggheatmap::ggheatmap(data = mat,
          color = colorRampPalette(c( "#0000ff","#fad541","#b60404"))(100),
          cluster_rows = clust_rows,
          cluster_cols = clust_cols,
          #scale = scale,
          legendName = legend_name,
          text_position_rows = "left",
          annotation_cols = col_metaData,
          annotation_color = col
)   %>% ggheatmap_theme(1:2,theme =list(
  
    theme(axis.title.x=element_blank(),
          text = element_text(size=20),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()),
    
    theme(
          legend.title = element_text(size = 18),
          legend.text = element_text(size= 18))
    )
  )
}

min_max_scale <- function(x, new_min = 0, new_max = 1) {
  (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
}

scale_matrix_rows <- function(mat, new_min = -1, new_max = 1) {
  t(apply(mat, 1, min_max_scale, new_min = new_min, new_max = new_max))
}

#scale_matrix_columns <- function(mat, new_min = -1, new_max = 1) {
#  apply(mat, 2, min_max_scale, new_min = new_min, new_max = new_max)
#}

scale_matrix_columns <- function(mat, new_min = -1, new_max = 1) {
  apply(mat, 2, function(x) {
    if (length(unique(x)) == 1) {
      rep(new_min, length(x))  # Set all values to new_min (0)
    } else {
      (x - min(x)) / (max(x) - min(x)) * (new_max - new_min) + new_min
    }
  })
}



