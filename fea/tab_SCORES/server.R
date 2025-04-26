
##-------------------------------------##
##             SCORES TAB              ##
##-------------------------------------##

get_score <- function(mat, gene_set, name){
  gene_set <- unlist(strsplit(gene_set, ' '))
 
  mat <- t(mat[gene_set,])
    
    
    mat <- scale(mat, center = colMedians(mat), scale = FALSE)
    score <- rowSums(mat)
    print(paste("Calculated scores for gene set", name))
  
    return(score)
}

get_score_pathways <- function(mat_name, mat, organism, database, pathways, name){
  get_categories <-  msigdbr_collections()
  
  print("Choosing gene sets")
  if (database %in% c("H", "C1", "C6", "C8")){  
    gene_sets <- msigdbr(species = organism, category = database)
    gene_sets <- gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
  }
  else{
    cat <- get_categories[get_categories$gs_subcat == database, "gs_cat"]
    gene_sets <- msigdbr(species = organism, category = cat, subcategory = database)
    gene_sets <- gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
  }
  
  scores <- list()
  report <- list()
  
  print("Calculating score")
  for (pathway in pathways){
    gene_set <- gene_sets[[pathway]]
    
    gene_set <- gene_set[gene_set %in% rownames(mat)]
    report[[paste(name, pathway, sep ="_")]][["genes"]] <- gene_set
    report[[paste(name, pathway, sep ="_")]][["matrix"]] <- mat_name
    report[[paste(name, pathway, sep ="_")]][["organism"]] <- organism
    report[[paste(name, pathway, sep ="_")]][["category"]] <- database
    
    scores[[paste(name, pathway, sep ="_")]] <- get_score(mat, gene_set, name)
  
  }
    return(list(scores, report))
  }


plot_gset_heatmap <- function(scores, ids, cluster, scale, selected_cols, groups){
 
  print(sessionInfo())
  
  df <- data.frame(cells = names(scores[[1]]))
  for (score in ids){
    #print(score)
    df[[score]] <- scores[[score]]
  }

  rownames(df) <- df$cells
  df$cells <- NULL
  
  mat <- as.matrix(t(df))
  print(dim(mat))
  
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
  
  print("Building heatmap...")

  colors <- c(selected_cols, "green", "red", "blue", "yellow", "orange",
              "pink", "cyan", "black", "grey", "brown", "purple")
  colors <- colors[1:length(unique(groups))]
  Group <- colors
  
  names(Group) <- unique(groups)
  
  col <- list(Group = Group)
  
  col_metaData <- data.frame(Group = groups)
  
  rownames(col_metaData) <- colnames(mat)
  
  col_metaData <- col_metaData[order(col_metaData$Group),, drop = FALSE]
  
  mat <- mat[, order(col_metaData$Group)]
  
  colnames(mat) <- order(as.character(seq(1, ncol(mat))))
  
  rownames(col_metaData) <- colnames(mat)
  
  if (scale == "none") {legend_name <- "log2(counts + 1)"
  } else {legend_name <- "Z-score"}
  
  p <- ggheatmap::ggheatmap(data = mat,
                       color = colorRampPalette(c( "#0000ff","#fad541","#b60404"))(100),
                       cluster_rows = clust_rows,
                       cluster_cols = clust_cols,
                       scale = scale,
                       legendName = legend_name,
                       text_position_rows = "left",
                       annotation_cols = col_metaData,
                       annotation_color = col
  ) 
  
  g <- ggheatmap_theme(ggheatmap = p,
                       plotlist = c(1,2), 
                       theme =list(
                         
                         theme(axis.title.x=element_blank(),
                               text = element_text(size=20),
                               axis.text.x=element_blank(), 
                               axis.ticks.x=element_blank()),
                         
                         theme(
                           legend.title = element_text(size = 18),
                           legend.text = element_text(size= 18))
                       )
  )
  g

}

plot_pathway <- function(selected_pathway, nrank, species){

  pathways <- readRDS(paste0(getwd(),"/",species, "_pathways.rds"))
  print(species)
  print(nrank)
  print(pathways[[selected_pathway]])
  plotEnrichment(pathways[[selected_pathway]],
                 nrank ) + labs(title = selected_pathway)
}

 


  
