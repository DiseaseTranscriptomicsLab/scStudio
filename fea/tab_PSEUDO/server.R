
##-------------------------------------##
##            Pseudobulk TAB           ##
##-------------------------------------##

majority_vote <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

aggregate_column <- function(x) {
  if (is.numeric(x)) {
    return(mean(x, na.rm = TRUE))
  } else {
    return(majority_vote(x))
  }
}


make_pseudobulk <- function(count_matrix, metadata, group){
   
  count_matrix <- count_matrix[[1]]

  count_matrix <- as.matrix(count_matrix)
  
  counts <- 2^count_matrix - 1
  counts[counts < 0] <- 0
  
  group <- paste("pseudo", group, sep = "_")
  
  pseudobulk_counts <- sapply(unique(group), function(g) {
    rowSums(counts[, group == g, drop = FALSE])
  })
  
  cells_by_group <- split(seq_len(nrow(metadata)), group)
  
  pseudo_meta <- do.call(rbind, lapply(names(cells_by_group), function(g) {
    idx <- cells_by_group[[g]]
    df <- metadata[idx, , drop = FALSE]
    agg <- lapply(df, aggregate_column)
    agg <- as.data.frame(agg, stringsAsFactors = FALSE)
    rownames(agg) <- g
    agg
  }))
  
  colnames(pseudobulk_counts) <- make.unique(colnames(pseudobulk_counts))
  pseudo_meta <- cbind(colnames(pseudobulk_counts), pseudo_meta)
  colnames(pseudo_meta)[1] <- "sampleID"
   
  return(list(pseudobulk_counts, pseudo_meta))
}


run_marker <- function(organism, database, pathways, count_matrix,
                       metadata, condition, method){
  # 1. Get the reference table
  get_categories <- msigdbr_collections()
  
  print("Choosing gene sets")
  
  # 2. Check if the database is a main collection (H, C1, etc.)
  if (database %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) {  
    # Fetch by collection only
    gene_sets_raw <- msigdbr(species = organism, category = database)
  } else {
    # 3. Pull the parent collection string correctly
    # Use [[1]] or $gs_cat to ensure we get a character string, not a dataframe
    main_cat <- get_categories$gs_cat[get_categories$gs_subcat == database]
    
    if (length(main_cat) == 0) stop("Database subcollection not found!")
    
    print(paste("Parent collection identified as:", main_cat))
    
    # Fetch using both category and subcategory
    gene_sets_raw <- msigdbr(species = organism, category = main_cat, subcategory = database)
  }
  
  # 4. Split into the named list format required by markeR or fgsea
  gene_sets <- split(gene_sets_raw$gene_symbol, gene_sets_raw$gs_name)
  gene_sets <- gene_sets[names(gene_sets) %in% pathways]
  print(pathways)
  
    
  PlotScores(data = count_matrix, 
             metadata = metadata, 
             gene_sets = gene_sets, 
             Variable = condition,  
             method = method,   
             compute_cohen = TRUE,
             nrow = 1,    
             pointSize=4,  
             title="", 
             widthTitle = 24,
             labsize=12, 
             titlesize = 12)
  
}

run_marker_roc <- function(organism, database, pathways, count_matrix,
                       metadata, condition, method){
  # 1. Get the reference table
  get_categories <- msigdbr_collections()
  
  print("Choosing gene sets")
  
  # 2. Check if the database is a main collection (H, C1, etc.)
  if (database %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")) {  
    # Fetch by collection only
    gene_sets_raw <- msigdbr(species = organism, category = database)
  } else {
    # 3. Pull the parent collection string correctly
    # Use [[1]] or $gs_cat to ensure we get a character string, not a dataframe
    main_cat <- get_categories$gs_cat[get_categories$gs_subcat == database]
    
    if (length(main_cat) == 0) stop("Database subcollection not found!")
    
    print(paste("Parent collection identified as:", main_cat))
    
    # Fetch using both category and subcategory
    gene_sets_raw <- msigdbr(species = organism, category = main_cat, subcategory = database)
  }
  
  # 4. Split into the named list format required by markeR or fgsea
  gene_sets <- split(gene_sets_raw$gene_symbol, gene_sets_raw$gs_name)
  gene_sets <- gene_sets[names(gene_sets) %in% pathways]
  print(pathways)
  
  
  ROC_Scores(data = count_matrix, 
             metadata = metadata, 
             gene_sets=gene_sets, 
             method = method, 
             variable = condition,
             colors = c(logmedian = "#3E5587", ssGSEA = "#B65285", ranking = "#B68C52"), 
             grid = TRUE, 
             spacing_annotation=0.3, 
             ncol=NULL, 
             nrow=1,
             mode = "simple",
             widthTitle = 28,
             titlesize = 10,  
             title="") 
  
}





  
