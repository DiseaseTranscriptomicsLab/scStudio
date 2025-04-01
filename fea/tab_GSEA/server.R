
##-------------------------------------##
##               GSEA TAB              ##
##-------------------------------------##

get_gsea <- function(dea, ordering, organism, subcategories){
  dea[["|SNR|"]] <- abs(dea$SNR)
  
  dea <- dea[order(dea[[ordering]], decreasing = TRUE),]
  
  rank <- as.data.frame(cbind(dea$Gene, as.numeric(dea[[ordering]])))
  colnames(rank) <- c("genes", "metric")
  rank$metric <- as.numeric(rank$metric)
  rank <- rank[!is.na(rank$metric),]
  
  #nrank <- as.numeric(as.character(rank[,2]))
  nrank <- as.numeric(rank$metric)
  names(nrank) <- rank$genes
  
  
  chosen_gene_sets <- data.frame()
  get_categories <-  readRDS(paste0(getwd(), "/","/msigdbr_collections.rds"))
  
  for (subcat in subcategories){
    
    if (subcat %in% c("H", "C1", "C6", "C8"))
    { chosen_gene_sets <-  rbind(chosen_gene_sets, msigdbr(species = organism, category = subcat))}
    else{
      cat <- get_categories[get_categories$gs_subcat == subcat, "gs_cat"]
    chosen_gene_sets <- rbind(chosen_gene_sets, msigdbr(species = organism, category = cat, subcategory = subcat))}
  }
  
  
  chosen_gene_sets <- chosen_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)
  
  gsea_result <- fgsea(pathways = chosen_gene_sets, eps = 0, stats = nrank)
  
  results_tidy <- gsea_result %>% as_tibble() %>% arrange(desc(NES))
  
  results_tidy[,2] <- round(results_tidy[,2], 4)
  results_tidy[,3] <- round(results_tidy[,3], 4)
  results_tidy[,4] <- round(results_tidy[,4], 4)
  results_tidy[,5] <- round(results_tidy[,5], 4)
  results_tidy[,6] <- round(results_tidy[,6], 4)
  
  
  
  return(list(results_tidy, nrank))
}

plot_gsea <- function(tab, pval_cut, NES_cut, ID){
  
  tab$pathway <- str_sub(tab$pathway, 1, 35)
  
ggplot(tab %>% filter(padj < pval_cut & abs(NES) >= NES_cut), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES > 0)) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score",
       title= paste(ID),
       caption = paste0("BH-adjusted p-value < ",pval_cut, " and abs(NES) >= ", NES_cut)) + 
  theme(text = element_text(size = 18), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
} 





 


  