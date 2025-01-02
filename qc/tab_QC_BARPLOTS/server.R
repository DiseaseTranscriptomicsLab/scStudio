##-------------------------------------##
##             BARPLOTS                ##
##-------------------------------------##

get_summary_barplot <- function(var1, var2, colors, label){
  
  var1 <- str_sub(var1, 1, 15)
  
  var1_terms <- unique(var1)
  var2_terms <- unique(var2)
  var1_size <- length(var1_terms)
  var2_size <- length(var2_terms)
  
  plotting_df <- data.frame(group1 = factor(rep(var1_terms, each = var2_size)),
                            group = factor(rep(var2_terms, var1_size)))
                              
  mat <- as.matrix(table(var1, var2))     
  mat_proportion <- mat/rowSums(mat)*100
  mat_proportion <- mat_proportion[var1_terms,var2_terms]
  p_v <- c()
  for (row in 1:nrow(mat_proportion)){
    p_v <- c(p_v, mat_proportion[row,])}
  
  plotting_df$proportion <- p_v
  print("Dataframe for plotting proportions:")
  print(plotting_df)
  
  p <- ggplot(data=plotting_df, aes(x = group1, y = proportion, fill = group, text = proportion)) +
    
    geom_bar(stat="identity") +
    
    theme(legend.position="bottom",
          legend.title = element_blank(),
          text = element_text(size=25),
          axis.text.x = element_text(size = 15,angle = 45, hjust = 1),
          axis.ticks.y=element_blank()) +
    
    ggtitle(paste0("Proportion of ", label)) +
    
    scale_fill_manual(values=colors) +
    
    xlab("") + ylab("Cell proportion (%)") + labs(fill = "Group")
    
    return(p)
}




