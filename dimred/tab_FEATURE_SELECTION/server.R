##-------------------------------------##
##         FEATURE SELECTION TAB       ##
##-------------------------------------##

calculate_hvg <- function(mat){
  hvg <- modelGeneVar(mat)
  hvg <- as.data.frame(hvg)
  hvg <- hvg[order(hvg$bio, decreasing = TRUE),]
  return(hvg)
  }

plot_varvsmean <- function(hvg, top, pval, gene){
  Gene <- rownames(hvg)
  
  gene_label <- c(rep(TRUE, top), rep(FALSE, dim(hvg)[1]-top)) &
    hvg$FDR <= pval
  
  hvg_highlight <- hvg[rownames(hvg) == gene,c("bio", "FDR")]
  hvg_highlight$FDR <- -log10(hvg_highlight$FDR)
  
  ggplot(hvg, aes(x=mean, y=total, label = Gene)) +  
  
    geom_point(color = ifelse(gene_label, "green3", "black"), size = 0.2, alpha = 0.5) +
  
    geom_smooth(method = "loess") +
    
  xlab("Mean of log-expression") +
  ylab("Variance of log-expression") +
  theme(text = element_text(size=20), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +

    annotate("point", x = hvg_highlight$bio, y = hvg_highlight$FDR, colour = "red")
  
  }


plot_pvalvsbio <- function(hvg, top, pval, gene){
  Gene <- rownames(hvg)
  gene_label <- c(rep(TRUE, top), rep(FALSE, dim(hvg)[1]-top)) &
    hvg$FDR <= pval
  
  hvg_highlight <- hvg[rownames(hvg) == gene,c("bio", "FDR")]
  hvg_highlight$FDR <- -log10(hvg_highlight$FDR)
  
  ggplot(hvg, aes(x=bio, y=-log10(FDR), label = Gene)) + 
    xlab("Biological component (variance)") +
    ylab("-log10(adjusted p-value)") +
    theme(text = element_text(size=20), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    
    geom_point(color = ifelse(gene_label, "green3", "black"), size = 0.2, alpha = 0.5) +
    geom_hline(yintercept=-log10(pval), linetype="dashed", color = "blue") +
    annotate("point", x = hvg_highlight$bio, y = hvg_highlight$FDR, colour = "red")
}



 
 


