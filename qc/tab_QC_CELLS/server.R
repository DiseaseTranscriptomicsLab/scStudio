##-------------------------------------##
##               QC - CELLS            ##
##-------------------------------------##

calculateQCmetrics <- function(countMatrix){
  
  print("Calculating CELL QC")
  
  mt_genes <- grep("^MT[-\\.]", rownames(countMatrix), ignore.case = TRUE, value = TRUE)

  lib_sizes <- colSums(countMatrix)
  total_features <- calculate_total_features(countMatrix)
  
  mt_sizes <- colSums(countMatrix[mt_genes,])
  
  mt_pct <- (mt_sizes/lib_sizes)*100
  
  df <- data.frame(library_size = lib_sizes,
                   total_features = total_features,
                   percentage_mt_genes = mt_pct)
  
  return (df)
}

calculate_total_features <- function(countMatrix){
  check_countMatrix <- countMatrix > 0
  total_features <- colSums(check_countMatrix)
  return(total_features)
}

select_scale <- function(df_qc, x, scale){
  if(scale == "Linear"){return(df_qc[[x]])}
  if(scale == "Log10"){return(log10(df_qc[[x]]))}
  if(scale == "Log2"){return(log2(df_qc[[x]]))}
}

select_vline <- function(x, scale){
  if(scale == "Linear"){return(x)}
  if(scale == "Log10"){return(log10(x))}
  if(scale == "Log2"){return(log2(x))}
}

plot_histogram <- function(df_qc, x ,bins, title, xlab, min, max, scale)
{ df_qc[[x]] <- select_scale(df_qc, x, scale)

  p<- gghistogram(data=df_qc, 
              x = x,
              y = "..count..",
              bins= bins,
              alpha = .2,
              color="black",
              fill="grey",
              title= title,
              ylab="#cells",
              xlab=xlab) + 
    theme_bw() +
    theme(text = element_text(size=20)) + 
    geom_vline(xintercept = min, linetype="dotted",
               color = "red", size=0.5) + 
    geom_vline(xintercept = max, linetype="dotted",
               color = "red", size=0.5)
  
} 


plot_scatterLibFeat <- function(df_qc, scale, minx, maxx, miny, maxy, leg, midpt)
{ #leg is not being used
   # future function: if radio button = mt_genes, plot this, if radio button = doublets, plot that
  #text <- paste("Condition: ", variables_df$condition, 
   #             "<br>Sample:", variables_df$sample,
    #            "<br>Cell type:", variables_df$cell_type,
     #           "<br>Cluster:", clusters)
  df_qc$library_size <- select_scale(df_qc, "library_size", scale)
  df_qc$total_features <- select_scale(df_qc, "total_features", scale)
  
  df_qc <- df_qc[order(df_qc$percentage_mt_genes, decreasing = FALSE),]
  
  if (scale == "Log2"){
    x_legend = "Log2(Total counts)"
    y_legend = "Log2(Total features)"
  }
  else if (scale == "Log10"){
    x_legend = "Log10(Total counts)"
    y_legend = "Log10(Total features)"
  }
  else({x_legend = "Total counts"
        y_legend = "Total features"})
  
  p <- ggplot(data = df_qc, aes(library_size, total_features, color = percentage_mt_genes)) + 
    geom_point( size = 1, alpha = 1) + #aes(text = text), size = size, alpha = alpha
    scale_color_gradient2(midpoint = midpt, low = "darkgreen", mid = "yellow", high="red") + 
    xlab(x_legend) +
    ylab(y_legend) +
    labs(color = "% MT-genes") +
    theme(text = element_text(size=20), 
          legend.title = element_text(size=14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
     geom_vline(xintercept = minx, linetype="dotted",
               color = "red", size=0.5) + 
     geom_vline(xintercept = maxx, linetype="dotted",
               color = "red", size=0.5) +
     geom_hline(yintercept = miny, linetype="dotted",
                color = "red", size=0.5) + 
     geom_hline(yintercept = maxy, linetype="dotted",
                color = "red", size=0.5)

}

plot_scatterLibMt <- function(df_qc, scale, minx, maxx, leg, midpt)
{ #leg is not being used
  df_qc$library_size <- select_scale(df_qc, "library_size", scale)
  
  df_qc <- df_qc[order(df_qc$percentage_mt_genes, decreasing = FALSE),]
  
  if (scale == "Log2"){
    x_legend = "Log2(Total counts)"

  }
  else if (scale == "Log10"){
    x_legend = "Log10(Total counts)"

  }
  else(x_legend = "Total counts")
  
  
  p <- ggplot(data = df_qc, aes(library_size, percentage_mt_genes, color = percentage_mt_genes)) + 
    geom_point( size = 1, alpha = 1) + #aes(text = text), size = size, alpha = alpha
    scale_color_gradient2(midpoint = midpt, low = "darkgreen", mid = "yellow", high="red") + 
    xlab(x_legend) +
    ylab("% MT-genes") +
    labs(color = "% MT-genes") +
    theme(text = element_text(size=20), 
          legend.title = element_text(size=14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
    geom_vline(xintercept = minx, linetype="dotted",
               color = "red", size=0.5) + 
    geom_vline(xintercept = maxx, linetype="dotted",
               color = "red", size=0.5)+
   geom_hline(yintercept = midpt, linetype="dotted",
             color = "red", size=0.5) 
  
}
 
      
plot_boxplot <- function(df,var, sample, title, ytitle, cols, scale, showCells){
  set.seed(1234)
  
  df[[sample]] <- str_sub(df[[sample]], 1, 15)
  
  if(scale =="Log2"){ytitle <- paste0("Log2(",ytitle,")")}
  else if(scale =="Log10"){ytitle <- paste0("Log10(",ytitle,")")}
  
  p <- ggplot(df, aes(x = df[[sample]], 
                      y = select_scale(df, var, scale),  
                      color = df[[sample]],
                      fill = df[[sample]],
                      text = "")) + 
    geom_violin(size = 0.2, colour = "black", scale = "width") +
    stat_summary(fun.y=median, geom="point", size=2, color="red") + 
    stat_summary(fun.y=mean, geom="point", size=2, color="blue") +
    scale_fill_manual(values = cols) + 
    ggtitle(title) +
    xlab("") + 
    ylab(ytitle) +
    theme(legend.position = "none", 
          axis.text.x = element_text(angle=45, hjust=1), 
          text = element_text(size=17), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  if (as.boolean(showCells)){
    p <- p + geom_jitter(shape= 20, position=position_jitter(0.3),  alpha = 0.5, size = 0.05, colour = "black")
  } #close if
  return(p)
}


