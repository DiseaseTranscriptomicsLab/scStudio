
##-------------------------------------##
##             SUBSET TAB              ##
##-------------------------------------##

# SEE DEA TAB WARNING

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



plot_dimRed <- function(coord, x, y, metadata, sel_var, title, 
                        pointSize, alpha = 1, cols)
{
  
  text <- paste("Condition: ", metadata[,sel_var],
                paste0("\nDim",x,": "), coord[,x],
                paste0("\nDim",y,": "), coord[,y])
  
  p <- ggplot(data = coord, 
              aes(x = coord[,x], 
                  y = coord[,y],
                  color = metadata[,sel_var],
                  text = text,
                  key = 1:nrow(coord))) +  
    geom_point(size = pointSize, alpha = alpha) +   
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
  

  
  return(p)
}


get_random_token <- function(){
  set.seed(as.numeric(Sys.time()))
  token <- stri_rand_strings(1, 10, pattern = "[a-z0-9]")
  print(token)
  return(token)
}





