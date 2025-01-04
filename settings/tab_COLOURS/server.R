##-------------------------------------##
##            COLOURS TAB              ##
##-------------------------------------##

make_show_color_plot <- function(df){
  
  text <- df$colors
  
  print(text)

    ggplot(df, aes(x = group, y = add, fill = group,  text = paste(colors))) +
    
    geom_point(size = 20, colour="black", pch=21) +
    
    scale_fill_manual(values= df$colors) +
    
    scale_y_continuous(limits = c(0, 5)) +
    
    scale_x_discrete(guide = guide_axis(n.dodge = 100)) +
    
      theme(text = element_text(size=20),
            axis.line=element_blank(), 
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()) 
}

get_ordered_colors <- function(list_of_colors){
  if (!is.null(list_of_colors)){
  df <- data.frame(colors = unlist(list_of_colors),
                   group = names(list_of_colors),
                   add = rep(1, length(names(list_of_colors))))
  df <- df[order(df$group),]
  return(df$colors)
  }
}






