##-------------------------------------##
#          VAREXPLAINED TAB            ##
##-------------------------------------##


plot_variance_explained <- function(mat, vars, df_metadata){
  
  df_metadata <- df_metadata[,vars]
  
  mat <- log2(mat+1)
  
  print("Calculating variance explained by each variable.")
  
  varMat <- getVarianceExplained(mat, variables = df_metadata)
  
  print("Calculation done.")
  
  p <- plotExplanatoryVariables(
                                varMat,
                                variables = variables)  + 
    ggtitle("") +
    
    theme(text = element_text(size=20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  
 return(p)
  
}


