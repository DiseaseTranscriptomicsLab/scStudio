##-------------------------------------##
##               BEC TAB               ##
##-------------------------------------##


correct_be <- function(countMatrix, method, batch_label1, batch_label2, keep_biological, covariate){

  if (method == "limma"){corrected_mat <- correct_limma(countMatrix, 
                                                        batch_label1,
                                                        batch_label2,
                                                        keep_biological,
                                                        covariate)}
 return(corrected_mat) 
}




correct_limma <- function(
  countMatrix, batch_label1, batch_label2, keep_biological, covariate){

  print(!is.null(keep_biological))
  
  if (!is.null(keep_biological)){
  
    design <- model.matrix(~ 0 + keep_biological)
  
  limma <- removeBatchEffect(countMatrix, 
                             batch = batch_label1, 
                             batch2 = batch_label2,
                             design = design,
                             covariates = covariate)}
  else {
  limma <- removeBatchEffect(countMatrix, 
                             batch = batch_label1, 
                             batch2 = batch_label2,
                             covariates = covariate)
  }

  return(limma)
}







