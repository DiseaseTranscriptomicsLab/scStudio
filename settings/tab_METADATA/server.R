
##-------------------------------------##
##            METADATA TAB             ##
##-------------------------------------##


replace_label <- function(var, old_label, new_label){
  var <- gsub(paste0("^", old_label,"$"), new_label, as.character(var))
  return(var)
}







