
# SELECT DATA TAB -----------------------------------------------------------------------------------------        

# UPLOAD PP -----------------------------------------------
upload_preprocessed <- function(name, type){
  dir_path <- paste0("./public_datasets/", name, "/", type, ".rds")
  obj <- readRDS(dir_path)
  return (obj)
}

# UPLOAD INHOUSE --------------------------------------------
upload_inhouse <- function(path){
  mat <- as.data.frame(vroom(path))
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]

  return(mat)
}

# GET GEO -----------------------------------------------------------------------------------------
get_geo <- function(geo_id){
  # Requires: geo_id is a string that contains a valid ncbi GEO ID for a scRNA-seq dataset entry
  # Returns: a list, for which the first entry is a list of count matrices and the second entry 
  # is a dataframe with metadata 
  
  print(paste("Retrieving",geo_id, "from GEO.", sep =" " ))
  
  # get metadata
  gse <- getGEO(geo_id)
  gse <- gse[[1]]
  metadata <- pData(gse)
  
  #folder_name <- paste(getwd(),geo_id, sep="/") #What is this for?
   
  # get the count matrices
  count_matrices <- get_GEOcountMatrix(geo_id)
  
  print("Download complete.")
  
  count_matrices_samples <-  str_extract(names(count_matrices), "^GS[EM][0-9]*")
  
  # metadata is usually stored at the sample level (tissue, patient, etc.) level, we need to 
  #convert at the cell-level
  df_sample_level <- data.frame(sample = count_matrices_samples,
                   number = unlist(lapply(count_matrices, ncol)))
  
  # we repeat the information vs the number of cells in each sample
  orig.ident <- c()
  for (original_sample in rownames(df_sample_level)){
    orig.ident <- c(orig.ident, rep(original_sample, df_sample_level[original_sample,"number"]))
  } #close for
  
  # deal with repeated information in the metadata annotation
  df_cell_level <- df_sample_level
  for(sample_name in unique(df_sample_level$sample)){
    if(length(grep(sample_name, df_sample_level$sample)) > 1){
      
      print("HERE")
      print(df_cell_level)
      print(sample_name)
      
      df_cell_level <-  dplyr::filter(df_cell_level, sample != sample_name)
      df_cell_level <- rbind(df_cell_level, c(sample_name, sum(df_sample_level[df_sample_level$sample == sample_name, "number"])))
      colnames(df_cell_level) <- c("sample", "number")
     } #close if
  } #close for loop
  
  # clean metadata - remove long and repetitive descriptions that can be viewed in the GEO entry
  column_names <- colnames(metadata)
  
  ## remove database information
  database_columns <- grep("status|taxid|submission|update|type|channel|relation", column_names)
  if (length(database_columns) != 0) {metadata <- metadata[,-database_columns]}
  column_names <- colnames(metadata)
  
  ## remove protocol information
  protocol_columns <- grep("protocol|molecule|platform|library_selection|library_source|library_strategy", column_names)
  if (length(protocol_columns) != 0) {metadata <- metadata[,-protocol_columns]}
  column_names <- colnames(metadata)
  
  ## remove contact information
  contact_columns <- grep("contact", column_names)
  if (length(contact_columns) != 0) {metadata <- metadata[,-contact_columns]}
  column_names <- colnames(metadata)
  
  ## remove supplementary_file information
  supplementary_columns <- grep("supplementary|data_row_count", column_names)
  if (length(supplementary_columns) != 0) {metadata <- metadata[,-supplementary_columns]}
  column_names <- colnames(metadata)
  
  ## remove processing information
  processing_columns <- grep("processing", column_names)
  if (length(processing_columns) != 0) {metadata <- metadata[,-processing_columns]}
  
  new_metadata <- list()
  
  # multiply metadata entries by the number of cells
  for (sample in df_cell_level$sample){ 
    new_metadata[[sample]] <- metadata[rep(grep(sample, metadata$geo_accession), df_cell_level[df_cell_level$sample == sample,"number"]),]
    } #close for loop
  
  if (nrow(new_metadata[[1]]) == 0){ #this happens when there is metadata for multiple samples  
                                       #but only one count matrix for the whole dataset
    new_metadata[[geo_id]] <- metadata[rep(1, ncol(count_matrices[[1]])),] 
    } #close if  
  
  new_metadata <- merge_matrices_by_row(new_metadata)
  
  # add the original sample ID and dataset information to metadata
  if (length(orig.ident) == nrow(new_metadata)){
    new_metadata$orig.ident <- orig.ident
    new_metadata$dataset <- rep(geo_id, dim(new_metadata)[1])
  } #close if 
  
  #if there is no annotation/metadata of samples in the count matrix files, just use orig.ident and dataset:
  else {
    new_metadata <- data.frame(orig.ident = orig.ident, dataset = rep(geo_id, length(orig.ident)))
  } #close else
  
  # save barcodes and create unique cell (key) identifier
  mats_barcodes <- c()
  mats_identifiers <- c()
  
  for (mat_name in names(count_matrices)){
  
    mat <- count_matrices[[mat_name]]
    
    print("Matrix size:")
    print(dim(mat))
    
    mats_barcodes <- c(mats_barcodes, colnames(mat))
    
    cell_ids <- define_cell_ids(ncol(mat))
    # for each cell, define a unique key
    #for (cell in 1:ncol(mat)){
    #  key <- gsub("[-:WET]","",Sys.time())
    #  key <- gsub(" ","", key)
    #  cell_ids <- c(cell_ids, paste0(key, stri_rand_strings(1, 6, pattern = "[a-z0-9]")))
    #}

    colnames(count_matrices[[mat_name]]) <- cell_ids 
    mats_identifiers <- c(mats_identifiers, cell_ids)
  } #close for loop
  
  
  new_metadata$barcodes <- mats_barcodes
  new_metadata$cell_id <- mats_identifiers
  
  print("Data preparation was successful.")
  
  unlink(paste0(getwd(),"/session_data/",geo_id), recursive = TRUE)

  return(list(count_matrices, new_metadata))
} #close get_geo function

# GET GEO matrix -----------------------------------------------------------------------------------------
get_GEOcountMatrix <- function(geo_id){
    # Requires: geo_id is a string that contains a valid ncbi GEO ID for a scRNA-seq dataset entry
    # Returns: matrices is a list of count matrices
  
    folder_name <- paste("./session_data/", geo_id, sep="")
    
    # scRNA-seq count matrices are stored in the supplementary files
    supp <- getGEOSuppFiles(geo_id, baseDir = "./session_data") 
    
    untar_everything(geo_id, folder_name)

    choose_workflow <- identify_workflow(folder_name)
    
    print(paste("Workflow identified:", choose_workflow))
    
      if (choose_workflow == 1){
        matrices <- workflow1(folder_name)
      }
      if (choose_workflow == 2){
        matrices <- workflow2(folder_name)
      }
      if (choose_workflow == 3){
        matrices <- workflow3(folder_name)
      }
      if (choose_workflow == 4){
        matrices <- workflow4(folder_name)
      }
      return(matrices)
} #close get_GEOcountMatrix

# WORKFLOW 1 -----------------------------------------------------------------------------------------
workflow1 <- function(dataset_path){
  # Requires: dataset_path is a string 
  # Returns: count_matrices is a list of count matrices
  # This workflow is applied when count matrices are stored in .csv or .txt files
  
  all_files <- list.files(path = dataset_path, 
                          #pattern = "matrix|counts",
                          ignore.case = TRUE)

  count_matrices <- list()
  
  for (file in all_files){
    
    print("Reading file:")
    print(file)
    
    tryCatch({
      
      file_path <- paste(dataset_path,file, sep = "/")
      
      sample_name <- str_extract(file, "^GS[EM][0-9]*_[aA-zZ1-9]*")
      
      # verify if the sample_name is not repeated
      if (sample_name %in% names(count_matrices)){sample_name <- paste0(sample_name, "_1")}
    
      #if sample name does not contain GEO ID as it should, replace with the file path
      if (is.na(sample_name)){sample_name <- file}  
      
      # check if file is xlsx
      if (identical(grep("xlsx", file), integer(1))){
        cm <- read_excel(file_path)
        rownames(cm) <- cm[,1]
      } #close if 
      
      else{
        cm <- as.data.frame(fread(file_path))
        rownames(cm) <- cm[,1]
      } #close else
      
      cm_clean <- remove_bad_columns(cm) 
      
      if (ncol(cm_clean) > 10 & nrow(cm_clean) > 100){ count_matrices[[sample_name]] <- cm_clean}
     
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #close tryCatch
  } #close for loop
  
  return(count_matrices)
} #close workflow1

# REMOVE BAD COLUMNS -----------------------------------------------------------------------------------------
remove_bad_columns <- function(mat){
  # Requires: mat is a count matrix
  # Returns: same count matrix without character/string columns 
  # (only keeps numeric values)
  
  mat <- as.data.frame(mat)
  mat_names <- colnames(mat)
  for (col in mat_names){
    if (is.na(as.numeric(mat[[col]]))){
      mat[[col]] <- NULL
    } #close if
    
    else{mat[[col]] <- as.numeric(mat[[col]])
    }#close else
  } #close for loop
  
  mat <- as.matrix(mat)
  
  return(mat)
} #close remove_bad_columns

# WORKFLOW 2 -----------------------------------------------------------------------------------------
workflow2 <- function(dataset_path){
  # Requires: dataset_path is a string 
  # Returns: count_matrices is a list of count matrices
  # This workflow is applied when count matrices are stored in the 10X format
  
  all_files <- list.files(path = dataset_path, pattern = "mtx")
  count_matrices <- list()
  subsample_names <- c()
  
  mat_names <- gsub("*\\.mtx.gz", "", all_files)
  
  other_names <- list.files(path = dataset_path ,pattern= ".+barcodes.[tc]sv.gz|.+cells.[tc]sv.gz")
  
  for (mat_name in mat_names){
  subsample_names <- c(subsample_names, 
                       find_sample_name(mat_name, other_names))
  } #close for loop
  
  subsample_names <- unique(subsample_names)
 
  for(sample in subsample_names){
    
    barcode_path <- paste(dataset_path,
      list.files(path = dataset_path ,pattern=paste(sample, ".+barcodes.[tc]sv|", sample ,".+cells.[tc]sv", sep ="")), sep ="/")
    
    features_path <- list.files(
      path = dataset_path, pattern=paste(sample, ".+features.[tc]sv|", sample,".+genes.[tc]sv", sep =""))
    
    # some datasets do not include a features files for each sample
    if (length(features_path) > 0){features_path <- paste(dataset_path,features_path, sep ="/")}
    else{
      # in that case, there should be a general features file for all count matrices
      print("No specific features file included.")
      features_path <- paste(dataset_path,
                             list.files(path = dataset_path ,pattern="genes|features"), sep ="/")
      } #close else
    
    matrix_path <- paste(dataset_path,
      list.files(path = dataset_path ,pattern=paste(sample, ".+mtx.gz", sep ="")), sep = "/")
    
    mtx <- readMM(file = matrix_path)
    print("Read matrix successfully.")
    
    # check if the barcodes file is in tsv or csv
    if (identical(grep("csv", barcode_path), integer(0))){
    print("reading tsv files")
    barcode.names = read.delim(barcode_path, 
                                header = FALSE,
                                stringsAsFactors = FALSE)
    } #close if 
    else {
      print("reading csv files")
      barcode.names = read.delim(barcode_path, 
                                     header = FALSE,
                                     stringsAsFactors = FALSE,
                                     sep = ",")
    } #close else
    print("Barcodes read successfully.")
    
    # check if the features file is in tsv or csv
    if (identical(grep("csv", features_path), integer(0))){
    feature.names = read.delim(features_path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    } #close if
    else {
      print("reading csv files")
      feature.names = read.delim(features_path, 
                                     header = FALSE,
                                     stringsAsFactors = FALSE,
                                     sep = ",")
    } #close else
    print("Features read successfully.")
    
    # attribute the barcodes to the cells
    step1 <- try(colnames(mtx) <- barcode.names$V1)
    
    # check if the features file is a list and number of rows matches the count matrix
    if(typeof(feature.names) == "list"){
      print("Feature.names is list.")
      
      if(dim(mtx)[1] == dim(feature.names)[1]){
        print("Dimension of feature.names matches.")
        
         if(!is.null(feature.names$V2)){step2 <- try(rownames(mtx) <- feature.names$V2)} #default case
        
         else {step2 <- try(rownames(mtx) <- feature.names[,1])}
           } #close second if 
      
      # check if the features file is a list and numbers don't match  
       else if(!is.null(feature.names$V2)){
         print("Dimension of feature.names doesn't match.")
         
         step2 <- try(rownames(mtx) <- feature.names$V2[-1])}
      
       else {step2 <- try(rownames(mtx) <- feature.names[-1,1])}  
        } #close first if 
    
    # if its a vector and numbers match
    if(typeof(feature.names) != "list"){
      print("Feature.names is not a list.")
      
      if(dim(mtx)[1] == length(feature.names)){
        print("Dimension of feature.names matches.")
        step2 <- try(rownames(mtx) <- feature.names)
        } #close second if 
      
      # if its a vector and numbers don't match  
      else{
        print("Dimension of feature.names doesn't match.")
        step2 <- try(rownames(mtx) <- feature.names[-1])}
       
    } #close first if 
    
    count_matrices[[sample]] <- mtx
  } #close for loop
  
  print("Retrieved count matrices sucessfully.")
  return(count_matrices)
} #close workflow2

# FIND SAMPLE NAME -----------------------------------------------------------------------------------------  
  find_sample_name <- function(target, query_strings){
    # Requires: target is string; query_strings is a vector of strings
    # Returns: longest common substring between target and any of the query_strings;
    #          if the last character of the string belongs to {-_.}, remove it.   
    
    results <- c()
    for (qs in query_strings){
      result <- ""
      
      for (cr in 1:str_length(target)){

        if (substr(target, cr,cr) == substr(qs, cr,cr)){
          result <- paste(result, substr(target, cr,cr), sep = "")
          if (cr == str_length(target)){results <- c(results, result)}

          } #close if 
        else {
          results <- c(results,result)
          break} #close else
      } #close second for
    }#close first for
    
    results_length <- str_length(results)
    
    results <- results[order(results_length, decreasing = TRUE)]
    
    bigger_string <- results[1]
    
    if (substr(bigger_string, str_length(bigger_string),
               str_length(bigger_string)) %in% c("_","-",".")){
      bigger_string <- substr(bigger_string, 1, str_length(bigger_string) - 1)
    } #close if 
    return(bigger_string)
} #close find_sample_name
  
# WORKFLOW 3 ----------------------------------------------------------------------------------------- 
workflow3 <- function(dataset_path){
  # Requires: dataset_path is a string 
  # Returns: count_matrices is a list of count matrices
  # This workflow is applied when count matrices are stored in the h5 format.
  
  all_files <- list.files(path = dataset_path, 
                          pattern = "h5",
                          ignore.case = TRUE)

  count_matrices <- list()
  for (file in all_files){
    
    tryCatch({
      
    file_path <- paste(dataset_path,file, sep = "/")
    sample_name <- str_extract(file, "^GS[EM][0-9]*_[aA-zZ1-9]*")
    
    cm <- Read10X_h5(file_path,use.names = TRUE, unique.features = TRUE)

    if (is.null(names(cm))){count_matrices[[sample_name]] <- cm}
    
    else (
      for (name in names(cm)){
      sub_cm <- cm[[name]]
      if (ncol(sub_cm) > 10 & nrow(sub_cm) > 100 ){
      count_matrices[[paste(sample_name, name, sep="_")]] <- sub_cm
      }
      } #close for loop
    ) #close else
    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #close tryCatch  
  } #close for loop
  return(count_matrices)
} #close workflow 3

# WORKFLOW 4 - DEPRECATED ----------------------------------------------------------------------------------------
workflow4 <- function(dataset_path){
  # Requires: dataset_path is a string 
  # Returns: count_matrices is a list of count matrices
  # This workflow is applied when count matrices are stored in the xlsx format.
  
  all_files <- list.files(path = dataset_path, 
                          pattern = "xlsx",
                          ignore.case = TRUE)
 
  count_matrices <- list()
  
  for (file in all_files){
    file_path <- paste(dataset_path,file, sep = "/")
    sample_name <- str_extract(file, "^GS[EM][0-9]*_[aA-zZ1-9]*")
    
    cm <- read_excel(file_path)

    count_matrices[[sample_name]] <- as.matrix(cm)
    } #close for loop
  return(count_matrices)
} #close workflow4

# WORKFLOW 5 -----------------------------------------------------------------------------------------
# workflow5: matlab
#cmat <- readMat(gzfile(
#"/mnt/scratch/home/marta/shiny/scStudio_2.0/session_data/GSE110823/GSM3017262.mat.gz"))

# DESCRIBE MATRICES -----------------------------------------------------------------------------------------
describe_matrices <- function(count_matrices){
  # Requires: count_matrices is a list of matrices; 
  # Returns: df is a data.frame with the following columns:
  #          names - names of the matrices in the list; 
  #          number_features - the number of rows/genes of each matrix;
  #          number_cells  - the number of columns/cells of each matrix.
  
  nr_columns <- c()
  nr_rows <- c()
  
  for (mtx in count_matrices){
    nr_columns <- c(nr_columns, dim(mtx)[2])
    nr_rows <- c(nr_rows, dim(mtx)[1])
  } #close for loop
  
  df <- data.frame(names = names(count_matrices),
                   number_features = nr_rows,
                   number_cells = nr_columns)
  return(df)
} #close describe_matrices

# IDENTIFY WORKFLOW -----------------------------------------------------------------------------------------
identify_workflow <- function(sample_path){
  # Requires: sample_path is a string with a file path to a GEO ID directory
  # Returns: a numeric from {1,2,3,4,5} that corresponds to the workflow needed to read the data
  #          into a list of count matrices.
  
  all_files <- list.files(path = sample_path ,pattern=".+mtx")
  
  # if we have mtx files, use workflow 2 
  if (sum(length(all_files),1)  > 1){ return(2)} 
  
  # if we have h5 files, use workflow 3
  else if (sum(length(list.files(path = sample_path ,pattern="h5")),1) > 1){return(3)}
  
  # if we have xlsx files, use workflow 4
  #else if (sum(length(list.files(path = sample_path ,pattern="xlsx")),1) > 1){return(4)}
  
  # otherwise, assume count matrix stored in a csv/tsv/txt/xlsx file and use workflow 1
  else{return(1)}
} #close identify_workflow

# MERGE MATRICES -----------------------------------------------------------------------------------------
merge_matrices <- function(countMatrices){
  # Requires: countMatrices is a list of count matrices
  # Returns: one count matrix which is the result of merging all count matrices 
  # This merge keeps every column of every matrix (cells) and every row (genes)
  
  print("Merging matrices...")
  
  # use the first count matrix as a seed for merging
  start <- countMatrices[[1]]
  countMatrices[[1]] <- NULL
  
  rownames(start) <- gsub("\\.", "-", rownames(start))
  
  cell_names <- c(colnames(start))
  
  for (mtx in countMatrices){
    rownames(mtx) <- gsub("\\.", "-", rownames(mtx))
    cell_names <-  c(cell_names, colnames(mtx))
    start <- merge(start, mtx, by = "row.names", all.x = TRUE, all.y = TRUE) #keep all genes
    rownames(start) <- start[,1]
    start <- start[,-1]
  } #close for loop
  
  # replace NAs with 0's
  start[is.na(start)] <- 0
  
  print("Merge of matrices complete.")
  
  return(list(start, cell_names))
} #close merge_matrices

# MERGE MATRICES BY ROW -----------------------------------------------------------------------------------------
merge_matrices_by_row <- function(metadata){
  # Requires: metadata a list of data frames
  # Returns: one data frame which is the result of merging all data frames
  # This merge assumes all data frames have same columns (variables) and adds only rows (cells)
  # It is intended to be use within dataset
  
  print("Merging metadata...")
  
  start <- metadata[[1]]
  metadata[[1]] <- NULL
  
  for (mtx in metadata){
    start <- rbind(start, mtx)
  } #close for loop
  
  print("Merging metadata complete.")
  
  return(start)
} #close merge_matrices_by_row

# get_panglaodb -----------------------------------------------------------------------------------------
get_panglaodb <- function(panglao_id){
  # Requires: panglao_id is a string of an SRA PanglaoDB ID, e.g., SRA553822
  # Returns: a list for which the first entry is a list of count matrices and the second a list of dataframes
  # with metadata information
  
  srts <- getSamples(sra = panglao_id, merge = FALSE)
  
  countMatrices <- list()
  metadata <- data.frame()

  for (srt in names(srts)){
    countMatrices[[srt]] <- as.matrix(Seurat::GetAssayData(srts[[srt]], slot = "counts"))
    metadata <- rbind(metadata, srts[[srt]][[]])
  } #close for loop
  
  metadata$dataset <- rep(panglao_id, dim(metadata)[1])
  cell_names <- c()
  
  for (mat in countMatrices){
    cell_names <- c(cell_names, colnames(mat))
  } #close for loop
  
  metadata$barcodes <- cell_names
  
  return(list(countMatrices, metadata))
} #close get_panglaodb

# SUBSET MATRICES COL -----------------------------------------------------------------------------------------
subset_matrices_col <- function(matrices, filter){
  for (mat in names(matrices)){
    matrices[[mat]] <- matrices[[mat]][,filter]}
}

# SUBSET MATRICES ROW -----------------------------------------------------------------------------------------
subset_matrices_row <- function(matrices, filter){
  for (mat in names(matrices)){
    matrices[[mat]] <- matrices[[mat]][filter,]}
}

# MERGE METADATA -----------------------------------------------------------------------------------------
merge_metadata <- function(meta1, meta2){
  # Requires: meta1 is a data.frame; meta2 is a data.frame; 
  # Ensures: if meta1 is not empty, returns a data.frame of meta1 + meta2; else returns meta2; 
  #          the merging is performed keeping all different columns and all different rows
  if (dim(meta1)[1] != 0){return(merge(meta1, meta2, all = TRUE))}
  else return(meta2)
  
}

# UNTAR EVERYTHING -----------------------------------------------------------------------------------------
untar_everything <- function(geo_id, folder_path){
  # Requires: geo_id is string; folder_path is a folder directory to a GEO dataset
  # Ensures: keeps unzipping .tar files contained inside each other; labels unzipped files with 
  #          folder/sample information to not lose sample identity 
  
  # First untar to unpack main directory
  files <- list.files(path = folder_path) 
  print("First round of files:")
  print(files) 
  
  for (file in files){
    if (sum(grep("tar", file), 1)  > 1){ 
      file_path <- paste(getwd(),"session_data",geo_id, file, sep = "/")
      print(file_path)
      untar(file_path,
            #extras = "--strip-components 1",
            exdir = folder_path)
      unlink(file_path)
    } #close if
  } # close for 
  
  # Second untar to resolve potential secondary directories
  files <- list.files(path = folder_path)
  print("Second round of files:")
  print(files)
  for (file in files){
    if (sum(grep("tar", file), 1)  > 1){
      file_path <- paste(getwd(),"session_data",geo_id, file, sep = "/")
      print("identified tar file:")
      print(file_path)
      untar(file_path, exdir = folder_path)
      folder <- list.dirs(folder_path)[2]
      print("created new folder:")
      print(folder)
      new_files <- list.files(path = folder)
      print("files inside folder:")
      print(new_files)
      
      for (new_file in new_files){
        print("changing name of file:")
        
        new_file_path <- paste(folder,new_file, sep = "/")
        print(new_file_path)
        new_file_name <- gsub("\\.tar\\.gz", "", file_path)
        new_file_name <- gsub("\\+", "_", new_file_name)
        new_file_name <- paste(new_file_name,new_file, sep ="_")
        print("to:")
        print(new_file_name)
        file.rename(new_file_path, new_file_name)
       
      } #close for
      unlink(folder, recursive = TRUE)
      unlink(file_path)
    } #close if 
  } # close for
  
  
}

#geo_id <- "GSE167597"
#getGEOSuppFiles(geo_id, baseDir = "./session_data") 
#folder_path <- "/mnt/scratch/home/marta/shiny/scStudio_2.0/session_data/GSE167597"
#untar_everything(geo_id, folder_path)

#meta1 <- vroom("/mnt/scratch/home/marta/shiny/scStudio_2.0/tokens/dgw7bm/testing_metadata.txt")
#test <- merge_metadata(meta1,meta2)

load_sce <- function(sce_path, name){
  
  sce <- readRDS(sce_path$datapath)
  counts <- counts(sce)
  metadata <- as.data.frame(colData(sce))
  metadata$orig.ident <- name
  metadata$barcodes <- colnames(counts)
  
  metadata$dataset <- colnames(counts)  
  
  cell_ids <- define_cell_ids(ncol(counts))

  metadata$cell_id <- cell_ids
  colnames(counts) <- cell_ids
  
  return(list(counts, metadata))
}


load_srt <- function(srt_path, name){
  srt <- readRDS(srt_path$datapath)
  counts <- GetAssayData(srt, slot = "counts")
  metadata <- srt@meta.data
  metadata$orig.ident <- name
  metadata$barcodes <- colnames(counts)
  
  metadata$dataset <- colnames(counts)  
  
  cell_ids <- define_cell_ids(ncol(counts))
  
  metadata$cell_id <- cell_ids
  colnames(counts) <- cell_ids
  
  return(list(counts, metadata))
}

load_10X_h5 <- function(h5_path, name){
  
  mat <- Read10X_h5(h5_path$datapath, 
             use.names = TRUE, 
             unique.features = TRUE)
  
  cell_ids <- define_cell_ids(ncol(mat))
  metadata <- data.frame(cell_ids = cell_ids,
                         barcodes = colnames(mat))
  
  metadata$orig.ident <- name
  metadata$dataset <- name
  
  colnames(mat) <- cell_ids
  
  return(list(mat, metadata))
  
}

load_10X_triplet <- function(barcodes_path, features_path, matrix_path, name){
  
  mat <- readMM(file = matrix_path$datapath)
  
  feature.names = read.delim(features_path$datapath, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  barcode.names = read.delim(barcodes_path$datapath, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V2
  feature.ids = feature.names$V1
  
  cell_ids <- define_cell_ids(ncol(mat))
  metadata <- data.frame(cell_ids = cell_ids,
                         barcodes = colnames(mat))
  
  metadata$orig.ident <- name
  metadata$dataset <- name
  
  colnames(mat) <- cell_ids
  
  return(list(mat, metadata))
}


define_cell_ids <- function(number_cells){
  cell_ids <- c()
  # for each cell, define a unique key
  for (cell in 1:number_cells){
    key <- gsub("[-:WET]","",Sys.time())
    key <- gsub(" ","", key)
    cell_ids <- c(cell_ids, paste0(key, stri_rand_strings(1, 6, pattern = "[a-z0-9]")))
  }
  return(cell_ids)
}








