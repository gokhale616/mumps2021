# read the processed data for analysis
path_processed_data <- "../processed_data/"

list.files(path = path_processed_data, 
           full.names = TRUE) %.>% 
  lapply(., 
         load, envir = .GlobalEnv)