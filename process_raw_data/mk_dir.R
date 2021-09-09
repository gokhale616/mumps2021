# this scripts produces a new directory to store processed data 
if(dir.create("../processed_data") == TRUE) {
  message("'../processed_data' already exists at root! Proceeeding to populate")
} else {
  dir.exists("../processed_data")  
  message("'../processed_data' created at root! Proceeeding to populate")
}

