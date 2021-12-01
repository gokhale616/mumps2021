# Assumes the use of processed data and no plots
lib_dep_cl <- c("lubridate", "tidyverse", 
                "pomp", "magrittr", "wrapr",
                "tictoc", "parallel", "subplex", 
                "LaplacesDemon",
                "DEoptim", "doParallel", "doRNG",
                "nloptr")

sapply(lib_dep_cl, library, character.only = TRUE)




# definitions some useful objects that can be used to quickly 


age_names <- c("[0,5)", "[5,15)", "[15,25)", "[25,40)", ">40") 

age_names_u <- c(age_names, "unknown") 


# reluctantly placing this function here - might not be the best thing!
age_class_lengths <- function(age_class = c("[00,05)")) {
  if(str_detect(age_class, ">") == TRUE) {
    return(80 - as.numeric(str_sub(age_class, start = 2L, end = -1L)))
  } else {
    return(as.numeric(str_sub(age_class, start = 5L, end = 6L)) - 
             as.numeric(str_sub(age_class, start = 2L, end = 3L)))
  }
}

age_class_duration <-  map_dbl(1:length(age_names), 
                               function(c, an = c("[00,05)", "[05,15)", "[15,25)", "[25,40)", ">40")) {
                                 age_class_lengths(age_class = an[c])
                               })


c_levels = c(sprintf("S_%d",1:5),  
             sprintf("E1_%d",1:5), sprintf("I1_%d",1:5),  
             sprintf("E2_%d",1:5), sprintf("I2_%d",1:5),
             sprintf("R_%d",1:5), sprintf("V_%d",1:5) , 
             sprintf("C_%d",1:5), sprintf("D_%d",1:5),
             sprintf("A_%d",1:5), sprintf("G_%d",1:5),
             sprintf("N_%d",1:5))



# percentile values
percentile_vals <- c(0.025, 0.5, 0.975)


# set up controls for deoptim
### controls list for DEoptim - for MLE estimation
np_val = 500
my_controls <- list(itermax = 1e5,
                    F = 0.6, CR = 0.9, 
                    strategy = 1,
                    steptol = 750, reltol = 1e-4)

