# Load packages ----------------------------------------------------------------------------------------------

lib_dep <- c("tidyverse", "pomp", "magrittr", "wrapr",
             "readxl", "tictoc", "parallel", "subplex", 
             "LaplacesDemon", "cowplot", 
             "rootSolve", "GGally", "RColorBrewer", 
             "ggpubr", "reshape2",
             "readxl", "socialmixr", "ggthemes", 
             "GA", "DEoptim", "doParallel", "doRNG", "xtable", "scales",
             "parallel", "nloptr")

sapply(lib_dep, require, character.only = TRUE)


# read the processed data for analysis
load("../processed_data/mumps_case_reports.rds")
load("../processed_data/mumps_covariates.rds")
load("../processed_data/contact_matrix.rds")


# Set the project theme and the base font size of all labels 
# definintion some useful objects that can be used to quickly 


age_names <- c("[0,5)", "[5,15)", "[15,25)", "[25,40)", "<40") 

age_names_u <- c(age_names, "unknown") 


## PNAS specs
## For details see https://blog.pnas.org/digitalart.pdf
## max height = 9in / 22.5cm
## knitr expects inches by default
## scale up by *1.56 here, down by 0.64 in final tex
#figwidth <- c(8.7, 11.4, 17.8)/2.54
figwidth <- c(13.6, 17.8, 27.7)/2.54
## Intended target dims
outwidth <- c('8.7cm', '11.4cm', '17.8cm')


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


# plots related definitions #
project_theme <- (theme_classic2(base_size = 13) +
                  # background_grid() +
                  theme(strip.background = element_blank(),
                        legend.position = "bottom",
                        legend.background = element_blank(), 
                        legend.key = element_blank(),
                        axis.ticks.length=unit(.25, "cm"))
                  )
#darkslateblue
#ivory3
# springgreen4
# age class colors
age_cols <- c(brewer_pal(palette = "Purples")(5)[5:1], "steelblue4") %.>% 
  setNames(., age_names_u)

# colour schemes for models and data 
mod_colours <- c("Black", "Green", "Purple", "Red")









