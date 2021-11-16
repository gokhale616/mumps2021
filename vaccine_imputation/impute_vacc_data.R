source("../00/src.R", chdir = TRUE)
source("./process_tycho_data.R")
source("./VSEI_1_step.R")
source("./filtering_utils.R")


# run the particle filter to impute missing vaccine values and save the the results in ../processed_data 
# takes roughly 3 hours to run


filtered_vacc <- pfilter_once(j_particles = 1000)


save(filtered_vacc, file = "../processed_data/filtered_vacc.rds")






