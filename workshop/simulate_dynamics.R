source("../00/src.R", chdir = TRUE)

# check if the model compiles 
# First check if the model compiles with data 
po_est_test <- make_pomp(covar = mumps_covariates)

# now check if the model compiles with mock data
po_sim_test <- make_pomp(covar = mumps_covariates, extrapolate_simulation = TRUE)

