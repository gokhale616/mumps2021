source("../00/src.R", chdir = TRUE)

# check if the model compiles 
# First check if the model compiles with data 
po_est_test <- make_pomp(covar = mumps_covariates)

# now check if the model compiles with mock data
po_sim_test <- make_pomp(covar = mumps_covariates, extrapolate_simulation = TRUE)


# use these pomp objects for a quick check
est_traj <- po_est_test %.>% 
  trajectory(., param = param_vals_est, format = "d", method = "ode45")
  
est_traj %.>% plot_dynamics(.)





