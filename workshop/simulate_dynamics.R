source("./treat_vacc_covar.R", chdir = TRUE)


# experiment with covariates
log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
                      t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)

mod_mumps_covariates_exptl <- treat_covar_data(par = log_par_defaults) 

plot_covars(covar_data = mod_mumps_covariates_exptl, filter_from = 1965)



# check if the model compiles 
# First check if the model compiles with data 
po_est_test <- make_pomp(covar = mod_mumps_covariates_slow_constant)
plot_covars(mod_mumps_covariates_slow_constant)

# # now check if the model compiles with mock data
# po_sim_test <- make_pomp(covar = mumps_covariates, extrapolate_simulation = TRUE)
# 
# 
# # use these pomp objects for a quick check
# est_traj <- po_est_test %.>% 
#   trajectory(., param = param_vals_est, format = "d", method = "ode45")
#   
# est_traj %.>% plot_dynamics(.)


################# these are some utility functions to make the plots as we want theme to be ##################

# this function allows to be looped over a range of parameter values 
multi_traj <- function(counter, p_name, p_vector, sim_p_vals = param_vals_est) {
  
  default_p_val <- sim_p_vals[p_name]
  sim_p_vals[p_name] <- p_vector[counter]
  # browser()
  
  po_est_test %.>% 
    trajectory(., sim_p_vals, format = "d", method = "ode23") %.>%
    # select(., -`.id`) %>%
    mutate(., 
           Value = p_vector[counter])
}  

# generate data frames looped over parameter values 
vlen <- 10


tic()
moving_dwan <- map_dfr(1:vlen, multi_traj, p_name = "dwan", 
                        p_vector = c(50000, seq(10, 500, length.out = 9))
                        )
toc()

plot_sweep_dynamics(data = moving_dwan, param_sweep = TRUE, lab_val = expression(1/delta))

