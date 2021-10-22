source("./treat_vacc_covar.R", chdir = TRUE)


# experiment with covariates
# log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
#                       t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)
# 
# mod_mumps_covariates_exptl <- treat_covar_data(par = log_par_defaults) 
# 
# plot_covars(covar_data = mod_mumps_covariates_exptl, filter_from = 1965)



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

common_pnames <- c(sprintf("q_age_%d", 1:5), sprintf("rho_age_%d", 1:5), "rho_age_u", 
                   sprintf("psi_%d", 1:5), "psi_u")

################# these are some utility functions to make the plots as we want theme to be ##################

mod_param_vals_est <- param_vals_est

mod_param_vals_est[c(common_pnames, "dwan")] <- (
  c(0.188365, 0.939224, 0.246510, 0.597938, 0.238201, 0.005039, 
    0.004913, 0.032407, 0.011305, 0.067853, 0.008251, 0.936209, 
    0.894261, 1.376325, 1.102456, 0.717972, 0.875927, 
    361.021130) %.>% 
    signif(., digits = 2)
  )
  


# this function allows to be looped over a range of parameter values 
multi_traj <- function(counter, p_name, p_vector, sim_p_vals = param_vals_est, 
                       sim_cases = FALSE, nsim = 1000) {
  
  default_p_val <- sim_p_vals[p_name]
  sim_p_vals[p_name] <- p_vector[counter]
  # browser()
  
  if(sim_cases == TRUE){
    
    sim_df <- (
      po_est_test %.>% 
        sim_obs_model(po_obj = ., params = sim_p_vals, 
                      times = time(po_est_test),
                      nsim = nsim)
      )
  } else {
      sim_df <- (
        po_est_test %.>% 
          trajectory(., sim_p_vals, format = "d", method = "ode23") %.>%
          select(., -`.id`)
        ) 
  }
  
  sim_df %.>% 
    mutate(., 
           Value = p_vector[counter])
  
}  

# generate data frames looped over parameter values 
vlen <- 10


tic()
moving_dwan <- map_dfr(1:vlen, multi_traj, p_name = "dwan", 
                        p_vector = seq(20, 400, length.out = 10),
                       sim_p_vals = mod_param_vals_est)
toc()

plot_sweep_dynamics(data = moving_dwan, param_sweep = TRUE, 
                    y_lab = "True Cases (C)",
                    colour_lab = expression(paste("Vaccine derived\nimmune duration (", 1/delta,", Years)"))
                    )


# test the simulator
# sim_obs_model(po_obj = po_est_test, params = param_vals_est, times = time(po_est_test), nsim = 5)
# it works!! 

# this function loops over parameter values and simulates from  the observation process
tic()
moving_dwan_sim <- map_dfr(1:vlen, multi_traj, sim_cases = TRUE, p_name = "dwan", 
                           p_vector = seq(20, 400, length.out = 10),
                           sim_p_vals = mod_param_vals_est)
toc()


moving_dwan_sim %.>%
  mutate_at(., 
            vars(`0.025`, `0.5`, `0.975`), 
            sqrt) %.>% 
  plot_sweep_dynamics(data = ., param_sweep = TRUE, 
                      sim_cases = TRUE,
                      y_lab = expression(sqrt(Simulated~Cases)),
                      colour_lab = expression(paste("Vaccine derived\nimmune duration (", 1/delta,", Years)"))
                      )






