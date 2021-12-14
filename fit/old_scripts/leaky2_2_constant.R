
# source the pre-reqs
source("./fit_prereqs.R")

# rest param for intro in the right age class
param_vals_est_l2_2 <- set_param(new_val = 1)

mle_leaky2_2_constant <- DE_traj_match(covar = mod_mumps_covariates_constant, 
                                       param_constraints = param_range_leaky2, 
                                       params = param_vals_est_l2_2, 
                                       ode_control = list(method = "ode23"),
                                       hypo_name = "leaky2_2_constant")



save(mle_leaky2_2_constant, file = "../result_data/mle/mle_leaky2_2_constant.rds")
