
# source the pre-reqs
source("./fit_prereqs.R")

# rest param for intro in the right age class
param_vals_est_l2_4 <- set_param(new_val = 3)

mle_leaky2_4_rapid <- DE_traj_match(covar = mod_mumps_covariates_rapid, 
                                      param_constraints = param_range_leaky2, 
                                      params = param_vals_est_l2_4, 
                                      ode_control = list(method = "ode23"),
                                      hypo_name = "leaky2_4_rapid")



save(mle_leaky2_4_rapid, file = "../result_data/mle/mle_leaky2_4_rapid.rds")
