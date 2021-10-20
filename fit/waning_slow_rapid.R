# source the pre-reqs
source("./fit_prereqs.R")


mle_waning_slow_rapid <- DE_traj_match(covar = mod_mumps_covariates_slow_rapid, 
                                       param_constraints = param_range_waning, 
                                       params = param_vals_est, 
                                       ode_control = list(method = "ode23"),
                                       hypo_name = "waning_slow_rapid")



save(mle_waning_slow_rapid, file = "../result_data/mle/mle_waning_slow_rapid.rds")






