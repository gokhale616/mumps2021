# source the pre-reqs
source("./fit_prereqs.R")


mle_waning_slow <- DE_traj_match(covar = mod_mumps_covariates_slow, 
                                 param_constraints = param_range_waning, 
                                 params = param_vals_est, 
                                 ode_control = list(method = "ode23"),
                                 hypo_name = "waning_slow")



save(mle_waning_slow, file = "../result_data/mle/mle_waning_slow.rds")
