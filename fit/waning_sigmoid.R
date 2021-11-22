# source the pre-reqs
source("./fit_prereqs.R")


mle_waning_sigmoid <- DE_traj_match(covar = mod_mumps_covariates_sigmoidal, 
                                         param_constraints = param_range_waning, 
                                         params = param_vals_est, 
                                         ode_control = list(method = "ode23"),
                                         hypo_name = "waning_sigmoid")



save(mle_waning_sigmoid, file = "../result_data/mle/mle_waning_sigmoid.rds")
