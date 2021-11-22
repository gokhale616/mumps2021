# source the pre-reqs
source("./fit_prereqs.R")


mle_waning_constant <- DE_traj_match(covar = mod_mumps_covariates_constant, 
                                          param_constraints = param_range_waning, 
                                          params = param_vals_est, 
                                          ode_control = list(method = "ode23"),
                                          hypo_name = "waning_constant")



save(mle_waning_constant, file = "../result_data/mle/mle_waning_consant.rds")



