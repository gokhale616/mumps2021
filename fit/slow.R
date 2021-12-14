# source the pre-reqs
source("./fit_prereqs.R")


# estimate parameters for the in sample (is) data (1977-2012)
mle_slow_is <- DE_traj_match(covar = mod_mumps_covariates_slow, 
                                 param_constraints = hypo_ind_params, 
                                 params = param_vals_est, 
                                 ode_control = list(method = "ode23"),
                                 hypo_name = "noloss_slow", 
                                 incidence_data = in_sample_mumps_case_reports)



save(mle_slow_is, file = "../result_data/mle/mle_slow_is.rds")



