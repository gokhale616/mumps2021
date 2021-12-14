# source the pre-reqs
source("./fit_prereqs.R")

# bring in a good guess (gg) from the past 
load("../result_data/mle/old/mle_waning_slow.rds")

# assign to a dataframe
gg_waning_slow <- mle_waning_slow$DEobj$optim$bestmem %.>% as.list(.) %.>% as_tibble(.)

# estimate parameters for the in sample (is) data (1977-2012)
mle_waning_slow_is <- DE_traj_match(covar = mod_mumps_covariates_slow, 
                                    param_constraints = param_range_waning, 
                                    params = param_vals_est, 
                                    ode_control = list(method = "ode23"),
                                    hypo_name = "waning_slow", 
                                    best_past_est = gg_waning_slow, 
                                    incidence_data = in_sample_mumps_case_reports)



save(mle_waning_slow_is, file = "../result_data/mle/mle_waning_slow_is.rds")
