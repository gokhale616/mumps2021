# source the pre-reqs
source("./fit_prereqs.R")

# bring in a good guess (gg) from the past 
load("../result_data/mle/old/mle_waning_constant.rds")

# assign to a dataframe
gg_waning_constant <- mle_waning_constant$DEobj$optim$bestmem %.>% as.list(.) %.>% as_tibble(.)

# estimate parameters for the in sample (is) data (1977-2012)
mle_gamma_waning_constant_is <- DE_traj_match(covar = mod_mumps_covariates_constant, 
                                              waning_distribution = "gamma",
                                              param_constraints = param_range_waning, 
                                              params = param_vals_est, 
                                              ode_control = list(method = "ode23"),
                                              hypo_name = "gwaning_constant", 
                                              best_past_est = gg_waning_constant, 
                                              incidence_data = in_sample_mumps_case_reports)



save(mle_gamma_waning_constant_is, file = "../result_data/mle/mle_gamma_waning_constant_is.rds")



