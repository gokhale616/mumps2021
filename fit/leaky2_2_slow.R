
# source the pre-reqs
source("./fit_prereqs.R")

# rest param for intro in the right age class
param_vals_est_l2_2 <- set_param(new_val = 1)

# bring in a good guess (gg) from the past 
load("../result_data/mle/old/mle_leaky2_2_slow.rds")

# assign to a dataframe
gg_leaky2_2_slow <- mle_leaky2_2_slow$DEobj$optim$bestmem %.>% as.list(.) %.>% as_tibble(.)

# estimate parameters for the in sample (is) data (1977-2012)
mle_leaky2_2_slow_is <- DE_traj_match(covar = mod_mumps_covariates_slow, 
                                      param_constraints = param_range_leaky2, 
                                      params = param_vals_est_l2_2, 
                                      ode_control = list(method = "ode23"),
                                      hypo_name = "leaky2_2_slow", 
                                      best_past_est = gg_leaky2_2_slow, 
                                      incidence_data = in_sample_mumps_case_reports)



save(mle_leaky2_2_slow_is, file = "../result_data/mle/mle_leaky2_2_slow_is.rds")
