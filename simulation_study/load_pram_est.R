# load all of the result objects in the global environment 
mle_result_path <- "../result_data/mle"

list.files(path = mle_result_path, 
           full.names = TRUE)[-25] %.>% 
  lapply(., load,  envir = .GlobalEnv)


# One of the waning hypothesis was wrongly named to noloss
# correct the hypothesis name 
mle_waning_constant_is$Hypothesis <- "waning_constant"

# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_gamma_waning_slow_is, mle_gamma_waning_sigmoid_is, 
    mle_gamma_waning_rapid_is, mle_gamma_waning_constant_is,
    mle_waning_slow_is, mle_waning_sigmoid_is, mle_waning_rapid_is, mle_waning_constant_is,
    mle_leaky2_2_slow_is, mle_leaky2_2_sigmoid_is, mle_leaky2_2_rapid_is, mle_leaky2_2_constant_is,
    mle_leaky2_3_slow_is, mle_leaky2_3_sigmoid_is, mle_leaky2_3_rapid_is, mle_leaky2_3_constant_is,
    mle_leaky2_4_slow_is, mle_leaky2_4_sigmoid_is, mle_leaky2_4_rapid_is, mle_leaky2_4_constant_is, 
    mle_slow_is, mle_sigmoid_is, mle_rapid_is, mle_constant_is
  )
)

# hard coding some parameter values for use downstream within estimation of R0s and Rps
# calculate the reproductive numbers 
params_for_R0 <- c(N = 100e6, nu = 1/80, p = 0, ad = age_class_duration)
params_for_Rp <- c(N = 100e6, nu = 1/80, p = 1, ad = age_class_duration)


mk_result_df <- function(c = 1, res = result_list) {
  # to look into the function at a specific iteration
  # if (c == 20) {browser()}
  
  # collect qualitative covariates
  extra_params <- res[[c]]$Hypothesis %.>% str_split(., pattern = "_") %.>% unlist(.)
  
  if(length(extra_params) == 3) {
    hypo_covar <- extra_params[1]
    p_intro    <- extra_params[2] %.>% as.numeric(.)
    vacc       <- extra_params[3]
  } else {
    hypo_covar <- extra_params[1]
    p_intro    <- NA
    vacc       <- extra_params[2]
  }
  
  if(hypo_covar == "gwaning") {
    R0 <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_R0) %.>% 
      calculate_R0_mq(.)$reprodutive_number 
    
    Rp <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_Rp) %.>% 
      calculate_R0_mq(.)$reprodutive_number   
  } else {
    R0 <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_R0) %.>% 
      calculate_R0_mq(.)$reprodutive_number 
    
    Rp <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_Rp) %.>% 
      calculate_R0_mq(.)$reprodutive_number   
  }
  
  
  # collate results in a dataframe
  res[[c]]$DEobj$optim$bestmem %.>% 
    as.list(.) %.>% 
    as_tibble(.) %>% 
    mutate(., 
           R0 = R0,
           Rp = Rp,
           impact = 1-Rp/R0,
           loglik = -res[[c]]$DEobj$optim$bestval, 
           npar = res[[c]]$DEobj$optim$bestmem %.>% length(.), 
           AIC = calculate_aic(loglik, npar), 
           hypothesis = hypo_covar, 
           p_intro    = p_intro,
           vacc_covariate = vacc) %.>% 
    select(., -c(loglik, npar))
  
}


all_result_df <- (
  map_dfr(1:length(result_list), mk_result_df) %.>% 
    mutate(., 
           d_AIC = AIC - min(AIC)) %.>% 
    mutate(., 
           best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC)
)
