# load all of the result objects in the global environment 
mle_result_path <- "../result_data/mle/"

list.files(path = mle_result_path, 
           full.names = TRUE) %.>% 
  lapply(., load,  envir = .GlobalEnv)



# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_waning_slow, mle_waning_sigmoid, mle_waning_rapid, mle_waning_constant,
    mle_leaky2_2_slow, mle_leaky2_2_sigmoid, mle_leaky2_2_rapid, mle_leaky2_2_constant,
    mle_leaky2_3_slow, mle_leaky2_3_sigmoid, mle_leaky2_3_rapid, mle_leaky2_3_constant,
    mle_leaky2_4_slow, mle_leaky2_4_sigmoid, mle_leaky2_4_rapid, mle_leaky2_4_constant
  )
)


mk_result_df <- function(c = 1, res = result_list) {
  
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
  
  # calculate the reproductive numbers 
  params_for_R0 <- c(N = 100e6, nu = 1/80, p = 0, ad = age_class_duration)
  params_for_Rp <- c(N = 100e6, nu = 1/80, p = 0.5, ad = age_class_duration)
  
  R0 <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_R0) %.>% 
    calculate_R0_mq(.)$reprodutive_number 
  
  Rp <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_Rp) %.>% 
    calculate_R0_mq(.)$reprodutive_number 
  
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
    group_by(., hypothesis) %.>% 
    mutate(., 
           best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC)
)






# table of estimates - process model 
#table_hypo_compare <- ()

all_result_df %.>% 
  select(., -c(best_fit_covar, starts_with("q_age"), 
               starts_with("rho_age"), starts_with("psi"), 
               p_intro)) %.>% 
  arrange(., d_AIC)
  group_by(., hypothesis) %.>% 
  filter(., d_AIC == min(d_AIC)) %.>% 
  ungroup(.) %.>%
  mutate_if(., is.numeric, function(x){round(x, digits = 3) %.>% format(., nsmall = 3)})  
  mutate(., 
         vacc_covariate = str_to_title(vacc_covariate),
         hypothesis = str_to_title(hypothesis)) %.>% 
  gather(., 
         key = "Parameter", value = "Estimate", 
         -c(hypothesis)) %.>% 
  spread(., key = hypothesis, value = Estimate) 
  


sim_from_these <- (
  all_result_df %.>%   
    filter(., best_fit_covar == TRUE) %.>% 
    select(., 
           -c(hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
    unlist(.) %.>% 
    sim_p_vals(.)
    
)

# make this happen nicely  
source("../workshop/treat_vacc_covar.R", chdir = TRUE)

# generate a po_object

po_est_test <- make_pomp(covar = mod_mumps_covariates_sigmoidal)


obs_sim <- sim_from_these %.>%   
  sim_obs_model(po_est_test, params = ., times = time(po_est_test), nsim = 1e3) %.>% 
  mutate(., hypothesis = "Waning")

  
obs_data <- mumps_case_reports %.>%  
  select(., -total) %.>% 
  gather(., key = "age_class", value = `0.5`, -year, factor_key = TRUE) %.>% 
  mutate(., 
         `0.5` = sqrt(`0.5`))
  
fit_plot <- (
  obs_sim %.>%   
  ggplot(., aes(x = year)) +
  geom_line(aes(y = `0.5`, colour = hypothesis), size = 0.8) +
  geom_line(data = obs_data, 
            aes(y = `0.5`), size = 0.8) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`, fill = hypothesis), alpha= 0.6) +
  labs(y = expression(sqrt(Cases)), 
       x = "Year", 
       fill = "Hypothesis", 
       colour = "Hypothesis") +
  facet_grid(rows = vars(age_class), cols = vars(hypothesis), 
             scales = "fixed") +
  scale_fill_manual(values = "#45a247") +
  scale_colour_manual(values = "#45a247") +
  scale_x_continuous(breaks = gen_x_breaks, expand = c(0, 0)) +
  project_theme +
  cap_axes
  )
  




