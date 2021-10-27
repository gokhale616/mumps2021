# load all of the result objects in the global environment 
mle_result_path <- "../result_data/mle/"

list.files(path = mle_result_path, 
           full.names = TRUE) %.>% 
  lapply(., load,  envir = .GlobalEnv)



# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_waning_slow, mle_waning_slow_sigmoid, mle_waning_slow_rapid, mle_waning_slow_constant
  )
)


mk_result_df <- function(c = 1, res = result_list) {
  
  hypo_covar <- res[[c]]$Hypothesis %.>% str_split(., pattern = "_") %.>% unlist(.)
  
  # add an empty space if the lengths are a mismatch
  if(length(hypo_covar) == 2) {
    hypo_covar <- c(hypo_covar[1:2], "")
  } 
  

  # collate results in a dataframe
  res[[c]]$DEobj$optim$bestmem %.>% 
    as.list(.) %.>% 
    as_tibble(.) %>% 
    mutate(., 
           loglik = -res[[c]]$DEobj$optim$bestval, 
           npar = res[[c]]$DEobj$optim$bestmem %.>% length(.), 
           AIC = calculate_aic(loglik, npar), 
           hypothesis = paste0(hypo_covar[1]), 
           vacc_covariate = paste0(hypo_covar[2], hypo_covar[3])) %.>% 
    select(., -c(loglik, npar))

}


all_result_df <- (
  map_dfr(1:length(result_list), mk_result_df) %.>% 
    mutate(., 
           d_AIC = AIC - min(AIC)) %.>% 
    group_by(., hypothesis) %.>% 
    mutate(., best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC)
)


# table of estimates 
all_result_df %.>% 
  select(., -best_fit_covar) %.>% 
  gather(., 
         key = "Parameter", value = "Estimate", -c(hypothesis, vacc_covariate, d_AIC)) %.>% 
  arrange(.)


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

po_est_test <- make_pomp(covar = mod_mumps_covariates_slow_sigmoid)


obs_sim <- sim_from_these %.>%   
  sim_obs_model(po_est_test, params = ., times = time(po_est_test), nsim = 1e3) %.>% 
  mutate(., hypothesis = "Waning")

  
obs_data <- mumps_case_reports %.>%  
  select(., -total) %.>% 
  gather(., key = "age_class", value = `0.5`, -year, factor_key = TRUE) %.>% 
  mutate(., 
         `0.5` = sqrt(`0.5`))
  
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
  




