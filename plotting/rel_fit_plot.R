# this script depends on ./prep_estm_tables.R

# decide on model colors 
# best fitting model will always be - "#11998e"
# scales::show_col(c("#11998e", "#003973", "#FF416C", "#654ea3"))


# generate a tibble of parameter vectors
sim_from_these_tibble <- (
  all_result_df %.>%   
    group_by(., hypothesis) %.>% 
    filter(., d_AIC == min(d_AIC)) %.>% 
    mutate(., 
           p_intro  = ifelse(is.na(p_intro) == TRUE, 6, p_intro), 
           t_intro  = ifelse(is.na(p_intro) == TRUE, 3000, p_intro), 
           dwan     = ifelse(is.na(dwan) == TRUE, Inf, dwan),  
           epsilon2 = ifelse(is.na(epsilon2) == TRUE, 0, epsilon2), 
           hypothesis = case_when(hypothesis == "waning"~ "Waning", 
                                  hypothesis == "leaky2"~ "Leaky", 
                                  hypothesis == "noloss"~ "No Loss")) %.>% 
    ungroup(.)
  
)


# generate the right time period of data for the relative plots (1977-2012)
source("../fit/prep_inc_data.R", chdir = TRUE)

# generate simulated dynamics for all three models 
sim_all_models <- map_dfr(1:nrow(sim_from_these_tibble), 
                          function(c, res = sim_from_these_tibble) {
                            # browser()
                            sliced_res <- res %.>% slice(., c)
                            
                            # collect some information for later  
                            hypo_name <-  sliced_res %.>% select(., hypothesis) %.>% unlist(.)
                            covar_name <- sliced_res %.>% select(., vacc_covariate) %.>% unlist(.)
                            
                            # generate the appropriate pomp  object 
                            if(covar_name == "constant") {
                              mod_po <- make_pomp(covar = mod_mumps_covariates_constant, 
                                                  incidence_data = in_sample_mumps_case_reports)
                            } else if(covar_name == "rapid") {
                                mod_po <- make_pomp(covar = mod_mumps_covariates_rapid, 
                                                    incidence_data = in_sample_mumps_case_reports)
                            } else if(covar_name == "sigmoid") {
                                mod_po <- make_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                                    incidence_data = in_sample_mumps_case_reports)
                            } else if(covar_name == "slow") {
                                mod_po <- make_pomp(covar = mod_mumps_covariates_slow, 
                                                    incidence_data = in_sample_mumps_case_reports)
                            }
                            
                            
                            # produce a vector of parameters to simulate from
                            sim_params <- (
                              sliced_res %.>% 
                                select(., -c(hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
                                unlist(.) %.>% 
                                sim_p_vals(.)
                              )
                            
                            # simulate the observation model and add the hypothesis name 
                            sim_params %.>% 
                              sim_obs_model(mod_po, params = ., times = time(mod_po), nsim = 1e3, 
                                            quantile_prob = c(0.1, 0.5, 0.90), 
                                            quantile_prob_names = c("0.1", "0.5", "0.90")) %.>% 
                              mutate(., 
                                     hypothesis = hypo_name)
                            
                          })



# anno_data for delta AIC
anno_d_AIC <- (
  sim_from_these_tibble %.>% 
    select(., hypothesis, d_AIC) %.>% 
    mutate(., age_class = "[0,5)")
  )


# prep_obs_data
in_sample_obs_data <- (
  in_sample_mumps_case_reports %.>%  
    select(., -total) %.>% 
    gather(., key = "age_class", value = `0.5`, -year, factor_key = TRUE) %.>% 
    mutate(., 
           `0.5` = sqrt(`0.5`)) %.>% 
    drop_na(.)
  )

# relative fits to the data 
rel_fit_plots <- (
  sim_all_models %.>% 
    mutate(., hypothesis = as_factor(hypothesis)) %.>% 
    ggplot(.) +
    geom_line(data = in_sample_obs_data, aes(x = year, y = `0.5`, colour = "Data"), 
              size = 1) +
    geom_line(aes(x = year, y = `0.5`, colour = hypothesis), size = 0.8) +
    geom_ribbon(aes(x = year, ymin = `0.1`, ymax = `0.90`, fill = hypothesis), alpha = 0.6) +
    geom_text(data = anno_d_AIC, aes(x = 1985, y = 150, label = paste0("Delta~AIC == ", 
                                                                       round(d_AIC, 2) %.>% 
                                                                         as.character(.))), 
              parse = TRUE) +
    facet_grid(rows = vars(age_class), cols = vars(hypothesis), scales = "fixed") +
    labs(x = "Year", 
         y = expression(sqrt(Cases)), 
         colour = "Cases", 
         fill = "Cases") +
    scale_colour_manual(values = c("Data" = "grey30", "Waning" = "#11998e", 
                                   "Leaky" = "#FF416C", "No Loss" = "#654ea3")) +
    scale_fill_manual(values = c("Leaky" = "#FF416C","Waning" = "#11998e", "No Loss" = "#654ea3")) +
    project_theme +
    cap_axes
  )














