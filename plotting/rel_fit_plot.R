# this script depends on ./prep_estm_tables.R

# decide on model colors 
# best fitting model will always be - "#11998e"
#scales::show_col(c("#11998e", "#003973", "#FF416C", "#654ea3"))


# generate a tibble of parameter vectors
sim_from_these_tibble <- (
   all_result_df %.>%   
    group_by(., hypothesis) %.>% 
    filter(., d_AIC == min(d_AIC)) %.>% 
    mutate(., 
           p_intro  = ifelse(is.na(p_intro) == TRUE, 6, p_intro-1), 
           t_intro  = ifelse(is.na(t_intro) == TRUE, 3000, t_intro), 
           dwan     = ifelse(is.na(dwan) == TRUE, Inf, dwan),  
           epsilon2 = ifelse(is.na(epsilon2) == TRUE, 0, epsilon2), 
           ) %.>% 
    ungroup(.)
  
)


# generate the right time period of data for the relative plots (1977-2012)
source("../fit/prep_inc_data.R", chdir = TRUE)

# generate simulated dynamics for all three models 
sim_all_models_pre <- (
  map_dfr(1:nrow(sim_from_these_tibble), 
                          function(c, res = sim_from_these_tibble) {
                            
                            sliced_res <- res %.>% slice(., c)
                            
                            # collect some information for later  
                            hypo_name <-  sliced_res %.>% select(., hypothesis) %.>% unlist(.)
                            covar_name <- sliced_res %.>% select(., vacc_covariate) %.>% unlist(.)
                            
                            # generate the appropriate pomp  object 
                            if(covar_name == "constant") {
                              if (hypo_name == "gwaning") {
                                mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_constant, 
                                                          incidence_data = in_sample_mumps_case_reports)
                              } else if (hypo_name == "gwaning2") {
                                mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_constant, 
                                                              incidence_data = in_sample_mumps_case_reports) 
                              } else {
                                  mod_po <- make_pomp(covar = mod_mumps_covariates_constant, 
                                                      incidence_data = in_sample_mumps_case_reports)
                              }
                            } else if(covar_name == "rapid") {
                              if (hypo_name == "gwaning") {
                                mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_rapid, 
                                                          incidence_data = in_sample_mumps_case_reports)
                              } else if (hypo_name == "gwaning2") {
                                mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_rapid, 
                                                              incidence_data = in_sample_mumps_case_reports)
                              } else {
                                  mod_po <- make_pomp(covar = mod_mumps_covariates_rapid, 
                                                      incidence_data = in_sample_mumps_case_reports)
                              }
                            } else if(covar_name == "sigmoid") {
                              if (hypo_name == "gwaning") {
                                mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                                          incidence_data = in_sample_mumps_case_reports)
                              } else if (hypo_name == "gwaning2") {
                                mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                                              incidence_data = in_sample_mumps_case_reports)  
                              } else {
                                  mod_po <- make_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                                      incidence_data = in_sample_mumps_case_reports)
                              }
                            } else if(covar_name == "slow") {
                              if (hypo_name == "gwaning") {
                                mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_slow, 
                                                          incidence_data = in_sample_mumps_case_reports)
                              } else if (hypo_name == "gwaning2") {
                                mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_slow, 
                                                              incidence_data = in_sample_mumps_case_reports)  
                              } else {
                                  mod_po <- make_pomp(covar = mod_mumps_covariates_slow, 
                                                      incidence_data = in_sample_mumps_case_reports)
                              }
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
                                            give_sys_soln = TRUE, 
                                            quantile_prob_names = c("0.1", "0.5", "0.90")) %.>% 
                              
                              mutate(., 
                                     hypothesis = hypo_name)
                            
                          })
)

hypo_factor_order <- c("No Loss", "Waning\n(Exponential)", "Waning\n(Erlang, x = 2)", 
                       "Waning\n(Erlang, x = 3)", "Leaky")

# add time-wise log-likelihood function to get a better idea of goodness of fit   

# prepare a combined covariates data to combine with the soln array
time_in_sample_data <- in_sample_mumps_case_reports$year  
  
covs_for_dmeasure <- (
  mod_mumps_covariates_constant %.>% 
    filter(., year %in% time_in_sample_data) %.>%
    select(., year, eta_a) %.>% 
    slice(., rep(1:n(), 5)) %.>% 
    mutate(., 
           hypothesis = rep(c("noloss", "leaky2", "waning", "gwaning2", "gwaning"), 
                            each = length(time_in_sample_data))
           ) %.>% 
    right_join(., 
               sim_from_these_tibble %.>% 
                 select(., hypothesis, starts_with("rho_age_"), starts_with("psi_")), 
               by = "hypothesis"
    ) %.>% 
    mutate(., 
           s_1 = rho_age_1*eta_a,
           s_2 = rho_age_2*eta_a,
           s_3 = rho_age_3*eta_a,
           s_4 = rho_age_4*eta_a,
           s_5 = rho_age_5*eta_a,
           s_u = rho_age_u*(1-eta_a),
           p_sq_1 = psi_1*psi_1,
           p_sq_2 = psi_2*psi_2,
           p_sq_3 = psi_3*psi_3,
           p_sq_4 = psi_4*psi_4,
           p_sq_5 = psi_5*psi_5,
           p_sq_u = psi_u*psi_u
           ) %.>% 
    select(., -c(eta_a, starts_with("rho_age_"), starts_with("psi_")))
    
  ) 
    
# treat the in-sample data to generate add as a covariate 
case_covs_for_dmeasure <- (
  in_sample_mumps_case_reports %.>% 
    select(., -total) %.>% 
    gather(., key = "age_class", value = "cases", -year) %.>% 
    mutate(., age_class = paste0("obs_", age_class)) %.>% 
    spread(., key = age_class, value = cases) %.>% 
    slice(., rep(1:n(), 5)) %.>% 
    mutate(., 
           hypothesis = rep(c("noloss", "leaky2", "waning", "gwaning", "gwaning2"), 
                               each = length(time_in_sample_data)))
  )


# generate time series of true cases under all models and add the covariate data
all_model_soln_data_w_covs <- (
  sim_all_models_pre %.>% 
    select(., year, age_class, hypothesis, sys_soln) %.>% 
    spread(., key = age_class, value = sys_soln) %.>% 
    right_join(., 
               covs_for_dmeasure,
               by = c("year", "hypothesis")) %.>% 
    right_join(., 
               case_covs_for_dmeasure,
               by = c("year", "hypothesis")) %.>% 
    select(., -unknown) %.>%
    group_by(., hypothesis) %.>% 
    mutate(., 
           total = `[0,5)`+`[5,15)`+`[15,25)`+`>40`, 
           m_1 = s_1*`[0,5)`, 
           m_2 = s_2*`[5,15)`, 
           m_3 = s_3*`[15,25)`, 
           m_4 = s_4*`[25,40)`, 
           m_5 = s_5*`>40`, 
           m_u = s_u*total, 
           v_1 = m_1*(1 - s_1 + p_sq_1*m_1), 
           v_2 = m_2*(1 - s_2 + p_sq_2*m_2), 
           v_3 = m_3*(1 - s_3 + p_sq_3*m_3), 
           v_4 = m_4*(1 - s_4 + p_sq_4*m_4), 
           v_5 = m_5*(1 - s_5 + p_sq_5*m_5), 
           v_u = m_u*(1 - s_u + p_sq_u*m_u), 
           `dmeas_[0,5)`   = dnorm(`obs_[0,5)` ,   mean = m_1, sd = (sqrt(v_1)+1e-18), log = TRUE), 
           `dmeas_[5,15)`  = dnorm(`obs_[5,15)`,   mean = m_2, sd = (sqrt(v_2)+1e-18), log = TRUE), 
           `dmeas_[15,25)` = dnorm(`obs_[15,25)`,  mean = m_3, sd = (sqrt(v_3)+1e-18), log = TRUE), 
           `dmeas_[25,40)` = dnorm(`obs_[25,40)`,  mean = m_4, sd = (sqrt(v_4)+1e-18), log = TRUE), 
           `dmeas_>40`     = dnorm(`obs_>40`   ,   mean = m_5, sd = (sqrt(v_5)+1e-18), log = TRUE), 
           `dmeas_unknown` = dnorm(`obs_unknown`,  mean = m_u, sd = (sqrt(v_u)+1e-18), log = TRUE), 
           `dmeas_total`   = `dmeas_[0,5)` + `dmeas_[5,15)` + `dmeas_[15,25)` + `dmeas_[25,40)` + 
             `dmeas_>40` + `dmeas_unknown`) %.>%
    select(., year, hypothesis, starts_with("dmeas")) #%.>% 
    #group_by(., hypothesis) %.>% 
    # transmute(., 
    #           year = year, 
    #           hypothesis = hypothesis, 
    #           `loglik_[0,5)`   = cumsum(coalesce(`dmeas_[0,5)`, 0)), 
    #           `loglik_[5,15)`  = cumsum(coalesce(`dmeas_[5,15)`, 0)),
    #           `loglik_[15,25)` = cumsum(coalesce(`dmeas_[15,25)`, 0)),
    #           `loglik_[25,40)` = cumsum(coalesce(`dmeas_[25,40)`, 0)),
    #           `loglik_>40`     = cumsum(coalesce(`dmeas_>40`, 0)),
    #           `loglik_unknown` = cumsum(coalesce(`dmeas_unknown`, 0)),
    #           `loglik`         = cumsum(coalesce(`dmeas_total`, 0))
    #           
    #           ) %.>% 
    # ungroup(.) %.>% 
    # mutate_all(., .funs = function(x) ifelse(x == 0, NA, x))
           
    
  ) 
  
# seprate age specific conditional logLiklihoods

age_mod_cond_loglik <- (
  all_model_soln_data_w_covs %.>%
  select(., -dmeas_total) %.>% 
  gather(., key = "age_class", value = "cond_loglik", -c(year, hypothesis)) %.>% 
  mutate(., 
         age_class = factor(age_class %.>% str_sub(., start = 7), 
                            levels = age_names_u), 
         hypothesis = case_when(hypothesis == "waning" ~ "Waning\n(Exponential)", 
                                hypothesis == "gwaning2" ~ "Waning\n(Erlang, x = 2)",
                                hypothesis == "gwaning" ~ "Waning\n(Erlang, x = 3)",
                                hypothesis == "leaky2" ~ "Leaky", 
                                hypothesis == "noloss" ~ "No Loss") %.>% 
           factor(., levels = hypo_factor_order)
         ) %.>% 
    group_by(., hypothesis, age_class) %.>% 
    ungroup(.)
  )
    

if(FALSE){
age_mod_cond_loglik %.>%
  ggplot(., aes(x = year, y = cond_loglik)) +
    geom_line() +
    facet_grid(rows = vars(age_class), cols = vars(hypothesis), scales = "free_y") +
    cap_axes()+
    project_theme
} 

  

sim_all_models <- (
  sim_all_models_pre %>%
    mutate(., hypothesis = case_when(hypothesis == "waning" ~ "Waning\n(Exponential)", 
                                     hypothesis == "gwaning2" ~ "Waning\n(Erlang, x = 2)",
                                     hypothesis == "gwaning" ~ "Waning\n(Erlang, x = 3)",
                                     hypothesis == "leaky2" ~ "Leaky", 
                                     hypothesis == "noloss" ~ "No Loss") %.>% 
             factor(., levels = hypo_factor_order)
           )
  )


# anno_data for delta AIC
anno_d_AIC <- (
  sim_from_these_tibble %.>% 
    select(., hypothesis, d_AIC) %.>% 
    mutate(., 
           hypothesis = case_when(hypothesis == "waning"~ "Waning\n(Exponential)", 
                                  hypothesis == "gwaning2" ~ "Waning\n(Erlang, x = 2)",
                                  hypothesis == "gwaning" ~ "Waning\n(Erlang, x = 3)",
                                  hypothesis == "leaky2" ~ "Leaky", 
                                  hypothesis == "noloss" ~ "No Loss") %.>% 
             factor(., levels = hypo_factor_order),
           age_class = "[0,5)")
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

# ratio for the secondary axis
sim_max <- sim_all_models$`0.90` %.>% max(., na.rm = TRUE)
cond_loglik_max <- -age_mod_cond_loglik$cond_loglik %.>% max(., na.rm = TRUE)

sec_axis_ratio <- cond_loglik_max/sim_max

# relative fits to the data 
rel_fit_plt <- (
  sim_all_models %.>% 
    mutate(., 
           hypothesis = as_factor(hypothesis)) %.>% 
    ggplot(.) +
    geom_line(data = in_sample_obs_data, aes(x = year, y = `0.5`, colour = "Observed Cases"), 
              size = 1) +
    #geom_line(aes(x = year, y = `0.5`, colour = hypothesis), size = 0.8) +
    geom_ribbon(aes(x = year, ymin = `0.1`, ymax = `0.90`, fill = hypothesis), alpha = 0.6) +
    geom_text(data = anno_d_AIC, aes(x = 1995, y = 175, label = paste0("Delta~AIC == ", 
                                                                       round(d_AIC, 2) %.>% 
                                                                         as.character(.))), 
              parse = TRUE) +
    geom_line(data = age_mod_cond_loglik %.>% na.omit(.), 
              aes(x = year, y = -cond_loglik/sec_axis_ratio, colour = hypothesis), size = 0.8) + 
    facet_grid(factor(age_class, levels = age_names_u)~hypothesis, scales = "fixed") +
    labs(x = "Year", 
         y = expression(sqrt(Cases)), 
         colour = "Case Trajectory/\nLog Density", 
         fill = "Prediction\nIntervals") +
    scale_y_continuous(sec.axis = sec_axis(~.*sec_axis_ratio, 
                                           name = expression(-log(P(D[t]~`|`~M[t])))
                                           )
                       )+
    scale_colour_manual(values = c("Observed Cases" = "grey30", 
                                   "No Loss" = "#654ea3", "Waning\n(Exponential)" = "#11998e", 
                                   "Waning\n(Erlang, x = 2)" = "blue", 
                                   "Waning\n(Erlang, x = 3)" = "#003973", 
                                   "Leaky" = "#FF416C"), 
                        breaks = c("Observed Cases", .$hypothesis %.>% levels(.))) +
    scale_fill_manual(values = c("No Loss" = "#654ea3", "Waning\n(Exponential)" = "#11998e", 
                                 "Waning\n(Erlang, x = 2)" = "blue", 
                                 "Waning\n(Erlang, x = 3)" = "#003973", 
                                 "Leaky" = "#FF416C")) +
    project_theme +
    cap_axes(right = "both") +
    guides(fill = guide_legend(nrow = 3, title.position = "top", order = 2), 
           colour = guide_legend(nrow = 3, title.position = "top", order = 1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  )



prop_obs <- (
  in_sample_obs_data %.>% 
    filter(., year > 2000) %.>% 
    group_by(., age_class) %.>%     
    summarize(., cases_obs = mean(`0.5`)) %.>%
    mutate(., prop_obs_cases = cases_obs/sum(cases_obs)) %.>% 
    ungroup(.) 
  )



# age distribution of cases for all models 
rel_case_distn_plt <- (
  sim_all_models %.>% 
  select(., year, age_class, `0.5`, hypothesis) %.>% 
  filter(., year > 2000) %.>% 
  group_by(., hypothesis, age_class) %.>% 
  summarize(., cases = mean(`0.5`)) %.>% 
  ungroup(.) %.>%
  group_by(., hypothesis) %.>% 
  mutate(., 
         prop_cases = cases/sum(cases), 
         age_class = factor(age_class, levels = age_names_u)) %.>% 
  ungroup(.) %.>% 
  right_join(., prop_obs, by = "age_class") %.>% 
  ggplot(.) +
  geom_bar(aes(x = age_class, y= prop_cases, fill = hypothesis), 
           stat = "identity", alpha = 0.5) +
  geom_point(aes(x = age_class, y = prop_obs_cases, colour = "Data")) +
  geom_line(aes(x = age_class, y = prop_obs_cases, group = hypothesis, colour = "Data")) +
  scale_fill_manual(values = c("No Loss" = "#654ea3", "Waning\n(Exponential)" = "#11998e", 
                               "Waning\n(Erlang, x = 2)" = "blue", 
                               "Waning\n(Erlang, x = 3)" = "#003973", 
                               "Leaky" = "#FF416C")
                    ) +
  scale_colour_manual(values = c("Data" = "grey30")) +
  scale_y_continuous(labels = function(x) scales::percent(x, suffix = "")) +
  facet_wrap(.~hypothesis, ncol = 2, scales = "fixed") +
  labs(x = "Age Class", 
       y = "Case Volume (%)",
       fill = "Hypothesis", 
       colour = "")+
  project_theme +
  cap_axes()+
  theme(legend.position = c(0.80, 0.1), 
        ) +
  guides(fill = guide_legend(title.position = "top", 
                             nrow = 3, 
                             order = 1), 
         colour = guide_legend(order = 2))  
  )
  
  
  

