# load all of the result objects in the global environment 
mle_result_path <- "../result_data/mle"

list.files(path = mle_result_path, 
           full.names = TRUE)[-21] %.>% 
  lapply(., load,  envir = .GlobalEnv)



# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_waning_slow_is, mle_waning_sigmoid_is, mle_waning_rapid_is, mle_waning_constant_is,
    mle_leaky2_2_slow_is, mle_leaky2_2_sigmoid_is, mle_leaky2_2_rapid_is, mle_leaky2_2_constant_is,
    mle_leaky2_3_slow_is, mle_leaky2_3_sigmoid_is, mle_leaky2_3_rapid_is, mle_leaky2_3_constant_is,
    mle_leaky2_4_slow_is, mle_leaky2_4_sigmoid_is, mle_leaky2_4_rapid_is, mle_leaky2_4_constant_is, 
    mle_slow_is, mle_sigmoid_is, mle_rapid_is, mle_constant_is
  )
)


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
    mutate(., 
           best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC)
)


# table of estimates - process model 
table_hypo_compare <- (
  all_result_df %.>% 
    select(., -c(best_fit_covar, starts_with("q_age"), 
                 starts_with("rho_age"), starts_with("psi"), 
                 p_intro)) %.>% 
    arrange(., d_AIC) %.>% 
    group_by(., hypothesis) %.>% 
    filter(., d_AIC == min(d_AIC)) %.>% 
    ungroup(.) %.>%
    mutate_if(., is.numeric, function(x){round(x, digits = 3)}) %.>%   
    mutate(., 
           vacc_covariate = str_to_title(vacc_covariate),
           hypothesis = str_to_title(hypothesis)) %.>% 
    gather(., 
           key = "Parameter", value = "Estimate", 
           -c(hypothesis)) %.>% 
    spread(., key = hypothesis, value = Estimate) 
)

# table_hypo_compare %.>% xtable(.)
  

sim_from_these <- (
  all_result_df %.>%   
    filter(., best_fit_covar == TRUE) %.>% 
    select(., 
           -c(hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
    mutate(., 
           p_intro = 6, 
           t_intro = 3000, 
           epsilon2 = 0) %.>% 
    unlist(.) %.>% 
    sim_p_vals(.)
    
)

# make this happen nicely  
source("../fit/treat_vacc_covar.R", chdir = TRUE)

# generate a po_object

po_est_test <- make_pomp(covar = mod_mumps_covariates_sigmoidal)


obs_sim <- sim_from_these %.>%   
  sim_obs_model(po_est_test, params = ., times = time(po_est_test), nsim = 1e3, 
                quantile_prob = c(0.10, 0.5, 0.90), 
                quantile_prob_names = c("0.025", "0.5", "0.975")) %.>% 
  mutate(., 
         hypothesis = "Waning", 
         sample = ifelse(year<2013, "in-sample", "our-sample"))

  
obs_data <- mumps_case_reports %.>%  
  select(., -total) %.>% 
  gather(., key = "age_class", value = `0.5`, -year, factor_key = TRUE) %.>% 
  mutate(., 
         sample = ifelse(year<2013, "in-sample", "our-sample"),
         `0.5` = sqrt(`0.5`)) %.>% 
  drop_na(.) 


model_col <- c("#11998e", "#fe8c00")

fit_plot <- (
  obs_sim %.>%   
    ggplot(., aes(x = year)) +
    geom_line(data = obs_data, 
              aes(y = `0.5`, colour = "data", group = sample), size = 0.8) +  
    geom_line(aes(y = `0.5`, colour = sample), size = 1.0) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`, fill = sample), alpha= 0.6) +  
    labs(y = expression(sqrt(Cases)), 
         x = "Year", 
         fill = "", 
         colour = "") +
    facet_grid(rows = vars(age_class), 
               scales = "fixed") +
    scale_fill_manual(values = model_col, 
                      labels = c("Within\nsample", "Out of\nsample"), 
                      name = "90% Prediction\nIntervals") +
    scale_colour_manual(values = c("grey30", model_col), 
                        labels = c("Observed", "Simulated-in", "Simulated-out"), 
                        name = "Cases") +
    scale_x_continuous(breaks = gen_x_breaks, expand = c(0, 0)) +
    project_theme +
    cap_axes +
    guides(colour = guide_legend(title.position = "top", ncol = 1), 
           fill = guide_legend(title.position = "top", ncol = 1))  
  )
  

# Rsq for the models 
data_to_compare <- (
  obs_sim %.>% 
  select(., year, age_class, `0.5`) %.>% 
  transmute(., 
            year = year, age_class = age_class, `0.5` = log(`0.5`^2, base = 10)) %.>% 
  left_join(., 
            by = c("year","age_class"), 
            obs_data %.>% 
              transmute(., 
                        year = year, 
                        age_class = age_class, 
                        `0.5_sim` = log(`0.5`^2, base = 10)))
  ) 
  
age_specific_Rsq <- (
  data_to_compare %.>%
    drop_na(.) %.>% 
    mutate(., 
           sample = ifelse(year<2013, "in-sample", "out-sample")) %.>% 
    group_by(., age_class, sample) %.>% 
    mutate(., 
           Rsq = cor(`0.5`, `0.5_sim`)^2*100, 
           )  
    ) 

Rsq_data <- (
  age_specific_Rsq %.>% 
    filter(., year %in% c(1977, 2014)) %.>% 
    select(., age_class, Rsq, sample) %.>% 
    spread(., key = sample, value = Rsq) %.>% 
    ungroup(.)
  )


Rsq_fit_plot <- (
  age_specific_Rsq %.>% 
    ggplot(.) +
    geom_point(aes(x = `0.5`, y = `0.5_sim`, colour = sample, fill = year), size = 1.5, pch = 21) +
    geom_smooth(aes(x = `0.5`, y = `0.5_sim`, colour = sample, fill = year), method = "lm", se = FALSE) +
    facet_grid(rows = vars(age_class), 
               scales = "fixed") +
    geom_text(data = Rsq_data, aes(x = 0.5, y = 3.5, label = paste0("R^2 == ", 
                                                                   round(`in-sample`, 0) %.>% 
                                                                     as.character(.), "*\`%\`")), 
              parse = TRUE, 
              colour = model_col[1]) +
    geom_text(data = Rsq_data, aes(x = 0.5, y = 2.3, label = paste0("R^2 == ", 
                                                                  round(`out-sample`, 0) %.>% 
                                                                    as.character(.), "*\`%\`")), 
              parse = TRUE, 
              colour = model_col[2]) +
    labs(x = expression(log[10](Observed~Cases)), 
         y = expression(log[10](Median~Simulated~Cases)), 
         fill = "Year") +
    scale_x_continuous(breaks  = c(0, 1, 2, 3, 4), limits = c(0,4)) +
    scale_fill_gradient(low = "white", high = "grey30", 
                        limits = c(1977, 2018), breaks = c(1977, 1997, 2018)) +
    scale_colour_manual(values = model_col, 
                        labels = c("Within\nsample", "Out of\nsample"), 
                        name = "Log-linear\nModel Fit") +
    project_theme + 
    cap_axes +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal"), 
           colour = guide_legend(title.position = "top", ncol = 1))
  )


fit_plot_grid <- plot_grid(fit_plot, Rsq_fit_plot, labels = c("A", "B"), align = "h")

# mixed effects model on the right panel
# gamma with infectious period and wanning duration
# 50% prediction intervals
# "wiretap" gradients on 'uigradients.com' has a good orange for out of sample prediction
