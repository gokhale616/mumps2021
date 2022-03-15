# this script depends on ./prep_estm_tables.R
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


# generate a po_object
po_est_best_mod <- make_pomp(covar = mod_mumps_covariates_sigmoidal)


obs_sim <- sim_from_these %.>%   
  sim_obs_model(po_est_best_mod, params = ., times = time(po_est_best_mod), nsim = 1e3, 
                quantile_prob = c(0.1, 0.5, 0.90), 
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
                      name = "Prediction\nIntervals") +
    scale_colour_manual(values = c("grey30", model_col), 
                        labels = c("Observed", "Within\nSample", "Out of\nSample"), 
                        name = "Case\nTrajectory") +
    scale_x_continuous(breaks = gen_x_breaks, expand = c(0, 0)) +
    project_theme +
    cap_axes() +
    guides(colour = guide_legend(title.position = "top", ncol = 1, order = 1), 
           fill = guide_legend(title.position = "top", ncol = 1, order = 2, 
                               override.aes = list(alpha = 1))) +
    theme(text = element_text(size = unit(n_size, "pt")), 
          axis.text.x = element_text(angle = 90, vjust = 0.5))
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
    geom_text(data = Rsq_data, aes(x = 0.95, y = 3.65, label = paste0("R^2 == ",
                                                                   round(`in-sample`, 0) %.>%
                                                                     as.character(.), "*\`%\`")),
              size = 6,
              parse = TRUE,
              colour = model_col[1]) +
    geom_text(data = Rsq_data, aes(x = 0.95, y = 2.7, label = paste0("R^2 == ",
                                                                  round(`out-sample`, 0) %.>%
                                                                    as.character(.), "*\`%\`")),
              size = 6,
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
    cap_axes() +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal"),
           colour = guide_legend(title.position = "top", ncol = 1)) +
    theme(text = element_text(size = unit(n_size, "pt")), 
          axis.text.x = element_text(angle = 90, vjust = 0.5)) 
  )


fit_plot_grid <- plot_grid(fit_plot, Rsq_fit_plot, labels = c("A", "B"), align = "h", 
                           label_size = grid_lab_size)

# mixed effects model on the right panel
# gamma with infectious period and wanning duration
# 50% prediction intervals
# "wiretap" gradients on 'uigradients.com' has a good orange for out of sample prediction
