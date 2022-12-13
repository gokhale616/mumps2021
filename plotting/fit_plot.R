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
                quantile_prob_names = c("0.1", "0.5", "0.9")) %.>% 
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
    geom_ribbon(aes(ymin = `0.1`, ymax = `0.9`, fill = sample), alpha= 0.6) +  
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
          axis.text.x = element_text(angle = 90, vjust = 0.5), 
          panel.spacing = unit(0.8, "lines"))
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
    geom_text(data = Rsq_data, aes(x = 0.95, y = 3.75, label = paste0("R^2 == ",
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
    scale_x_continuous(breaks  = c(0, 1, 2, 3, 4), limits = c(0,4.5)) +
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
          axis.text.x = element_text(angle = 90, vjust = 0.5), 
          panel.spacing = unit(0.8, "lines")
          ) 
  )


fit_plot_grid <- plot_grid(fit_plot, Rsq_fit_plot, labels = c("A", "B"), align = "h", 
                           label_size = grid_lab_size)

# supplementary plot to test why out why out of fit samples are better 
calculate_Rsq_temp <- function(data) {
  cor(data$`0.5`,data$`0.5_sim`)^2*100
}

rollingRsq <- (
  age_specific_Rsq %.>% 
    group_by(., age_class) %.>%   
    mutate(., 
           roll_or = slider::slide_dbl(.x = cur_data(), .f = ~calculate_Rsq_temp(.x), 
                                       .before = 5,
                                       .complete = TRUE)
           )
  )
    
    
Rsq_data %.>% 
  gather(., key = "sample", value = "Rsq", -age_class)
  


Rsq_roll_plot <- (
  rollingRsq %.>%   
    mutate(., 
           Rsq = Rsq/100,
           roll_or = roll_or/100,
           sample = ifelse(sample == "in-sample", "Within\nsample", "Out of\nsample")) %.>% 
    ggplot(.)+
    geom_line(aes(x = year, y = roll_or, colour = "6 Year\nRolling\nWindow"), size = 0.8)+
    geom_line(aes(x = year, y = Rsq, colour = sample), linetype = "dashed", size = 0.8) +
    facet_grid(rows= vars(age_class), scales = "free_x")+
    labs(x = "Year", 
         y = expression(atop(Coefficent~of, Determination~(R^2))))+
    scale_colour_manual(name = "", 
                        values = c("6 Year\nRolling\nWindow" = "#1f77b4", 
                                   "Within\nsample" = model_col[1], 
                                   "Out of\nsample" = model_col[2]))+
    scale_y_continuous(breaks = c(0,.50,1), labels = scales::percent)+
    scale_x_continuous(breaks = c(1977,1984,1991,1998,2005,2012,2018), limits = c(1977, 2018))+
    project_theme+
    cap_axes()
  )


# another plot to evaluates the generate a boot strapped distribution of R^2

data_distn_to_compare <- sim_from_these %.>%   
  sim_obs_model(po_est_best_mod, params = ., times = time(po_est_best_mod), nsim = 1e3, 
                summarise_obs_process = FALSE) %.>% 
  transmute(., 
            `.id` = `.id`, 
            year = year, 
            age_class = age_class, 
            sim = cases+1,
            hypothesis = "Waning", 
            sample = ifelse(year<2013, "in-sample", "out-sample")) %.>% 
  left_join(., 
            by = c("year","age_class"), 
            obs_data %.>% 
              transmute(., 
                        year = year, 
                        age_class = age_class, 
                        obs = (`0.5`+1)^2
                        )
            )


point_error_bar_data <- (
  data_distn_to_compare %.>% 
    group_by(., year, age_class) %.>% 
    summarise(., 
              qs = quantile(sim, c(0.1, 0.5, 0.9), na.rm = TRUE), 
              prob = c("0.01", "0.5", "0.9"), 
              .groups = 'drop') %.>%
    spread(., key = prob, value = qs) %.>% 
    left_join(., 
              by = c("year","age_class"), 
              obs_data %.>% 
                transmute(., 
                          year = year, 
                          age_class = age_class, 
                          obs = `0.5`+1
                )
    ) %.>% 
    gather(., key = "stat", value = "cases", -c(year, age_class)) %.>% 
    mutate(., 
           cases = log(cases, base = 10)) %.>% 
    spread(., key = stat, value = cases) %.>% 
    mutate(., 
           hypothesis = "Waning", 
           sample = ifelse(year<2013, "in-sample", "out-sample")
           )
)


if(FALSE){
# calculate and plot bootstrapped distribution of R^2
Rsq_distn_data <- (
  data_distn_to_compare %.>% 
    drop_na(.) %.>% 
    mutate(., 
           log_sim = log(sim, base = 10), 
           log_obs = log(obs, base = 10)) %.>% 
    select(., -c(sim, obs)) %.>% 
    group_by(., `.id`, age_class, sample) %.>% 
    mutate(., 
           Rsq = cor(log_sim, log_obs)^2*100
           )
  )


Rsq_distn_summary_data <- (
  Rsq_distn_data %.>% 
    group_by(., age_class, sample) %.>%
    summarise(.,
             qs = quantile(Rsq, c(0.025, 0.5, 0.975), na.rm = TRUE),
             prob = c("Rsq_0.025", "Rsq_0.5", "Rsq_0.975"),
             .groups = 'drop',
             ) %.>%
    ungroup(.) %.>%
    mutate(., 
           sample = ifelse(sample == "in-sample", "is", "os")) %.>% 
    pivot_wider(., names_from = c(sample, prob), values_from = qs) %.>% 
    mutate(., 
           `Age Cohort` = age_class,
           `In-Sample` = paste0(round(is_Rsq_0.5, 1), "\\% (", round(is_Rsq_0.025, 1), "\\%, ", round(is_Rsq_0.975,1), "\\%)"), 
           `Out-Of-Sample` = paste0(round(os_Rsq_0.5, 1), "\\% (", round(os_Rsq_0.025, 1), "\\%, ", round(os_Rsq_0.975,1), "\\%)")) %.>% 
    select(., `Age Cohort`, `In-Sample`, `Out-Of-Sample`)
)  
  




Rsq_kbl <- (
  Rsq_distn_summary_data %.>%   
    kbl(., 
        align = "c", 
        digits = 4, 
        linesep = "",
        booktabs = T, 
        format = "latex", 
        caption = "Age Specific Coeficent of Determinition", 
        escape = FALSE) %.>% 
    kable_styling(.,
                  position='left', full_width = F,
                  latex_options=c('striped', 'HOLD_position', "scale_down"))
)


Rsq_fit_plot_2 <- (
  point_error_bar_data %.>% 
    ggplot(.) +
    geom_errorbar(aes(x = obs, ymin = `0.01`, ymax = `0.9`, colour = sample), size = 0.2)+
    geom_point(aes(x = obs, y = `0.5`, colour = sample, fill = year), size = 1.5, pch = 21) +
    geom_smooth(aes(x = obs, y = `0.5`, colour = sample, fill = year), method = "lm", se = FALSE) +
    facet_grid(rows = vars(age_class), 
               scales = "fixed") +
    labs(x = expression(log[10](Observed~Cases+1)), 
         y = expression(log[10](Simulated~Cases+1)), 
         fill = "Year") +
    scale_x_continuous(breaks  = c(0, 1, 2, 3, 4), limits = c(0,4.5)) +
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
          axis.text.x = element_text(angle = 90, vjust = 0.5), 
          panel.spacing = unit(0.8, "lines")
    ) 
)

Rsq_distn_plot <- (
  Rsq_distn_data %.>% 
    ggplot(., aes(x = sample, y = Rsq, colour = sample))+
    geom_violin() +
    geom_boxplot(width=0.05)+
    facet_grid(rows = vars(age_class), 
               scales = "fixed") +
    labs(y = expression(paste("Coefficent Of Variation (", R^2, ", %)")))+
    scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100))+
    scale_colour_manual(values = model_col, 
                        labels = c("Within\nsample", "Out of\nsample"), 
                        name = "Log-linear\nModel Fit") +
    project_theme +
    cap_axes() +
    guides(colour = guide_legend(title.position = "top", ncol = 1)) +
    theme(axis.title.x =element_blank(), 
          axis.text.x =element_blank(), 
          #axis.ticks.x = element_blank()
          )
  )
  
  


fit_plot_grid_2 <- (  
  fit_plot +
  Rsq_fit_plot_2 + 
    Rsq_distn_plot +      
    plot_layout(
      width = c(1, 1, 0.5),
      #heights = c(1, 1, 0.8),
      guides = "collect", 
      design = "
      ABC
      " 
    ) +
    plot_annotation(tag_levels = "A")&
    theme(plot.tag = element_text(size = 10, face = "bold"), 
          text = element_text(size = unit(9, "pt")), 
          legend.position = "bottom")
  )
    
}


