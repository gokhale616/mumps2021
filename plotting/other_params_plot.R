# this script depends on ./prep_estm_tables.R

aux_params <- (
  all_result_df %.>% 
    select(., c(d_AIC, hypothesis, starts_with("q_age"), 
                starts_with("rho_age"), starts_with("psi"))) %.>% 
    arrange(., d_AIC) %.>% 
    group_by(., hypothesis) %.>% 
    filter(., d_AIC == min(d_AIC)) %.>% 
    ungroup(.) %.>% 
    select(., -d_AIC) %.>% 
    gather(., key = "param", value ="value", -hypothesis) %.>% 
    mutate(., param = case_when(param == "psi_1"~"psi_age_1", 
                                param == "psi_2"~"psi_age_2", 
                                param == "psi_3"~"psi_age_3", 
                                param == "psi_4"~"psi_age_4", 
                                param == "psi_5"~"psi_age_5", 
                                param == "psi_u"~"psi_age_u", 
                                TRUE ~ param)) %.>% 
    mutate(., 
           aux_param = str_split_fixed(param, pattern = "_", 3)[,1], 
           age_cohort = str_split_fixed(param, pattern = "_", 3)[,3], 
           age_cohort = case_when(age_cohort == "1"~ age_names_u[1], 
                                  age_cohort == "2"~ age_names_u[2],
                                  age_cohort == "3"~ age_names_u[3],
                                  age_cohort == "4"~ age_names_u[4],
                                  age_cohort == "5"~ age_names_u[5],
                                  age_cohort == "u"~ age_names_u[6]
                                  ) %.>% factor(., levels = age_names_u), 
           aux_param =case_when(aux_param == "q" ~ "q", 
                                aux_param == "rho" ~ "rho", 
                                aux_param == "psi" ~ "psi") %.>% as_factor(.), 
           hypothesis =  case_when(hypothesis == "noloss" ~ "No Loss", 
                                   hypothesis == "waning" ~ "Waning\n(Exponential)", 
                                   hypothesis == "gwaning" ~ "Waning\n(Erlang, N = 3)", 
                                   hypothesis == "leaky2" ~ "Leaky") %.>% 
             factor(., levels = c("No Loss", "Waning\n(Exponential)", 
                                  "Waning\n(Erlang, N = 3)", "Leaky"))) %.>% 
    select(., -param)
  )


poi_plot <- (
  aux_params %.>% 
  filter(., aux_param == "q") %.>% 
  ggplot(., aes(x = age_cohort, 
                y  = value, 
                fill = age_cohort)) +
  labs(x = "", y = expression(atop("Probability Of", paste("Infection (", q, ")")))) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9)) +
  facet_grid(cols = vars(hypothesis),
             scales = "free_y") +
  scale_fill_manual(values = c(brewer_pal(palette = "Blues", direction = -1)(5)), 
                    guide= "none") +
  project_theme +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()) +
  coord_capped_cart(left = "both")
  )
  
  
rho_plot <- (
  aux_params %.>% 
    filter(., aux_param == "rho") %.>% 
    ggplot(., aes(x = age_cohort, 
                  y  = value, 
                  fill = age_cohort)) +
    labs(x = "", y = expression(atop("Reporting", paste("Probability (", rho, ")")))) +
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.9)) +
    facet_grid(cols = vars(hypothesis),
               scales = "free_y") +
    scale_fill_manual(values = c(brewer_pal(palette = "Blues", direction = -1)(5), 
                                 "darkturquoise"), 
                      guide= "none") +
    project_theme +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(), 
          strip.text.x = element_blank()) +
    coord_capped_cart(left = "both")
)


psi_plot <- (
  aux_params %.>% 
    filter(., aux_param == "psi") %.>% 
    ggplot(., aes(x = age_cohort, 
                  y  = value, 
                  fill = age_cohort)) +
    labs(x = "", 
         y = expression(atop("Dispersion", paste("Parameter (", psi, ")"))), 
         fill = "Age\nCohort") +
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.9)) +
    facet_grid(cols = vars(hypothesis),
               scales = "free_y") +
    scale_fill_manual(values = c(brewer_pal(palette = "Blues", direction = -1)(5), 
                                 "darkturquoise")) +
    project_theme +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(), 
          strip.text.x = element_blank()) +
    coord_capped_cart(left = "both")
)



aux_para_plot <- (
    poi_plot + rho_plot + psi_plot +
      plot_layout(
        guides = "collect",
        design = "
    A
    B
    C
    "
      ) + 
      plot_annotation(tag_levels = "A") &
      theme(legend.position='bottom') 
  )
