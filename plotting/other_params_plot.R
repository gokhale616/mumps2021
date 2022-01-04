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
           aux_param =case_when(aux_param == "q" ~ "Probability~Of~Infection~(q)", 
                                aux_param == "rho" ~ "Reporting~Probabilty~(rho)", 
                                aux_param == "psi" ~ "Dispersion~Parameter~(psi)") %.>% as_factor(.), 
           hypothesis =  case_when(hypothesis == "noloss" ~ "No Loss", 
                                   hypothesis == "waning" ~ "Waning\n(Exponential)", 
                                   hypothesis == "gwaning" ~ "Waning\n(Erlang, N = 3)", 
                                   hypothesis == "leaky2" ~ "Leaky") %.>% 
             factor(., levels = c("No Loss", "Waning\n(Exponential)", 
                                  "Waning\n(Erlang, N = 3)", "Leaky"))) %.>% 
    select(., -param)
  )


aux_para_plot <- (
  aux_params %.>% 
    ggplot(., aes(x= hypothesis, y  = value, fill = age_cohort)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(rows = vars(aux_param), scales = "fixed", 
               labeller = label_parsed) +
    geom_text(aes(y = 1e1, label = round(value, 3)), 
               position = position_dodge(width = 0.9), 
               angle = 90) +
    labs(., 
         x = "Model", y = "Estimate Value", fill = "Age\nCohort") +
    scale_y_continuous(trans = "log10", limits = c(1e-4, 5e1)) +
    scale_fill_manual(values = c(brewer_pal(palette = "Blues", direction = -1)(5), "darkturquoise")) +
    project_theme +
    cap_axes())


  


