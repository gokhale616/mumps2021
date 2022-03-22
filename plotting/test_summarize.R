best_model_p_vec <- (
  all_result_df %.>% 
    filter(., best_fit_covar == 1) %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar, 
                 epsilon2, p_intro, t_intro, impact)) %.>% 
    unlist(.) %.>% sim_p_vals(., default_p_vals = param_vals_vacc_effective)
)


# define a pomp obect to summarize the epidemiology
po_to_sum <- (
  make_pomp_vacc_effective(., 
                           covar = mod_mumps_covariates_sigmoidal, 
                           temp_scale = 1/1, 
                           extra_start_t = 1967-1/52,
                           extra_end_t = 2020)
)

# Generate trajectories from all the compartments 
comp_traj <- (
  trajectory(po_to_sum, 
             param = best_model_p_vec, 
             method = "ode23", format = "d") 
) %.>% 
  slice(., -1)


V_prop_traj <- (
  comp_traj %>% 
    mutate(., 
           Vs = Vs_1 + Vs_2 + Vs_3 + Vs_4 + Vs_5, 
           V = V_1 + V_2 + V_3 + V_4 + V_5,
           N = Nv_1 + Nv_2 + Nv_3 + Nv_4 + Nv_5, 
           Vprop = (V)/N, 
           Vprop_1 = (V_1)/Nv_1,
           Vprop_2 = (V_2)/Nv_2, Vprop_3 = (V_3)/Nv_3,
           Vprop_4 = (V_4)/Nv_4, Vprop_5 = (V_5)/Nv_5,
           Vsprop_1 = (Vs_1)*1e5/Nv_1,
           Vsprop_2 = (Vs_2)*1e5/Nv_2, Vsprop_3 = (Vs_3)*1e5/Nv_3,
           Vsprop_4 = (Vs_4)*1e5/Nv_4, Vsprop_5 = (Vs_5)*1e5/Nv_5,
    ) %.>% 
    select(., `.id`, 
           year, starts_with("Vprop"), starts_with("Vsprop"),
           starts_with("Nv"), starts_with("V")
           )
  
)


V_prop_traj_for_plot <- (
  V_prop_traj %.>% 
    select(., -Vprop) %.>% 
    prep_as_plot_trajectory(., init_year = 1967-1/52) 
  
)


anno_Vs_tot_prop_data <- (
  V_prop_traj %.>% 
    select(., year, Vprop) #%.>% 
  #mutate(., year = floor(year)) %.>% 
  #group_by(., year) %.>% 
  #summarize_all(., sum)
)


year_break_x <- seq(1970, fin_year, 15)

Vs_plot <- (
  V_prop_traj_for_plot %.>%
    filter(., comp_exp == "Vsprop") %.>% 
    ggplot(.) +
    geom_line(aes(x = year, y = count, 
                  colour = age_cohort), size = 0.8) +
    labs(y = expression(Immunity~Lost~Per~10^5), x = "Year", 
         color = "Age\nCohort") +
    scale_y_continuous() +
    scale_x_continuous(breaks = year_break_x) +
    scale_colour_brewer(palette = "Oranges", direction = -1, guide = "none") +    
    project_theme +
    cap_axes() +
    theme() 
)

orange_age_cohort <- brewer_pal(palette = "Oranges", direction = -1)(5)

V_prop_plot <- (
  V_prop_traj_for_plot %.>%
    filter(., comp_exp == "Vprop") %.>% 
    ggplot(.) +
    geom_line(aes(x = year, y = count, 
                  colour = age_cohort), size = 0.8) +
    geom_line(data = anno_Vs_tot_prop_data, 
              aes(x = year, y = Vprop, colour = "total"), size = 0.8) +
    annotate(geom = "segment", 
             x = 1967, xend = 2020, 
             y = critic_lev_vacc, yend = critic_lev_vacc,
             colour = "#6be585", linetype = "dashed", size = 0.9) +
    annotate(geom = "text", 
             label = "Critical Vaccination Level",
             y = critic_lev_vacc+0.065, x = 1992, 
             colour = "grey30", size = 2.5) +
    labs(x = "", y = "Vaccinated (%)") +
    #scale_y_continuous(labels = function(x) scales::percent(x, suffix = ""), 
                       #limits = c(0, 1), 
                       #breaks = seq(0, 1,by = 0.2)) +
    scale_x_continuous(breaks = year_break_x) +
    scale_colour_manual(values = c(orange_age_cohort, "#4b6cb7"), name = "Age\nCohort", 
                        breaks = c(age_names, "total")) +
    project_theme +
    cap_axes() +
    guides(colour = guide_legend(nrow = 2, title.position = "left"))
) 



# incidence
comp_traj %.>% 
  mutate(., 
         Incs_1 = (Is_1 + Iw_1)*1e5/Nv_1, 
         Incs_2 = (Is_2 + Iw_2)*1e5/Nv_2, 
         Incs_3 = (Is_3 + Iw_3)*1e5/Nv_3, 
         Incs_4 = (Is_4 + Iw_4)*1e5/Nv_4, 
         Incs_5 = (Is_5 + Iw_5)*1e5/Nv_5) %.>% 
  select(., `.id`, year, starts_with("Incs")) %.>% 
  prep_as_plot_trajectory(., init_year = 1967-1/52) %.>% 
  ggplot(., aes(x = year)) +
  geom_line(aes(y = count, colour = age_cohort), size = 0.8) +
  labs(y = expression(Prevalence~per~10^5), 
       x = "Year") +
  scale_x_continuous(breaks = year_break_x) +
  scale_y_continuous(trans = "log10", breaks = c(1e-3, 1e-2, 1e-1, 1e0, 
                                                 1e1, 1e2, 1e3), 
                     limits = c(1e-3, 1e3)) +
  annotation_logticks(sides = "l")  +
  scale_colour_brewer(palette = "Oranges", direction = -1, guide = "none") +
  project_theme +
  cap_axes()


comp_traj %.>% 
  mutate(., 
         Sp_1 = (Ss_1 + Sw_1)/Nv_1, 
         Sp_2 = (Ss_2 + Sw_2)/Nv_2, 
         Sp_3 = (Ss_3 + Sw_3)/Nv_3, 
         Sp_4 = (Ss_4 + Sw_4)/Nv_4, 
         Sp_5 = (Ss_5 + Sw_5)/Nv_5) %.>% 
  select(., `.id`, year, starts_with("Sp_")) %.>% 
  prep_as_plot_trajectory(., init_year = 1967-1/52) %.>% 
  mutate(., count = ifelse(count > 1, 1, count)) %.>% 
  ggplot(., aes(x = year)) +
  geom_line(aes(y = count, colour = age_cohort), size = 0.8) +
  labs(y = "Susceptible (%)", 
       x = "Year") +
  scale_x_continuous(breaks = year_break_x) +
  scale_y_continuous(labels = function(x) scales::percent(x, suffix = "")) +
  scale_colour_brewer(palette = "Oranges", direction = -1, guide = "none") +
  project_theme +
  cap_axes()




