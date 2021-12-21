# this script is dependent on ./fit_plot.R
# it requires the object ./sim_from_these which contain mle of the best fitting model

# define a pomp obect to suumarize the epidemiology
po_to_sum <- make_pomp(., 
                       covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1965-1/52, extra_end_t = 2065, 
                       temp_scale = 1/52)

# Generate trajectories from all the compartments 
comp_traj <- trajectory(po_to_sum, 
                        param = sim_from_these, 
                        method = "ode23", format = "d")



# generate an interpolated dataframe with the population sizes

interp_range <- c(1965, 2065)

xnew <- seq(interp_range[1], interp_range[2], by = 1/52)

interpolated_covs <- (
  mod_mumps_covariates_sigmoidal %.>% 
    bind_rows(., 
              mod_mumps_covariates_sigmoidal %.>% 
                slice(., rep(n(), interp_range[2]-2018))) %.>% 
    mutate(., 
           year = seq(1910, interp_range[2], by = 1)) %.>% 
    filter(., year > 1949) %.>%   
    select(., year, starts_with("N_")) %.>% 
    interp.dataset(y = ., x=.$year, 
                   xout = xnew, 
                   method = "linear") %.>% 
    as_tibble(.)
  )  


# Join the two data frames together for compartmental goodness!

comp_traj_w_cov <- (
  comp_traj %.>% 
    select(., -`.id`) %.>% 
    filter(., year > 1950-1/52) %.>% 
    right_join(., interpolated_covs, by ="year") %.>% 
    drop_na(.) %.>% 
    mutate(., 
           R_1 = N_1 - (S_1 + V_1 + E1_1 + E2_1 + I1_1 + I2_1), 
           R_2 = N_2 - (S_2 + V_2 + E1_2 + E2_2 + I1_2 + I2_2), 
           R_3 = N_3 - (S_3 + V_3 + E1_3 + E2_3 + I1_3 + I2_3), 
           R_4 = N_4 - (S_4 + V_4 + E1_4 + E2_4 + I1_4 + I2_4), 
           R_5 = N_5 - (S_5 + V_5 + E1_5 + E2_5 + I1_5 + I2_5))
  )


# only compartments that are not infection related
comp_traj_F <- (
  comp_traj_w_cov %.>% 
  mutate(., 
         S_1 = S_1/N_1, S_2 = S_2/N_2, S_3 = S_3/N_3, S_4 = S_4/N_4, S_5 = S_5/N_5, 
         V_1 = V_1/N_1, V_2 = V_2/N_2, V_3 = V_3/N_3, V_4 = V_4/N_4, V_5 = V_5/N_5, 
         R_1 = R_1/N_1, R_2 = R_2/N_2, R_3 = R_3/N_3, R_4 = R_4/N_4, R_5 = R_5/N_5, 
         `.id` = 1) %.>% 
  select(., year, `.id`, starts_with("S_"), starts_with("V_"), starts_with("R_")) %.>% 
  prep_as_plot_trajectory(., init_year = 1950-1/52) %.>%  
  mutate(., 
         comp_exp = case_when(comp_exp == "S" ~ "Susceptible", 
                              comp_exp == "V" ~ "Vaccinated", 
                              comp_exp == "R" ~ "Recovered", 
                              TRUE ~ comp_exp) %.>% factor(., 
                                                           level = c("Susceptible", "Recovered", 
                                                                      "Vaccinated")))  
    
)


year_break_x <- seq(1965, 2065, 20)

comp_traj_F_plot <- (
  comp_traj_F %.>%
  mutate(., count = case_when(count > 1 ~ 1,
                              count < 0 ~ 0,
                              TRUE ~ count)
         ) %.>%
  ggplot(., aes(x = year, y = count, colour = age_cohort)) +
  geom_line(size = 1.2) +
  facet_wrap(vars(comp_exp), scales = "fixed") +
  project_theme +
  labs(y = "Percent of Total", 
       x = "", 
       colour = "Age\nCohort") +
  scale_color_brewer(palette = "Oranges", direction = -1) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = year_break_x) +  
  cap_axes +
  theme(legend.position = c(0.2, 0.7)) +
  guides(colour = guide_legend(nrow = 3))
  )



# only compartment that are infectious 
comp_traj_V <- (
  comp_traj_w_cov %.>% 
    mutate(.,
           Ns_1 = 1e5/N_1, Ns_2 = 1e5/N_2, Ns_3 = 1e5/N_3, Ns_4 = 1e5/N_4, Ns_5 = 1e5/N_5, 
           E_1 = (E1_1 + E2_1)*Ns_1, E_2 = (E1_2 + E2_2)*Ns_2, E_3 = (E1_3 + E2_3)*Ns_3, 
           E_4 = (E1_4 + E2_4)*Ns_4, E_5 = (E1_5 + E2_5)*Ns_5, 
           I_1 = (I1_1 + I2_1)*Ns_1, I_2 = (I1_2 + I2_2)*Ns_2, I_3 = (I1_3 + I2_3)*Ns_3, 
           I_4 = (I1_4 + I2_4)*Ns_4, I_5 = (I1_5 + I2_5)*Ns_5, 
           C_1 = C_1*Ns_1, C_2 = C_2*Ns_2, C_3 = C_3*Ns_3, C_4 = C_4*Ns_4, C_5 = C_5*Ns_5, 
           `.id` = 1) %.>% 
    select(., year, `.id`, starts_with("E_"), starts_with("I_"), starts_with("C_")) %.>% 
    prep_as_plot_trajectory(., init_year = 1950-1/52) %.>%  
    mutate(., 
           comp_exp = case_when(comp_exp == "E" ~ "Exposed", 
                                comp_exp == "I" ~ "Infectious", 
                                comp_exp == "C" ~ "True Cases", 
                                TRUE ~ comp_exp) %.>% factor(., 
                                                             level = c("Exposed", "Infectious", 
                                                                       "True Cases")))  
  
    
    )
  
         
  
comp_traj_V_plot <- (
  comp_traj_V %.>%
    ggplot(., aes(x = year, y = count, colour = age_cohort)) +
    geom_line(size = 1.2) +
    facet_wrap(vars(comp_exp), scales = "fixed") +
    project_theme +
    labs(y = expression(log[10](Count~Per~100000)), 
         x = "Year", 
         colour = "Age\nCohort") +
    scale_color_brewer(palette = "Oranges", direction = -1) +
    scale_y_continuous(trans =  "log10", breaks = c(1e-2, 1e-1, 1e0, 1e1, 1, 1e2, 1e3)) +
    scale_x_continuous(breaks = year_break_x) +
    cap_axes +
    theme(legend.position = "none") +
    guides(colour = guide_legend(nrow = 3))
)


plot_grid(comp_traj_F_plot, comp_traj_V_plot, nrow = 2, align = "v", labels = c("A", "B"))






