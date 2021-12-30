# this script is dependent on ./fit_plot.R
# it requires the object ./sim_from_these which contain mle of the best fitting model

fin_year <- 2020

# define a pomp obect to suumarize the epidemiology
po_to_sum <- make_pomp(., 
                       covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1965-1/52, extra_end_t = fin_year, 
                       temp_scale = 1/52)

# Generate trajectories from all the compartments 
comp_traj <- trajectory(po_to_sum, 
                        param = sim_from_these, 
                        method = "ode23", format = "d")



# generate an interpolated dataframe with the population sizes

interp_range <- c(1965, fin_year)

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
    mutate(., year = floor(year)) %.>%  
    group_by(., year) %.>%   
    summarize_all(., mean) %.>%   
    ungroup(.) %.>% 
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

# only compartment that are infectious 
comp_traj_V <- (
  comp_traj_w_cov %.>% 
    mutate(., year = floor(year)) %.>%  
    group_by(., year) %.>%   
    summarize_all(., mean) %.>%   
    ungroup(.) %.>% 
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

# ### .. calculate the mean age at infection
# mean_age_inf_data <- (
#   comp_traj_w_cov %.>% 
#     filter(., year > 1965 & year < 2020) %.>% 
#     mutate(., `.id` = 1) %.>% 
#     select(., year, `.id`, starts_with("I1_")) %.>% 
#     prep_as_plot_trajectory(., init_year = 1965-1/52) %.>% 
#     group_by(., year) %.>% 
#     mutate(., wt = count/sum(count)) %.>% 
#     ungroup(.) %.>% 
#     mutate(., mean_cohort_age = case_when(age_cohort == age_names[1] ~ (4+0)/2, 
#                                           age_cohort == age_names[2] ~ (14+5)/2, 
#                                           age_cohort == age_names[3] ~ (24+15)/2, 
#                                           age_cohort == age_names[4] ~ (39+25)/2, 
#                                           age_cohort == age_names[5] ~ (40+80)/2)) %.>% 
#     select(., year, mean_cohort_age, wt) %.>% 
#     group_by(., year) %.>% 
#     mutate(., mean_age_at_inf  = sum(mean_cohort_age*wt)) %.>% 
#     select(., year, mean_age_at_inf) %.>% 
#     ungroup()
# )

year_break_x <- seq(1970, fin_year, 15)

plot_comp_summary <- function(traj_data, y_label, legend_position = c(0.2, 0.7), ...) {
  
  traj_data %.>% 
    ggplot(., aes(x = year, y = count, colour = age_cohort)) +
    geom_line(size = 1) +
    project_theme +
    labs(y = y_label, 
         x = "", 
         colour = "Age\nCohort") +
    scale_color_brewer(palette = "PuBuGn", direction = -1) +
    scale_y_continuous(...) +
    scale_x_continuous(breaks = year_break_x) +  
    cap_axes +
    theme(legend.position = legend_position) +
    guides(colour = guide_legend(nrow = 3))
    
  
}


pSusc_plot <-(
  comp_traj_F %.>% 
    mutate(., count = case_when(count > 1 ~ 1,
                                count < 0 ~ 0,
                                TRUE ~ count)
    ) %.>% 
    filter(., comp_exp == "Susceptible") %.>%  
    plot_comp_summary(., 
                      y_label = "Percent Susceptible", 
                      legend_position = c(0.5, 0.7), 
                      labels = scales::percent, limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1))
  )

pRecv_plot <-(
  comp_traj_F %.>% 
    mutate(., count = case_when(count > 1 ~ 1,
                                count < 0 ~ 0,
                                TRUE ~ count)
    ) %.>%
    filter(., comp_exp == "Recovered") %.>%  
    plot_comp_summary(., y_label = "Percent Recovered", legend_position = "none", 
                      labels = scales::percent, limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1))
)


pVacc_plot <-(
  comp_traj_F %.>% 
    filter(., comp_exp == "Vaccinated") %.>%  
    plot_comp_summary(., y_label = "Percent Vaccinated", legend_position = "none", 
                      labels = scales::percent, limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1))
)


nExpo_plot <- (
  comp_traj_V %.>% 
    filter(., comp_exp == "Exposed") %.>% 
    plot_comp_summary(., 
                      y_label = expression(log[10](Exposed~Per~100000)), legend_position = "none", 
                      trans =  "log10",  breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2), 
                      limits = c(1e-3, 1e3))
)

nInfc_plot <- (
  comp_traj_V %.>% 
    filter(., comp_exp == "Infectious") %.>% 
    plot_comp_summary(., 
                      y_label = expression(log[10](Infectious~Per~100000)), legend_position = "none", 
                      trans =  "log10", breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2), 
                      limits = c(1e-3, 1e3))
)

  
nTrCs_plot <- (
  comp_traj_V %.>% 
    filter(., comp_exp == "True Cases") %.>% 
    plot_comp_summary(., 
                      y_label = expression(log[10](True~Cases~Per~100000)), legend_position = "none", 
                      trans =  "log10", breaks = c(1e-2, 1e-1, 1e0, 1e1, 1e2), 
                      limits = c(1e-3, 1e3))
)

         




# estimate the average transmission rate and the corresponding mean age at infection! - Fingers crossed!!!
append_rp_with_these <- c(N = 1e8, p = 1, nu = 1/80, ad = age_class_duration)

tic()
summarized_traj <- summarize_epidemiology(traj_cov_data = comp_traj_w_cov %.>% mutate(., `.id` = 1),
                                          p_vals = c(sim_from_these,
                                                append_rp_with_these))
toc()


#load("./summarized_traj.rds")

# find annual averages
treated_summarized_traj <- (
  summarized_traj %.>% 
    mutate(., year_val = floor(year)) %.>% 
    select(., year_val, Reff, average_beta, mean_age_at_infection) %.>% 
    group_by(., year_val) %.>% 
    summarise_all(., mean) %.>% 
    ungroup(.) %.>% 
    mutate(., mean_age_at_infection = ifelse(mean_age_at_infection > 80, 80, mean_age_at_infection)) %.>%
    gather(., "summary", "value", -year_val)
  ) 



average_beta_plot <- (
  treated_summarized_traj %.>% 
    filter(., summary == "average_beta") %.>% 
    ggplot(., aes(x = year_val, y = value)) +
    labs(x = "Year", 
         y = "Average Transmission Rate") +
    geom_line(size = 1, color = "grey30") +
    scale_x_continuous(breaks = year_break_x) +
    scale_y_continuous(breaks = c(80, 120, 160, 200, 240), limits = c(80, 240)) + 
    project_theme +
    cap_axes
    )

mean_age_inf_plot <- (
  treated_summarized_traj %.>% 
    filter(., summary == "mean_age_at_infection") %.>%
    ggplot(., aes(x = year_val, y = value)) +
    labs(x = "Year", 
         y = "Mean Age At Infection") +
    geom_line(size = 1, color = "grey30") +
    scale_x_continuous(breaks = year_break_x) +
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80)) +
    project_theme +
    cap_axes
)


Reff_plot <- (
  treated_summarized_traj %.>% 
  filter(., summary == "Reff") %.>% 
  ggplot(., aes(x = year_val, y = value)) +
  labs(x = "Year", 
       y = "Effective Reproductive Number") +
  geom_line(size = 1, color = "grey30") +
  geom_hline(yintercept = 1, size = 0.8, linetype = "dotdash", colour = "red")  +
  scale_x_continuous(breaks = year_break_x) +
  scale_y_continuous(breaks = c(0.92, 0.96, 1, 1.04, 1.08), limits = c(0.92, 1.08)) +
  project_theme +
  cap_axes
  )



F_plot_grid <- (
  plot_grid(pSusc_plot, pRecv_plot, pVacc_plot, nrow = 1,
            labels = c("A", "B", "C"), align = "hv")
)

V_plot_grid <- (
  plot_grid(nExpo_plot, nInfc_plot, nTrCs_plot, nrow = 1,
            labels = c("C", "D", "E"), align = "hv")
)


stat_plot_grid <- (
  plot_grid(average_beta_plot, mean_age_inf_plot, Reff_plot, nrow = 1,
          labels = c("F", "G", "H"), align = "hv")
  )



summary_plot_grid <- (
  plot_grid(F_plot_grid, V_plot_grid, 
            stat_plot_grid,   
            nrow = 3, align = "h", axis = "lr")
  )






  







