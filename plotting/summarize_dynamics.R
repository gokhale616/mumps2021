# get some values ready for this plot 
# simulate dynamics to the check simulation match to theoretical results 
fin_year <- 2020

# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
  select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
  mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
  unlist(.) %.>% sim_p_vals(.)
  )


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



# test how the pdf looks for different values 
time_tibble <- tibble(time = seq(0, 700, by = 1/1e3))


prob_exp <- function(time, rate) {
  exp(-rate*time)
}


log_t <- function(prob,rate) {
  -log(prob)/rate
}

time_tibble %<>% mutate(., p_immunity_lost = prob_exp(time, rate = (1/best_model$dwan)))

# base on the the estimated R0 what is the critical level of vaccination?
critic_lev_vacc <- 1-1/best_model$R0

t_critic_lev_vacc_reached <- log_t(prob = critic_lev_vacc, rate = (1/best_model$dwan))

anno_segment <- (
  time_tibble %.>% 
    filter(., time %in% c(25, 50, 75)) %.>% 
    transmute(., 
              x_c = time, 
              y_c = p_immunity_lost, zero = 0) %.>% 
    bind_rows(., 
              tibble(x_c = t_critic_lev_vacc_reached, 
                     y_c = critic_lev_vacc, 
                     zero = 0)) %.>% 
    arrange(., x_c)
  )


y_low_lim <- prob_exp(100, rate = (1/best_model$dwan))


# this is the plot of immune distribution taking values from the best fitting model

immune_distbn_plot <- (
  time_tibble %.>% 
    ggplot(.) +
    geom_segment(aes(x = 0, xend = time, y = p_immunity_lost, yend = p_immunity_lost, 
                     colour = p_immunity_lost)) +
    geom_point(data = anno_segment, aes(x = x_c, y = y_c, colour = y_c), shape = 21, fill = "white", 
               size = 3) +
    geom_label(data = anno_segment, 
              aes(label = paste0("(", round(x_c, 1), " yrs, " , round(y_c, 2)*100, "%)"), 
                  x = x_c, y = y_c, fill = y_c), 
              nudge_x = 13.5, nudge_y = 0.065, colour = "white", size = 1.75) +
    annotate(geom = "segment", y = 1, yend = 1, x = 42.5, xend = 35.5, 
            arrow = arrow(length = unit(1.75, "mm")), colour = "grey30", size = 0.5) +
    annotate(geom = "text", y = 1, x = 70, 
             label= "Critical Vaccination Level",
             colour = "grey30", size = 2.5) +
    labs(x = "Time Since Immunization (Years)", 
         y = "Percent Immune\nPost Vaccination") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(round(0, 1), 1, by = 0.25), 
                       labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, 100, by = 25)) +
    scale_colour_gradient(low = "grey30", high = "#6be585", 
                          limits = c(0, 1), 
                          guide = "none") +
    scale_fill_gradient(low = "grey30", high = "#6be585",
                          limits = c(0, 1), 
                          guide = "none") +
    project_theme +
    cap_axes(xlim = c(0,100), expand = FALSE)
  )


# define a pomp obect to summarize the epidemiology
po_to_sum <- make_pomp(., 
                       covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1965-1/52, extra_end_t = fin_year, 
                       temp_scale = 1/52)

# Generate trajectories from all the compartments 
comp_traj <- (
  trajectory(po_to_sum, 
             param = best_model_p_vec, 
             method = "ode23", format = "d") 
  )


comp_traj_w_cov <- (
  comp_traj %.>% 
    right_join(., interpolated_covs, by = "year") %.>% 
    as_tibble(.) %.>% 
    mutate(., `.id` = as.numeric(`.id`)) %.>% 
    drop_na(.)
)

##############################################################################################################  
# from here on we will be looking at what is happening at the vaccinated population and how they lose immunity
##############################################################################################################
V_prop_traj <- (
  comp_traj_w_cov %>% 
    mutate(., 
           Vs = Vs_1 + Vs_2 + Vs_3 + Vs_4 + Vs_5, 
            V = V_1 + V_2 + V_3 + V_4 + V_5,
            N = N_1 + N_2 + N_3 + N_4 + N_5, 
           Vprop = (V)/N, 
           Vprop_1 = (V_1)/N_1,
           Vprop_2 = (V_2)/N_2, Vprop_3 = (V_3)/N_3,
           Vprop_4 = (V_4)/N_4, Vprop_5 = (V_5)/N_5,
           Vsprop_1 = (Vs_1)*1e5/N_1,
           Vsprop_2 = (Vs_2)*1e5/N_2, Vsprop_3 = (Vs_3)*1e5/N_3,
           Vsprop_4 = (Vs_4)*1e5/N_4, Vsprop_5 = (Vs_5)*1e5/N_5,
           ) %.>% 
    select(., `.id`, year, starts_with("Vprop"), starts_with("Vsprop"))
           
)

V_prop_traj_for_plot <- (
  V_prop_traj %.>% 
    select(., -Vprop) %.>% 
    prep_as_plot_trajectory(., init_year = 1965-1/52) 
  
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
    labs(y = expression(Immunity~Lost~Per~10^5~phantom(1000)), x = "Year", 
         color = "Age\nCohort") +
    scale_y_continuous(labels = scales::scientific, 
                       limits = c(0, 3e4), 
                       breaks = seq(0, 3e4,by = 1e4)
                       ) +
    scale_x_continuous(breaks = year_break_x) +
    scale_colour_brewer(palette = "Oranges", direction = -1) +    
    project_theme +
    cap_axes() +
    theme(legend.position = c(0.33, 0.55)) +
    guides(colour = guide_legend(nrow = 3))
  )

  


V_prop_plot <- (
  V_prop_traj_for_plot %.>%
    filter(., comp_exp == "Vprop") %.>% 
    ggplot(.) +
    geom_area(aes(x = year, y = count, 
                  fill = age_cohort), position = position_dodge(width = 0), alpha = 0.7) +
    geom_line(data = anno_Vs_tot_prop_data, 
              aes(x = year, y = Vprop, colour = "total"), size = 1.2) +
    annotate(geom = "segment", 
             x = 1965, xend = 2020, 
             y = critic_lev_vacc, yend = critic_lev_vacc,
             colour = "#6be585", linetype = "dotdash", size = 1.2) +
    annotate(geom = "text", 
             label = "Critical Vaccination Level",
             y = (critic_lev_vacc)-0.1, x = 1978, 
             colour = "grey30", size = 2.5) +
    labs(x = "", y = "Percent Effectively           \nVaccinated           ", fill = "Age\nCohort") +
    scale_y_continuous(labels = scales::percent, 
                       limits = c(0, 1), 
                       breaks = seq(0, 1,by = 0.2)) +
    scale_x_continuous(breaks = year_break_x) +
    scale_fill_brewer(palette = "Oranges", direction = -1, guide = "none") +    
    scale_colour_manual(values = ("total" = "#FFE000"), labels = "Total", name= "") +    
    project_theme +
    theme(legend.position = c(0.15, 0.45), 
          legend.spacing.y = unit(0.1, "cm")
          ) +
    cap_axes() 
    ) 


Vs_prop_plot <- plot_grid(V_prop_plot, Vs_plot, nrow = 2, rel_heights = c(1, 1), 
                          labels = c("B", "C"), align = "v", axis = "lb")


##############################################################################################################
############# now we look at what is happening to the effective R0 and the infectious individuals ############
##############################################################################################################

# check if the epi_summary exits if not make it 

path_dir <- "../result_data/epi_summary"

if(dir.exists(path_dir) == FALSE) {
  dir.create(path_dir)
  message("Directory 'epi_summary' has been created, proceeding to populate!")
  
  # generate summary measures for for compartments 
  
  summarized_traj <- summarize_epidemiology(traj_cov_data = comp_traj_w_cov,
                                            p_vals = c(best_model_p_vec,
                                                       params_for_Rp)) %.>% 
    select(., `.id`, year, starts_with("Is_"), Is, Reff)
  
  
  save(summarized_traj, file = paste0(path_dir, "/summarized_traj.rds"))
  
} else {
  message("Directory 'epi_summary' already exists, moving on!")
}


load(paste0(path_dir, "/summarized_traj.rds"))


Is_anno_data <- (
  summarized_traj %.>% 
    mutate(., year = floor(year)) %.>% 
    select(., year, Is) #%.>% 
    # group_by(., year) %.>% 
    # summarize_all(., sum) %.>% 
    # ungroup(.)
    )
  

prevalence_plot <- (
  summarized_traj %.>% 
  select(., `.id`, year, starts_with("Is_")) %.>% 
  prep_as_plot_trajectory(., init_year = 1965-1/52) %.>% 
  ggplot(., aes(x = year)) +
  geom_line(aes(y = count, colour = age_cohort), size = 0.8) +
  labs(y = expression(Infectious~per~10^5), 
       x = "Year") +
  scale_x_continuous(breaks = year_break_x) +
  scale_y_continuous(trans = "log10", breaks = c(1e-3, 1e-2, 1e-1, 1e0, 
                                                 1e1, 1e2, 1e3), 
                     limits = c(1e-3, 1e3)) +
  annotation_logticks(sides = "l")  +
  scale_colour_brewer(palette = "Oranges", direction = -1, guide = "none") +
  project_theme+
  cap_axes(expand = FALSE)
  )
  
  
  
  
Reff_plot <- (
  summarized_traj %.>% 
    select(., year, Reff) %.>% 
    mutate(., gt1 = case_when(Reff > 1 ~ "yes!", 
                              Reff < 1 ~ "nein!")
           ) %.>% 
    ggplot(., aes(x = year, y = Reff)) +
    geom_line(aes(colour = gt1, group = 1), size = 1) +
    labs(x = "Year", y = "Effective Reproductive    \nNumber    ") +
    scale_colour_manual(values = c("yes!" = "#f64f59", "nein!" = "darkseagreen4"), ##009FFF
                        labels = c("yes!" = "Super Critical", "nein!" = "Sub-critical"), 
                        name = "") +
    scale_y_continuous(limits = c(0.90, 1.10), breaks = c(0.90, 0.95, 1, 1.05, 1.10)) +
    scale_x_continuous(breaks = year_break_x) +
    project_theme +
    cap_axes() +
    theme(legend.position = c(0.35, 0.8)) +
    guides(colour = guide_legend(nrow = 2))
)
    







immune_summary_plot <- plot_grid(immune_distbn_plot, Vs_prop_plot, nrow = 2, 
                                 align = "h", axis = "l", rel_heights = c(0.5, 1), 
                                 labels = c("A", ""))



prev_summary_plot <- plot_grid(Reff_plot, prevalence_plot, nrow = 2, 
                               align = "v", rel_heights = c(0.5, 1), 
                               labels = c("D", "E"))



epi_summary_plt <- (
  plot_grid(immune_summary_plot, prev_summary_plot, nrow = 1, align = "hv", axis = "bl"))










