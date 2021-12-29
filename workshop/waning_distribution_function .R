source("../00/src.R", chdir = TRUE)
source("../plotting/prep_estm_tables.R", chdir = TRUE)

# get some values ready for this plot 
# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
  select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
  mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
  unlist(.) %.>% sim_p_vals(.)
  )


# test how the pdf looks for different values 
time_tibble <- tibble(time = seq(0, 300, by = 1/1e3))


prob_exp <- function(time, rate) {
  exp(-rate*time)
}

log_t <- function(prob,rate) {
  -log(prob)/rate
}

time_tibble %<>% mutate(., p_immunity_lost = prob_exp(time, rate = (1/best_model$dwan+1/80)))

# base on the the estimated R0 what is the critical level of vaccination?
critic_lev_vacc <- 1-1/best_model$R0

t_critic_lev_vacc_reached <- log_t(prob = critic_lev_vacc, rate = (1/best_model$dwan+1/80))

anno_segment <- (
  time_tibble %.>% 
    filter(., time %in% c(5, 10, 20, 40, 80)) %.>% 
    transmute(., 
              x_c = time, 
              y_c = p_immunity_lost, zero = 0) %.>% 
    bind_rows(., 
              tibble(x_c = t_critic_lev_vacc_reached, 
                     y_c = critic_lev_vacc, 
                     zero = 0)) %.>% 
    arrange(., x_c)
  )


y_low_lim <- prob_exp(100, rate = (1/best_model$dwan+1/80))


# this is the plot of immune distribution taking values from the best fitting model

immune_distbn_plot <- (
  time_tibble %.>% 
    ggplot(.) +
    geom_segment(aes(x = 0, xend = time, y = p_immunity_lost, yend = p_immunity_lost, 
                     colour = p_immunity_lost)) +
    geom_point(data = anno_segment, aes(x = x_c, y = y_c, colour = y_c), shape = 21, fill = "white", 
               size = 4) +
    geom_label(data = anno_segment, 
              aes(label = paste0("(", round(x_c, 1), " yrs, " , round(y_c, 2)*100, "%)"), 
                  x = x_c, y = y_c, fill = y_c), 
              nudge_x = 6.5, nudge_y = 0.015, colour = "white") +
    annotate(geom = "segment", y = 0.945, yend = 0.945, x = 20.5, xend = 16.5, 
             arrow = arrow(length = unit(2, "mm")), colour = "grey30") +
    annotate(geom = "text", y = 0.945, x = 30, label= "Critical vaccination level", 
             arrow = arrow(length = unit(2, "mm")), colour = "grey30") +
    labs(x = "Time Since Immunization (Years)", 
         y = "Percent Immune\nPost Vaccination") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(round(0, 1), 1, by = 0.25), 
                       labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, 100, by = 25)) +
    scale_colour_gradient(low = "grey30", high = "#6be585", 
                          limits = c(0, 1), 
                          guide = FALSE) +
    scale_fill_gradient(low = "grey30", high = "#6be585",
                          limits = c(0, 1), 
                          guide = FALSE) +
    project_theme +
    cap_axes(xlim = c(0,100))
  )


# simulate dynamics to the check simulation match to theoretical results 
fin_year <- 2020

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


##############################################################################################################  
# from here on we will be looking at what is happening at the vaccinated population and how they lose immunity
##############################################################################################################
V_traj <- (
  comp_traj %.>% 
    mutate(., 
          `.id` = as.numeric(`.id`),
           year = floor(year)) %.>% 
    group_by(., year) %.>% 
    summarise_all(., sum) %.>% 
    ungroup(.) %.>% 
    select(., year, `.id`, starts_with("V_")) %.>% 
    prep_as_plot_trajectory(., init_year = 1965-1/52) %.>% 
    mutate(., V_count = count) %.>% 
    select(., -c(comp_exp, count)) %.>% 
    group_by(., year) %.>% 
    mutate(., V_tot_count = sum(V_count)) %.>% 
    ungroup(.)
)



Vs_traj <- (
  comp_traj %.>% 
    mutate(., 
           `.id` = as.numeric(`.id`),
           year = floor(year)) %.>% 
    group_by(., year) %.>% 
    summarise_all(., sum) %.>% 
    ungroup(.) %.>% 
    select(., year, `.id`, starts_with("Vs_")) %.>% 
    prep_as_plot_trajectory(., init_year = 1965-1/52) %.>% 
    mutate(., Vs_count = count) %.>% 
    select(., -c(comp_exp, count)) %.>% 
    group_by(., year) %.>% 
    mutate(., Vs_tot_count = sum(Vs_count)) %.>% 
    ungroup(.)  
)


Vs_prop_traj <- (
  right_join(V_traj, Vs_traj, by = c("year", "age_cohort")) %.>% 
  mutate(., 
         prop_newly_susc  = ifelse(V_count == 0,0, Vs_count/V_count), 
         prop_newly_susc_tot  = ifelse(V_tot_count == 0,0, Vs_tot_count/V_tot_count))
  )

year_break_x <- seq(1970, fin_year, 15)

anno_Vs_tot_prop_data <- (
  Vs_prop_traj %.>% 
    filter(., age_cohort == ">40") %.>% 
    select(., year, prop_newly_susc_tot)
    )

Vs_prop_plot <- (
  Vs_prop_traj %.>%
    ggplot(.) +
    geom_area(aes(x = year, y = prop_newly_susc, 
                  fill = age_cohort), position = position_dodge(width = 0), alpha = 0.9) +
    geom_line(data = anno_Vs_tot_prop_data, 
              aes(x = year, y = prop_newly_susc_tot, colour = "total")) +
    geom_point(data = anno_Vs_tot_prop_data, 
              aes(x = year, y = prop_newly_susc_tot, colour = "total"), shape = 21, fill = "white") +
    annotate(geom = "rect", xmin = 1965, xmax = 2020, ymin = 0, ymax = 1-critic_lev_vacc, 
             fill = "#6be585", alpha = 0.4) +
    annotate(geom = "text", label = "Permissible Immune Loss Level", 
             y = (1-critic_lev_vacc)+0.005, x = 1972, 
             colour = "grey30") +
    labs(x = "Year", y = "Percent Vaccine\nImmunity Lost", fill = "Age\nCohort") +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = year_break_x) +
    scale_fill_brewer(palette = "Oranges", direction = -1) +    
    scale_colour_manual(values = ("total" = "#F9D423"), labels = "Total", name= "") +    
    project_theme +
    theme(legend.position = c(0.2, 0.8)) +
    cap_axes() +
    guides(fill = guide_legend(nrow = 2, override.aes = list(alpha = 1), order = 2), 
           colour = guide_legend(order = 1))
    ) 

##############################################################################################################
############# now we look at what is happening to the effective R0 and the infectious individuals ############
##############################################################################################################

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

comp_traj_w_cov <- (
  comp_traj %.>% 
    right_join(., interpolated_covs, by = "year") %.>% 
    as_tibble(.) %.>% 
    mutate(., `.id` = as.numeric(`.id`)) %.>% 
    drop_na(.)
  )


# generate summary measures for for compartments 
tic()
summarized_traj <- summarize_epidemiology(traj_cov_data = comp_traj_w_cov,
                                          p_vals = c(best_model_p_vec,
                                                     params_for_Rp)) %.>% 
  select(., `.id`, year, starts_with("Is_"), Is, Reff)
toc()


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
  labs(y = "Infectious per 10,0000", 
       x = "Year") +
  scale_x_continuous(breaks = year_break_x) +
  scale_y_continuous(trans = "log10") +
  scale_colour_brewer(palette = "Oranges", direction = -1, guide = FALSE) +
  project_theme+
  cap_axes()
  )
  
  
  
  
Reff_plot <- (
  summarized_traj %.>% 
    select(., year, Reff) %.>% 
    mutate(., gt1 = case_when(Reff > 1 ~ "yes!", 
                              Reff < 1 ~ "nein!",
                              Reff == 1 ~ "uno!")
           ) %.>% 
    ggplot(., aes(x = year, y = Reff)) +
    geom_line(aes(colour = gt1, group = 1), size = 1.2) +
    labs(x = "Year", y = "Effective Reproductive\nNumber") +
    scale_colour_manual(values = c("yes!" = "#f64f59", "nein!" = "#009FFF", "uno!" = "grey30"), 
                        labels = c("yes!" = ">1", "nein!" = "<1", "uno!" = "1"), 
                        name = "") +
    scale_y_continuous(limits = c(0.90, 1.10), breaks = c(0.90, 0.95, 1, 1.05, 1.10)) +
    scale_x_continuous(breaks = year_break_x) +
    project_theme +
    cap_axes() +
    theme(legend.position = c(0.2, 0.8)) +
    guides(colour = guide_legend(nrow = 2))
)
    















plot_grid(immune_distbn_plot, Reff_plot, Vs_prop_plot, prevalence_plot, nrow = 2, 
          align = "v", rel_heights = c(0.65, 0.65, 1, 1), labels = c("A", "B", "C", "D"))









if(FALSE) {
# lets see how the next generation matrix looks

eigen_analysis <- calculate_R0_mq(p_vals = c(best_model_p_vec, params_for_R0))


ngm <- eigen_analysis$K_mat[1:5, 6:10] 
dimnames(ngm) <- list(age_names, age_names)

ngm_plot <- (
  ngm %.>% 
    plot_contact_matrix(., col_max = "#ec2F4B",
                        limits = c(0, 15), 
                        breaks = seq(0, 15, 5), 
                        fill_lab = "Average Number\nOf New Cases") +
    project_theme + 
    cap_axes +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal", 
                                 barheight = 0.5))
  )

eigen_vector_tibble <- (
  eigen_analysis$other_output$vectors[1:5, 1] %.>% 
    tibble(.) %.>%
    transmute(., age_class = age_names, reproductive_contri = .)
  )
  
eigen_vector_tibble %.>%   
  ggplot(., aes(x = age_class, y = reproductive_contri)) +
  geom_bar(stat = "identity", fill = "#ec2F4B" , width = 0.45) +
  labs(x = "Age Cohort", 
       y = "Proportion of Cases") +
  project_theme +
  cap_axes

}