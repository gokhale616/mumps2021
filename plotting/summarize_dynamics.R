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
    annotate(geom = "segment", y = critic_lev_vacc, yend = critic_lev_vacc, x = 0, xend = 200, 
             size = 0.9, colour = "#6be585", linetype = "dashed") +
    geom_segment(aes(x = 0, xend = time, y = p_immunity_lost, yend = p_immunity_lost, 
                     colour = p_immunity_lost)) +
    geom_point(data = anno_segment, aes(x = x_c, y = y_c, colour = y_c), shape = 21, fill = "white", 
               size = 3) +
    geom_label(data = anno_segment, 
              aes(label = paste0("(", round(x_c, 1), " yrs, " , round(y_c, 2)*100, "%)"), 
                  x = x_c, y = y_c, fill = y_c), 
              nudge_x = 17, nudge_y = 0.065, colour = "white", size = 1.25) +
    annotate(geom = "text", y = 1, x = 72, 
             label= "Critical Vaccination Level",
             colour = "grey30", size = 2.5) +
    labs(x = "Time Since Immunization\n(Years)", 
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
    cap_axes(xlim = c(0,100), ylim = c(0.25, 1))
  )


# define a pomp obect to summarize the epidemiology
po_to_sum <- make_pomp(., 
                       covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1967-1/52, extra_end_t = fin_year, 
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
    labs(y = expression(Immunity~Lost~Per~10^5), x = "", 
         color = "Age\nCohort") +
    scale_y_continuous(labels = scales::scientific, 
                       limits = c(0, 1.5e1), 
                       #breaks = seq(0, 3e4,by = 1e4)
                       ) +
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
             y = critic_lev_vacc+0.065, x = 1983, 
             colour = "grey30", size = 2.5) +
    labs(x = "", y = "Percent Effectively   \nVaccinated    ") +
    scale_y_continuous(labels = scales::percent, 
                       limits = c(0, 1), 
                       breaks = seq(0, 1,by = 0.2)) +
    scale_x_continuous(breaks = year_break_x) +
    scale_colour_manual(values = c(orange_age_cohort, "#4b6cb7"), name = "Age\nCohort", 
                        breaks = c(age_names, "total")) +
    project_theme +
    cap_axes() +
    guides(colour = guide_legend(nrow = 2, title.position = "left"))
    ) 


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
    select(., `.id`, year, 
           starts_with("Ss_"), starts_with("Sp_"), 
           starts_with("Is_"), starts_with("Cs_"), Is, Reff)
  
  
  save(summarized_traj, file = paste0(path_dir, "/summarized_traj.rds"))
  
} else {
  message("Directory 'epi_summary' already exists, moving on!")
}


load(paste0(path_dir, "/summarized_traj.rds"))


Is_anno_data <- (
  summarized_traj %.>% 
    mutate(., year = floor(year)) %.>% 
    select(., year, Is) 
    )
  

Susc_plot <- (
  summarized_traj %.>% 
    select(., `.id`, year, starts_with("Sp_")) %.>% 
    prep_as_plot_trajectory(., init_year = 1967-1/52) %.>% 
    mutate(., count = ifelse(count > 1, 1, count)) %.>% 
    ggplot(., aes(x = year)) +
    geom_line(aes(y = count, colour = age_cohort), size = 0.8) +
    labs(y = "Percent Susceptible", 
         x = "Year") +
    scale_x_continuous(breaks = year_break_x) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_brewer(palette = "Oranges", direction = -1, guide = "none") +
    project_theme +
    cap_axes()
)

prevalence_plot <- (
  summarized_traj %.>% 
  select(., `.id`, year, starts_with("Is_")) %.>% 
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
  )
  

# also plot the vaccine uptake on the reffective plot
vacc_plot_data <- (
  mod_mumps_covariates_sigmoidal %.>% 
    bind_rows(., 
              mod_mumps_covariates_sigmoidal %.>% 
                slice(., rep(n(), interp_range[2]-2018))
              ) %.>% 
    select(., year, p1, p2) %.>% 
    filter(., 
           year > 1966) %.>% 
    mutate(., year =  seq(1967, 2020, by = 1)) %.>% 
    gather(., key = "dose", value = "coverage", -year)
  )



Reff_plot <- (
  summarized_traj %.>% 
    filter(., year > 1967-1/52) %.>% 
    select(., year, Reff) %.>% 
    mutate(., gt1 = case_when(Reff > 1 ~ "yes!", 
                              Reff < 1 ~ "nein!")
    ) %.>% 
    ggplot(.) +
    geom_bar(data = vacc_plot_data, 
             aes(x = year, y = coverage*2, fill = dose), 
             position="dodge", stat="identity", 
             alpha = 0.5) +
    geom_line(aes(x = year, y = Reff, colour = gt1, group = 1), size = 0.8) +
    labs(x = "", y = "Effective Reproductive    \nNumber    ") +
    scale_colour_manual(values = c("yes!" = "#f64f59", "nein!" = "darkseagreen4"), ##009FFF
                        labels = c("yes!" = "Super-critical", "nein!" = "Sub-critical"), 
                        name = "Epidemic\nSignature") +
    scale_y_continuous(sec.axis = sec_axis(~./2,
                                           labels = scales::percent,
                                           name = "Vaccine\nCoverage")) +
    scale_x_continuous(breaks = year_break_x) +
    scale_fill_manual(name = "Dose\nType", 
                      values = c("grey50", "black"), 
                      labels = c("Neonatal", "Booster")) +
    project_theme +
    coord_capped_cart(bottom='both', left = "both",  right='both') +
    guides(colour = guide_legend(nrow = 2), 
           fill = guide_legend(nrow = 2, 
                               override.aes = list(alpha = 1)))
)



epi_summary_plt <- (
  immune_distbn_plot + 
    V_prop_plot + Susc_plot +
    Reff_plot + Vs_plot + 
    prevalence_plot +
    plot_layout(
      guides = "collect",
      design = "
    AD
    BE
    CF
    "
    ) + 
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", 
          legend.key.size = unit(0.75, "lines"), 
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10))
)




reporting_probs <- best_model_p_vec[sprintf("rho_age_%s", c(1:5, "u"))] 

# calculate the observed incidence per 1e5 for some years - here, the unknown age class ignored 
obs_age_incidence <- (
  mumps_case_reports_l %.>%
    select(., year, age_cohort, cases) %.>% 
    filter(., age_cohort != "unknonwn") %.>% 
    mutate(., tcases = case_when(age_cohort == age_names[1] ~ cases/reporting_probs[1], 
                                 age_cohort == age_names[2] ~ cases/reporting_probs[2], 
                                 age_cohort == age_names[3] ~ cases/reporting_probs[3], 
                                 age_cohort == age_names[4] ~ cases/reporting_probs[4], 
                                 age_cohort == age_names[5] ~ cases/reporting_probs[5])) %.>% 
    right_join(., 
               mumps_covariates %.>% 
                 filter(., year > 1975) %.>% 
                 select(., year, starts_with("N_")) %.>% 
                 mutate(., `.id` = 1) %.>% 
                 prep_as_plot_trajectory(., init_year = 1965) %.>% 
                 transmute(., year = year, age_cohort = age_cohort, popn_factor = 1e5/count),
                 by = c("year", "age_cohort")) %.>% 
    transmute(., year = year, age_cohort = age_cohort, inc = tcases*popn_factor) %.>% 
    group_by(., year) %.>% 
    mutate(., obs_prop_inc = inc/sum(inc))
  )


expctd_age_incidence <- (
  summarized_traj %.>% 
    select(., `.id`, year, starts_with("Cs_")) %.>% 
    prep_as_plot_trajectory(., init_year = 1965-1/52) %.>% 
    mutate(.,
           in_sample = ifelse(year<2012+20/52, "ja", "nein"),
           cmplx_grdnt = paste0(age_cohort, "_", in_sample) %.>% 
             factor(., levels = c(paste0(age_names, "_", "ja"), 
                                  paste0(age_names, "_", "nein")
                                  )
                    )
           )
  )



# simulate from the observation process 
simulated_annual_case_data <- (
  sim_obs_model(po_obj = po_to_sum, params = best_model_p_vec, 
                times = time(po_to_sum), 
                nsim = 1000, 
                root_transform = FALSE, 
                summarise_obs_process = FALSE) %.>% 
    mutate(., 
           year = floor(year)) %.>%   
    filter(., year > 1976 & year < 2020) %.>% 
    # calculate the annual number of cases recorded
    group_by(., year, `.id`, age_class) %.>% 
    summarise(., cases = sum(cases)) %.>% 
    ungroup(.)
)


# calculate true incidence per 10^5
simulated_annual_inc_data <- (
  simulated_annual_case_data %.>%
    # add scaling factor to calculate true annual cases 
    right_join(., 
               mumps_covariates %.>% 
                 filter(., year > 1976 & year < 2020) %.>% 
                 select(., year,eta_a),
               by = "year") %.>% 
    mutate(., 
           rp_sf = case_when(age_class == age_names_u[1] ~ eta_a*reporting_probs[1], 
                             age_class == age_names_u[2] ~ eta_a*reporting_probs[2], 
                             age_class == age_names_u[3] ~ eta_a*reporting_probs[3], 
                             age_class == age_names_u[4] ~ eta_a*reporting_probs[4], 
                             age_class == age_names_u[5] ~ eta_a*reporting_probs[5], 
                             age_class == age_names_u[6] ~ (1-eta_a)*reporting_probs[6]
           ),
           true_cases = cases/rp_sf
    ) %.>% 
    # ignoring the "unknown" age class from here on
    filter(., age_class != "unknown") %.>% 
    # add population sizes to calculate and calculate true incidence
    right_join(., 
               mumps_covariates %.>% 
                 mutate(., `.id` = 1) %.>% 
                 select(., year, `.id`, starts_with("N_")) %.>% 
                 prep_as_plot_trajectory(., init_year = 1977-1/52) %.>% 
                 transmute(., 
                           year = year, age_class = age_cohort, 
                           popn = count, popn_sf = 1e5/count),
               by = c("year", "age_class")) %.>% 
    mutate(., 
           sim_inc = true_cases*popn_sf) %.>% 
    # calculate incidence age distribution
    group_by(., year, `.id`) %.>% 
    mutate(., 
           sim_prop_inc = sim_inc/sum(sim_inc)) %.>% 
    ungroup(.)
)


# make one data with observed v/s simlulated proportion
obs_sim_annual_inc_data <- (
  simulated_annual_inc_data %.>%   
    right_join(., 
               by = c("year", "age_class"),
               mumps_case_reports_l %.>%
                 # ignoring the "unknown" age class
                 filter(., age_cohort != "unknown") %.>% 
                 transmute(., 
                           year = year, age_class = age_cohort, 
                           obs_cases = cases)
    ) %.>% 
    # calculate the true observed incidence
    mutate(., 
           true_obs_cases = obs_cases/rp_sf, 
           obs_inc = true_obs_cases*popn_sf) %.>% 
    # calculate incidence distribution 
    group_by(., year, `.id`) %.>% 
    mutate(., 
           obs_prop_inc = obs_inc/sum(obs_inc)) %.>% 
    ungroup(.)
)


# calculate the mean age of infection
# make data frame with age classes
age_df <- (
  expand_grid(year = seq(1977, 2018, by = 1), 
              age = seq(0, 80, by = 1)) %.>% 
    mutate(., 
           age_class = case_when(age < 5  ~ "[0,5)", age < 16 ~ "[5,15)", 
                                 age < 25 ~ "[15,25)", age < 40 ~ "[25,40)", 
                                 TRUE ~ ">40") %.>% 
             as_factor(.)
           )
    )
  

# extract weights for the simulation
mean_age_data <- (
  obs_sim_annual_inc_data %.>% 
    select(., 
           year, age_class, sim_prop_inc) %.>% 
    right_join(., 
               obs_age_incidence %.>% 
                 transmute(., 
                           year = year, age_class = age_cohort, obs_prop_inc = obs_prop_inc), 
               by = c("year", "age_class")) %.>% 
    right_join(., 
               age_df, 
               by= c("year", "age_class")) %.>% 
    group_by(., year) %.>% 
    mutate(., 
           exp_mean_age_of_infection = sum(age*sim_prop_inc)/sum(sim_prop_inc), 
           obs_mean_age_of_infection = sum(age*obs_prop_inc)/sum(obs_prop_inc)) %.>% 
    ungroup(.) %.>%   
    select(., year, exp_mean_age_of_infection, obs_mean_age_of_infection) %.>% 
    distinct(.) %.>% 
    gather(., "stat", "age", -year)
)

# compare age distribution of observed vs expected 
  
wis_fill_gradnt <- c(model_col[1], "#3EABA2", "#6BBDB6", "#98CECA", "#C5E0DE")

oos_fill_gradnt <- c(model_col[1], "#FE8C00", "#FCA030", "#F9B561","#F7C991", "#F4DEC2")

grey30_gradnt <- c("grey30", "#6E6E6E", "#8F8F8F", "#B0B0B0", "#D1D1D1")

age_distribution_plot <- (
  expctd_age_incidence %.>% 
    filter(., year > 1977-20/52 & year < 2018+20/52) %.>% 
    ggplot(., aes(x = year, y = count)) +
    geom_area(aes(fill = cmplx_grdnt), stat = "identity", 
              position = position_fill(reverse = TRUE)) +
    geom_bar(data = obs_age_incidence, 
             aes(x = year, y = inc, colour = age_cohort), 
             stat = "identity", position = position_fill(reverse = TRUE), fill = NA, size = 0.8, width =0.75) +
    geom_line(data = mean_age_data, 
              aes(x = year, y = age/100, linetype = stat),
              size = 1.0, colour = "grey30", alpha = 0.75) +
    labs(x = "", y = "Observed v/s Expected \nIncidence Age Distribution", 
         colour = "Age\nCohort", 
         pattern_fill = "Age\nCohort")+
    scale_y_continuous(labels = scales::percent, 
                       sec.axis = sec_axis(~.*100, breaks = seq(0, 100, 25), 
                                           name = "Age (Years)")) +
    scale_x_continuous(breaks = c(1977,1984, 1991, 1998, 2005, 2012, 2018)) +
    scale_fill_manual(values = c(wis_fill_gradnt, oos_fill_gradnt), 
                      breaks = c(paste0(age_names, "_", "ja"), 
                                 paste0(age_names, "_", "nein")),
                      guide = "none") +    
    scale_colour_brewer(palette = "Purples", direction = -1) +
    scale_linetype_manual(values = c(1, 2), labels = c("Expected", "Observed"), 
                          name = "Mean Age Of\nFirst Infection") +
    project_theme +
    theme(legend.position = "top") +
    cap_axes(right = "both") +
    theme() +
    guides(colour = guide_legend(title.position = "top", nrow = 2, 
                                 override.aes=list(fill = grey30_gradnt)), 
           linetype = guide_legend(title.position = "top", nrow = 2, 
                                   override.aes = list(alpha = 1)))
  
  
)



# check if the Kl_summary exits if not make it 
path_dirKL <- "../result_data/KL_summary"

if(dir.exists(path_dirKL) == FALSE) {
  
  dir.create(path_dirKL)
  message("Directory 'KL_summary' has been created, proceeding to populate!")

  
  # prepare data for plotting
  obs_sim_annual_KL_data <- (
    obs_sim_annual_inc_data %.>% 
      select(., 
             year, `.id`, age_class, 
             sim_prop_inc, obs_prop_inc)
  )
  
  
  cal_age_distn_KL <- function(y, i) {
    
    df <- obs_sim_annual_KL_data
    year_vec <- df$year %.>% unique(.)
    # browser()
    sim_distn <- (
      df %.>% 
        filter(., year == year_vec[y] &`.id` == i) %.>% 
        select(., age_class, sim_prop_inc) %.>%
        spread(., key = age_class, value = sim_prop_inc) %.>% 
        unlist(.)
      )
  
    obs_distn <- (
      df %.>% 
        filter(., year == year_vec[y] &`.id` == i) %.>% 
        select(., age_class, obs_prop_inc) %.>%
        spread(., key = age_class, value = obs_prop_inc) %.>% 
        unlist(.)
      )
    
    
    if(is.na(obs_distn) == FALSE) {
      accomodate_NA <- FALSE
      KLD_res <- KLD(sim_distn, obs_distn, base = 2)  
    } else {
      accomodate_NA <- TRUE
    }
    
    # browser()
    
    tibble(year = year_vec[y],
           `.id` = i,
           ms_KLD = ifelse(accomodate_NA == TRUE, NA, KLD_res$mean.sum.KLD), 
    )
    
  
  }


  # calculate KL div distribution
  KL_distbn <- map_dfr(1:42, function(y) {map_dfr(1:1000, function(i) {cal_age_distn_KL(y, i)})})

  save(KL_distbn, file = paste0(path_dirKL, "/KL_distbn.rds"))
  
} else {
  message("Directory 'KL_summary' already exists, moving on!")
}


load(paste0(path_dirKL, "/KL_distbn.rds"))


obs_inc_data <- (
  obs_sim_annual_inc_data %.>% 
    filter(., `.id` == 1) %.>% 
    group_by(., year) %.>%   
    summarise(., 
              sqrt_tot_inc = ((sum(true_obs_cases)/sum(popn))*1e5) %.>% sqrt(.))  
)
  

# this is a supplementary plot looking at the discrepancy between observed and expected age distribution.
KL_divergence_plot <- (
  KL_distbn %.>% 
    mutate(., in_or_out  = ifelse(year <= 2012, "in", "out")) %.>% 
    ggplot(., aes(x = year)) +
    geom_boxplot(aes(y = ms_KLD, group = year, colour = in_or_out),
                 outlier.shape = 21, outlier.fill = "white") +
    geom_area(data = obs_inc_data, aes(y = sqrt_tot_inc/100), fill = "#734b6d", alpha = .2) +
    labs(x = "Year", y = "Kullback-Leibler\nDivergence") +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1), 
                       sec.axis = sec_axis(~.*100, breaks = seq(0, 50, 10), 
                                           name = expression(sqrt(Incidence~Per~10^5)))) +
    scale_x_continuous(breaks = c(1977,1984, 1991, 1998, 2005, 2012, 2018)) +
    scale_colour_manual(values = model_col, 
                      labels = c("Within\nSample", "Out of\nSample"), 
                      name = "Prediction\nEpoch") +
    project_theme +
    coord_capped_cart(bottom='both', left = "both",  right='both') +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title.position = "top", nrow = 2))
  
)



compare_age_dstbn_plt <- (
  age_distribution_plot + KL_divergence_plot + 
    plot_layout(
      guides = "collect",
      design = "
    A
    B
    "
    ) + 
    plot_annotation(tag_levels = "A") &
    theme(legend.position='bottom') 
)

  

