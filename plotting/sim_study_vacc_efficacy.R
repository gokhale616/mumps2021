# requires 
# source("./00/src.R", chdir = TRUE)
# source("./plotting/inc_plot.R", chdir = TRUE)
# source("./plotting/prep_estm_tables.R", chdir = TRUE)
# source("./plotting/fit_plot.R", chdir = TRUE)
# source("./plotting/rel_fit_plot.R", chdir = TRUE)
# source("./plotting/summarize_dynamics.R", chdir = TRUE)



# prepare the vector of the values used to in this sim study
best_model_p_vec_vacc_eff <- (
  all_result_df %.>% 
    filter(., best_fit_covar == 1) %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar, 
                 epsilon2, p_intro, t_intro, impact)) %.>% 
    unlist(.) %.>% sim_p_vals(., default_p_vals = param_vals_vacc_effective)
)


# simulate dynamics

# define a pomp object for this simulation study
# with booster
po_vacc_eff_with_booster <- (
  make_pomp_vacc_effective(., 
  covar = mod_mumps_covariates_sigmoidal, 
  temp_scale = 1/1, 
  extra_start_t = 2019,
  extra_end_t = 2300)
  )


# test_vacc_traj <- (
#   po_vacc_eff_with_booster  %.>%  
#     trajectory(., 
#                param = best_model_p_vec_vacc_eff, 
#                method = "ode23", format = "d") 
# )


pdwan_18 <- 1-exp(-18/best_model_p_vec_vacc_eff["dwan"]) %.>% unname(.)

dwan_grid <- tibble(dwan = c(seq(1, 200, by = 1), best_model_p_vec_vacc_eff["dwan"])) %.>% arrange(., dwan)


vacc_eff_grid <- function(x = 11, dwan_data = dwan_grid, 
                          po_obj
                          ) {
  
  # browser()
  # collect the control parameter value
  dwan_int <- dwan_grid %.>%  slice(., x) %.>% unlist(.) %.>% unname(.)
  
  # replace in the defaults
  p_for_sim <- best_model_p_vec_vacc_eff
  p_for_sim["dwan"] <- dwan_int
  
  # Generate trajectories from all the compartments 
  comp_vacc_eff_traj <- (
    trajectory(po_obj, 
               param = p_for_sim, 
               method = "ode23", format = "d") %.>% 
    transmute(., 
              year = year,
              `.id`= 1, 
              PR_1 = PR_1,
              PR_2 = PR_2,
              PR_3 = PR_3,
              PR_4 = PR_4,
              PR_5 = PR_5,
              Itp_1 = It_1/Nv_1*1e5, 
              Itp_2 = It_2/Nv_2*1e5, 
              Itp_3 = It_3/Nv_3*1e5, 
              Itp_4 = It_4/Nv_4*1e5, 
              Itp_5 = It_5/Nv_5*1e5, 
              H_1 = lambdaw_1/lambdas_1, 
              H_2 = lambdaw_2/lambdas_2, 
              H_3 = lambdaw_3/lambdas_3, 
              H_4 = lambdaw_4/lambdas_4, 
              H_5 = lambdaw_5/lambdas_5)
    ) %.>% 
    prep_as_plot_trajectory(., init_year = 2019) %.>% 
    mutate(., dwan = dwan_int) %>% 
    filter(., year > 2019)
  
  
# join the trajectory to the co-variate data to evaluate the vaccine efficacy and return data frame
    
comp_vacc_eff_traj 
      
      
}

#vacc_eff_grid(po_obj = po_vacc_eff_with_booster)

vacc_eff_res_path <- "../result_data/"

if(dir.exists(paste0(vacc_eff_res_path, "vacc_effect")) == FALSE) {
  dir.create(paste0(vacc_eff_res_path, "vacc_effect"))
   message("Directory 'vacc_effect' has been created, proceeding to populate!")


  vaccine_eff <- (
    pbmclapply(1:(nrow(dwan_grid)), 
               vacc_eff_grid, po_obj = po_vacc_eff_with_booster, mc.cores = detectCores()) %.>% 
      form_dfr(.) 
  )
  

  save(vaccine_eff, file = paste0(vacc_eff_res_path, "vacc_effect/vaccine_eff.rds"))

} else {

  message("Directory 'vacc_effect' already exists, moving on!")
}

load(paste0(vacc_eff_res_path, "vacc_effect/vaccine_eff.rds"))

# prepare age-specific trajectroies grouped by age-classes
# for [0,5)
vaccine_eff_1 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[1] & year < 2025 & comp_exp %in% c("PR", "H")
           ) %.>% 
    mutate(., year = year - 2020)
  )

# for [5,15)
vaccine_eff_2 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[2] & year > 2024 & year < 2035 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025)
)  
  
# for [15,25)
vaccine_eff_3 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[3] & year > 2034 & year < 2045 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025)
)  

# for [25,40)
vaccine_eff_4 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[4] & year > 2044 & year < 2060 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025)
)  

# for >40
vaccine_eff_5 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[5] & year > 2059 & year < 2101 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025)
)  


# combine the age data
vaccine_eff_plot_df <- (
  vaccine_eff_1 %.>% 
    bind_rows(., vaccine_eff_2) %.>% 
    bind_rows(., vaccine_eff_3) %.>% 
    bind_rows(., vaccine_eff_4) %.>% 
    bind_rows(., vaccine_eff_5) 
  ) %.>% 
  mutate(., pdwan = 1-exp(-18/dwan))


  
# subset from the vaccine eff data for background annotation
rect_anno_data <- (
  vaccine_eff_plot_df %.>%
    mutate(., neonatal_dose = ifelse(age_cohort == age_names[1], "si", "nein!")) %.>% 
    select(., age_cohort, neonatal_dose) %.>% 
    distinct(.)
    
  ) 


# subset from the vaccine eff for line annotation
line_anno_data <- (
  vaccine_eff_plot_df %.>% 
    filter(., comp_exp == "PR") %.>% 
    filter(., pdwan == 1-exp(-18/111)) %.>% 
    select(., year, count, age_cohort)
)


# make the PR plot with annotations
vaccine_eff_PR_plot <- (
  vaccine_eff_plot_df %.>% 
    filter(., comp_exp == "PR") %.>% 
    ggplot(.) +
    geom_rect(data = rect_anno_data,
              aes(fill = neonatal_dose), alpha = 0.3, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_line(aes(x = year, y = count, colour = pdwan, group = pdwan), size = 0.8) +
    geom_line(data = line_anno_data, aes(x = year, y = count), 
              colour = "#642B73", size = 0.8) +
    geom_hline(yintercept = 1, colour = "grey30", linetype = "dotted", size = 0.8) +
    labs(x = "Years Since Last Dose", 
         y = "Relative Prevalence", 
         colour = "P(Immune Loss By Age 18)", 
         fill = "Dose Type") +
    scale_y_continuous(trans = "log10", 
                       breaks = c(1e-1, 1, 1e1, 1e2), 
                       limits = c(1e-1, 1.3e2), 
                       labels = label_scientific()
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    scale_fill_manual(breaks = c("si", "nein!"), values = c("grey50", "black"), 
                      labels = c("Neonatal", "Booster")) +
    annotation_logticks(sides = "l") +  
    facet_grid(cols = vars(age_cohort), scales  = "free_x") +
    continuous_scale(
      "color", "modified_palette",
      mod_palette(target = pdwan_18*c(0.95, 1.05), range = vaccine_eff_plot_df$pdwan %.>% range(.)),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1),
      guide = guide_colourbar(nbin = 100, 
                              barheight = unit(7, "pt"),
                              barwidth = 9.5,
                              frame.colour = "black", 
                              ticks.colour = "black", 
                              title.position = "top", 
                              direction = "horizontal") 
    ) +
    project_theme +
    cap_axes() +
    theme(axis.text.x=element_text(angle=90, vjust = 0.5),
          legend.key.size = unit(7, 'pt')  
          ) +
    guides(fill = guide_legend(nrow = 1, 
                               title.position = "top", 
                               override.aes = list(alpha = 1))
           )
  
)

vaccine_eff_data_for_plot <- (
  vaccine_eff %.>% 
  filter(., year == max(year) & comp_exp == "Itp") %.>%
  mutate(., pdwan = 1-exp(-18/dwan))
  ) 


anno_point_vacc_eff_itp_data <- (
  vaccine_eff_data_for_plot %.>% 
    filter(., pdwan == pdwan_18)
)

# make incidence plot as a response to 
vaccine_eff_Itp_plot <- (
    vaccine_eff_data_for_plot %.>% 
    ggplot(.) +
    geom_area(aes(x = pdwan, y = count, fill = age_cohort), 
              position="stack", stat="identity") +
    geom_line(aes(x = pdwan, y = dwan*0.60), size = 1, colour = "grey30", 
              linetype = "dotdash") +
    geom_segment(data = anno_point_vacc_eff_itp_data, 
                 aes(x = pdwan, xend = pdwan, y = 0, yend = 33), 
                 size = 0.4, colour = "#642B73", lineend = "round") +
    geom_point(data = anno_point_vacc_eff_itp_data, 
               aes(x = pdwan, y = dwan*0.60), 
               size = 2, pch = 21, fill = "white", colour = "#642B73") +
    labs(y = expression(Equilibrium~Prevalence~per~10^5), 
         fill = "Age Cohort", 
         x = "P(Immune Loss By Age 18)") +
    scale_y_continuous(breaks = c(0, 40, 80, 120), limits = c(0, 130), 
                       sec.axis = sec_axis(trans =  ~./0.60, expression(Immune~Duration~(Years)))) +
    scale_x_reverse(breaks = seq(0.1,1, by = 0.3), limits = c(1, 0.06)) +
    scale_fill_brewer(palette = "Oranges", direction = -1) +
    project_theme +
    cap_axes(right = "both") +
    guides(fill = guide_legend(nrow = 2, title.position = "top")) +
    theme(legend.key.size = unit(7, 'pt'))  
  )
  


# calculate the impact over various waning values 
cal_var_impact <- function(count = 30, dwan_data = dwan_grid) {
  # browser()
  
  dwan_val <- dwan_data[count,] %.>% unlist(.) %>% unname(.)
  vec_int <- best_model_p_vec
  vec_int["dwan"] <- dwan_val
  
  
  R0 <- c(vec_int, params_for_R0) %.>% calculate_R0_mq(.)$reprodutive_number  
  
  Rp <- c(vec_int, params_for_Rp) %.>% calculate_R0_mq(.)$reprodutive_number
  
  impact <- 1-Rp/R0
  
  tibble(impact = impact, 
         Rp = Rp,
         dwan = dwan_val, 
         pdwan = 1-exp(-18/dwan))
  
}



var_vacc_imp <- map_dfr(1:nrow(dwan_grid), cal_var_impact)

anno_point_var_vacc_data <- (
  var_vacc_imp %.>% 
    filter(., pdwan == pdwan_18)
)
  

vaccine_imp_plot <- (
  var_vacc_imp %.>% 
    ggplot(., aes(x = pdwan)) +
    geom_line(aes(y = impact, linetype = "impact"), size = 1.0, colour = "grey30") +
    geom_line(aes(y = Rp/20, linetype = "Rp"), size = 1.0, colour = "grey30") +
    geom_point(data = anno_point_var_vacc_data, aes(x = pdwan, y = impact), 
               pch = 21, fill = "white", colour = "#642B73", size = 2) +
    geom_point(data = anno_point_var_vacc_data, aes(x = pdwan, y = Rp/20), 
               pch = 21, fill = "white", colour = "#642B73", size = 2) +
    labs(y = expression(Vaccine~Impact~(xi)~phantom(10)), 
         x = "P(Immune Loss By Age 18)") +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*20, expression(Reproductive~Number~(R[p]))),
                       breaks = c(0, 0.25, 0.5, 0.75),
                       limits = c(0,0.75),
                       labels = scales::percent) +
    scale_x_reverse(breaks = seq(0.1,1, by = 0.3), limits = c(1, 0.06)) +
    scale_linetype_manual(breaks = c("impact", "Rp"), 
                          values = c(1, 2), 
                          name = "Quantity", 
                          labels = c(expression(xi), expression(R[p]))) +
    project_theme +
    cap_axes(right = "both") +
    guides(linetype = guide_legend(title.position = "top", 
                                   nrow = 1, position = "none", 
                                   override.aes = list(size = 0.5)))
  )



vaccine_eff_panel_plt <- (
    vaccine_eff_PR_plot + 
    vaccine_imp_plot +      
    guide_area() +  
    vaccine_eff_Itp_plot +
    plot_layout(
      width = c(1, 0.6),
      heights = c(1, 0.8),
      guides = "collect", 
      design = "
      AC
      BD
      " 
    ) +
    plot_annotation(tag_levels = "A")&
      theme(plot.tag = element_text(size = grid_lab_size, face = "bold"))
  )





