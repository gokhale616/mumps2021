# requires 
source("../00/src.R", chdir = TRUE)
source("../plotting/prep_estm_tables.R", chdir = TRUE)


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
  extra_end_t = 2100)
  )

test_vacc_traj <- (
  po_vacc_eff_with_booster  %.>%  
    trajectory(., 
               param = best_model_p_vec_vacc_eff, 
               method = "ode23", format = "d") 
)




dwan_grid <- tibble(dwan = c(seq(1, 200, by = 5), best_model_p_vec_vacc_eff["dwan"])) %.>% arrange(.)


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
              Itp_1 = It_1/Nv_1, 
              Itp_2 = It_2/Nv_2, 
              Itp_3 = It_3/Nv_3, 
              Itp_4 = It_4/Nv_4, 
              Itp_5 = It_5/Nv_5, 
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

#vacc_eff_res_path <- "../result_data/"

# if(dir.exists(paste0(vacc_eff_res_path, "vacc_effect")) == FALSE) {
#   dir.create(paste0(vacc_eff_res_path, "vacc_effect"))
#   message("Directory 'vacc_effect' has been created, proceeding to populate!")
  

  vaccine_eff <- (
    pbmclapply(1:(nrow(dwan_grid)), 
               vacc_eff_grid, po_obj = po_vacc_eff_with_booster, mc.cores = detectCores()) %.>% 
      form_dfr(.) 
  )
  

#   save(vaccine_eff, file = paste0(vacc_eff_res_path, "vacc_effect/vaccine_eff.rds"))
#   
# } else {
#   
#   message("Directory 'vacc_effect' already exists, moving on!")
# }
# 
# load(paste0(vacc_eff_res_path, "vacc_effect/vaccine_eff.rds"))

# prepare age-specific trajectroies grouped by age-classes
# for [0,5)
vaccine_eff_1 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[1] & year < 2025 & comp_exp %in% c("PR", "H")
           ) %.>% 
    mutate(., year = year - 2020)
  )

vaccine_eff_2 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[2] & year > 2024 & year < 2035 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025)
)  
  

vaccine_eff_3 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[3] & year > 2034 & year < 2045 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025 + 10)
)  

vaccine_eff_4 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[4] & year > 2044 & year < 2060 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025 + 15)
)  

vaccine_eff_5 <- (
  vaccine_eff %.>% 
    filter(., 
           age_cohort == age_names[5] & year > 2059 & comp_exp %in% c("PR", "H")
    ) %.>% 
    mutate(., year = year - 2025 + 35 )
)  

  
plot_vacc_eff <- function(vacc_eff_df, 
                          ylab = "Relative Hazard",
                          breaks = c(1e-2, 1, 1e2, 1e4, 1e6), 
                          limits = c(1e-2, 1e6), 
                          xbr = seq(0,4,1), 
                          legpos = "none") {
  
  # log_a_tic <- annotation_logticks(sides = "l")  
  # log_a_tic$data <- tibble(x = NA, 
  #                          age_cohort = age_names[1])
  
  # browser()
  
  vacc_eff_df %.>% 
    ggplot(., aes(x = year, y = count)) +
    geom_line(size = 0.8, aes(colour = dwan, group = dwan)) +
    geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed", size = 1) +
    labs(x = "Years Since Immunizaion", 
         y = ylab, 
         colour = "Immune\nDuration") +
    scale_y_continuous(trans = "log10", 
                       breaks = breaks, 
                       limits = limits
    ) +
    scale_x_continuous(breaks = xbr) +
    annotation_logticks(sides = "l") +  
    continuous_scale(
      "color", "modified_palette",
      modify_palette_by_target(plot_var = unique(.$dwan), 
                               target = best_model_p_vec_vacc_eff["dwan"]),
      guide = guide_colourbar(nbin = 100, 
                              frame.colour = "black", 
                              ticks.colour = "black", 
                              title.position = "top", 
                              direction = "horizontal") 
    ) +
    project_theme +
    cap_axes() +
    theme(legend.position = legpos)
  
}



vaccine_eff_1_plot <- (
    vaccine_eff_1 %.>% 
    filter(., 
           comp_exp == "H") %.>% 
    plot_vacc_eff(., 
                  breaks = c(1e-2, 1e-1, 1, 1e1, 1e2), 
                  limits = c(1e-2, 1e2)
                  )
  )


vaccine_eff_2_plot <- (
  vaccine_eff_2 %.>% 
    filter(., 
           comp_exp == "H") %.>% 
    plot_vacc_eff(., 
                  breaks = c(1e-2, 1e-1, 1, 1e1, 1e2), 
                  limits = c(1e-2, 1e2), 
                  xbr = seq(0, 3, by = 2)
    )
)


vaccine_eff_3_plot <- (
  vaccine_eff_3 %.>% 
    filter(., 
           comp_exp == "H") %.>% 
    plot_vacc_eff(., 
                  breaks = c(1e-2, 1e-1, 1, 1e1, 1e2), 
                  limits = c(1e-2, 1e2), 
                  xbr = seq(20, 29, by = 2)
    )
)





vaccine_eff %.>% 
  ggplot(., aes(x = year, y = count)) +
  geom_line(size = 0.8, aes(colour = dwan, group = dwan)) +
  # geom_line(data = vaccine_rel_prev_sub,
  #           colour = "#642B73", 
  #           size = 0.8) +
  # geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed", size = 1) +
  labs(x = "", 
       y = "Relative Infection Prevalence", 
       colour = "Immune\nDuration") +
  facet_grid(cols = vars(factor(age_cohort, levels = age_names[1:5])), 
             rows = vars(comp_exp), scales = "free") +
  scale_y_continuous(trans = "log10"#, 
                     #breaks = c(1e-2, 1, 1e2, 1e4, 1e6), 
                     #limits = c(1e-2, 1e6)
                     ) +
  log_a_tic +
  continuous_scale(
    "color", "modified_palette",
    modify_palette_by_target(plot_var = unique(.$dwan), 
                             target = best_model_p_vec_vacc_eff["dwan"]),
    guide = guide_colourbar(nbin = 100, 
                            frame.colour = "black", 
                            ticks.colour = "black", 
                            title.position = "top", 
                            direction = "horizontal") 
  ) +
  project_theme +
  cap_axes()
  
  

if(FALSE) {

vacc_prev <- (
  vaccine_eff %.>% 
    #filter(., year > 21) %.>% 
    select(., dwan, year, sprintf("Itp_%s", 1:3)) %.>% 
    gather(., key = "age_cohort", value = "prev", -c(dwan, year)) %.>% 
    mutate(., 
           year = year - 22, 
           age_cohort = case_when(age_cohort == "Itp_1" ~ age_names[1], 
                                     age_cohort == "Itp_2" ~ age_names[2], 
                                     age_cohort == "Itp_3" ~ age_names[3]))
  
)

# generating subset to plot over the plot
vaccine_rel_prev_sub <- (
  vacc_rel_prev %.>%
  filter(., 
         age_cohort %in% age_names[1:3] & dwan == 111)
  )


vaccine_prev_sub <- (
  vacc_prev %.>%
    filter(., 
           age_cohort %in% age_names[1:3] & dwan == 111)
)




vacc_rel_prev_plot <- (
  vacc_rel_prev %.>% 
    ggplot(., aes(x = year, y = rel_prev)) +
    geom_line(size = 0.8, aes(colour = dwan, group = dwan)) +
    geom_line(data = vaccine_rel_prev_sub,
                colour = "#642B73", 
              size = 0.8) +
    geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed", size = 1) +
    labs(x = "", 
         y = "Relative Infection Prevalence", 
         colour = "Immune\nDuration") +
    facet_grid(cols = vars(factor(age_cohort, levels = age_names[1:3])), scales = "fixed") +
    scale_y_continuous(trans = "log10", 
                       breaks = c(1e-2, 1, 1e2, 1e4, 1e6), 
                       limits = c(1e-2, 1e6)) +
    log_a_tic +
    continuous_scale(
      "color", "modified_palette",
      modify_palette_by_target(plot_var = unique(.$dwan), 
                               target = best_model_p_vec_vacc_eff["dwan"]),
      guide = guide_colourbar(nbin = 100, 
                              frame.colour = "black", 
                              ticks.colour = "black", 
                              title.position = "top", 
                              direction = "horizontal") 
    ) +
    project_theme +
    cap_axes()
  ) 



vacc_prev_plot <- (
  vacc_prev %.>% 
    ggplot(., aes(x = year, y = prev)) +
    geom_line(size = 0.8, aes(colour = dwan, group = dwan)) +
    geom_line(data = vaccine_prev_sub,
              colour = "#642B73", 
              size = 0.8) +
    labs(x = "Years Since Last Vaccination", 
         y = "Infection Prevalence", 
         colour = "Immune\nDuration") +
    facet_grid(cols = vars(factor(age_cohort, levels = age_names[1:3])), scales = "fixed") +
    scale_y_continuous(trans = "log10", breaks = c(1e-3, 1e-2, 1e-1, 1e0, 1e1), limits = c(1e-3, 1e1)) +
    log_a_tic +
    continuous_scale(
      "color", "modified_palette",
      modify_palette_by_target(plot_var = unique(.$dwan), 
                               target = best_model_p_vec_vacc_eff["dwan"]),
      guide = guide_colourbar(nbin = 100, 
                              frame.colour = "black", 
                              ticks.colour = "black", 
                              title.position = "top", 
                              direction = "horizontal") 
    ) +
    project_theme +
    cap_axes() +
    theme(legend.position = "none")
) 




vacc_RR_plt <- (
  vacc_rel_prev_plot  + vacc_prev_plot  +
  plot_layout(
    guides = "collect", 
    design = "
    A
    B
    " 
  ) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
  )
  

}




