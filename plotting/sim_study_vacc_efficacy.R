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
  temp_scale = 1/365.25)
  )

test_vacc_traj <- (
  po_vacc_eff_with_booster  %.>%  
    trajectory(., 
               param = best_model_p_vec_vacc_eff, 
               method = "ode23", format = "d") 
)




dwan_grid <- tibble(dwan = c(seq(1, 200, by = 10), best_model_p_vec_vacc_eff["dwan"])) %.>% arrange(.)


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
    right_join(., 
               mod_mumps_covariates_sigmoidal %.>% 
                 filter(., year > 1966) %.>% 
                 select(., year, starts_with("N")), 
               by = "year"
               ) %.>% 
    mutate(., 
           N = N_1 + N_2 + N_3 + N_4 + N_5, 
           Is = Is_1 + Is_2 + Is_3 + Is_4 + Is_5, 
           Iw = Iw_1 + Iw_2 + Iw_3 + Iw_4 + Iw_5,
           It = It_1 + It_2 + It_3 + It_4 + It_5)
    )
  
  
  # join the trajectory to the co-variate data to evaluate the vaccine efficacy and return data frame
    
comp_vacc_eff_traj %.>% 
      transmute(., 
                dwan = dwan_int,    
                year = year,#seq(0, n()-1, by = 1),
                Itp_1 = It_1/N_1, 
                Itp_2 = It_2/N_2, 
                Itp_3 = It_3/N_3,
                Itp_4 = It_4/N_4,
                Itp_5 = It_5/N_5, 
                Itp   = It/N,  
                Ef_I_1 = (Iw_1)/(Is_1), 
                Ef_I_2 = (Iw_2)/(Is_2),
                Ef_I_3 = (Iw_3)/(Is_3),
                Ef_I_4 = (Iw_4)/(Is_4),
                Ef_I_5 = (Iw_5)/(Is_5),
                Ef_I   = (Iw)/(Is)
      )
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

# interpolate to generate smoother trajectories 

# xnew2 <- seq(min(vaccine_eff$year), max(vaccine_eff$year), by = 1/52)
# 
# interpolated_vaccine_eff <- (
#   vaccine_eff %.>% 
#     interp.dataset(y = ., x=.$year, 
#                    xout = xnew2, 
#                    method = "linear") %.>%
#     as_tibble(.)) 



log_a_tic <- annotation_logticks(sides = "l")  
log_a_tic$data <- tibble(x = NA, 
                         age_cohort = age_names[1])

vacc_rel_prev <- (
  vaccine_eff %.>% 
    #filter(., year > 21) %.>% 
    select(., dwan, year, sprintf("Ef_I_%s", 1:3)) %.>% 
    gather(., key = "age_cohort", value = "rel_prev", -c(dwan, year)) %.>% 
    mutate(., 
           year = year - 22,
           age_cohort = case_when(age_cohort == "Ef_I_1" ~ age_names[1], 
                                     age_cohort == "Ef_I_2" ~ age_names[2], 
                                     age_cohort == "Ef_I_3" ~ age_names[3]))
    
)


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
  






