
# after completion add this to the './00/src.R'
source("../00/cov_make_pomp_model_vacc_effectiveness.R", chdir = TRUE)




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
  temp_scale = 1)
  )
# without booster
po_vacc_eff_without_booster <- (
  make_pomp_vacc_effective(., 
                           covar = mod_mumps_covariates_sigmoidal %.>% mutate(., p2 = 0), 
                           temp_scale = 1)
)

dwan_grid <- tibble(dwan = c(seq(1, 200, by = 10), best_model_p_vec_vacc_eff["dwan"])) %.>% arrange(.)


vacc_eff_grid <- function(x = 112, dwan_data = dwan_grid, 
                          po_obj
                          ) {
  
  # collect the control parameter value
  dwan_int <- dwan_grid %.>%  slice(., x) %.>% unlist(.) %.>% unname(.)
  
  # replace in the defaults
  p_for_sim <- best_model_p_vec_vacc_eff
  p_for_sim["dwan"] <- dwan_int
  
  # Generate trajectories from all the compartments 
  comp_vacc_eff_traj <- (
    trajectory(po_obj, 
               param = p_for_sim, 
               method = "ode23", format = "d")) %.>% 
    mutate(., 
           Is = Is_1 + Is_2 + Is_3 + Is_4 + Is_5, 
           Iw = Iw_1 + Iw_2 + Iw_3 + Iw_4 + Iw_5,
           It = It_1 + It_2 + It_3 + It_4 + It_5)
  
  
  # join the trajectory to the co-variate data to evaluate the vaccine efficacy and return data frame
    
comp_vacc_eff_traj %.>% 
      select(., 
             year, 
             starts_with("Is"), 
             starts_with("Iw"), 
             starts_with("It")
      ) %.>%
      transmute(., 
                dwan = dwan_int,    
                year = seq(0, n()-1, by = 1),
                Ef_I_1 = (Iw_1)/(Is_1), 
                Ef_I_2 = (Iw_2)/(Is_2),
                Ef_I_3 = (Iw_3)/(Is_3),
                Ef_I_4 = (Iw_4)/(Is_4),
                Ef_I_5 = (Iw_5)/(Is_5),
                Ef_I   = (Iw)/(Is)
      )
}

tic()
vaccine_eff_without_booster <- (
  pbmclapply(1:(nrow(dwan_grid)), 
             vacc_eff_grid, po_obj = po_vacc_eff_without_booster, mc.cores = 8) %.>% 
    form_dfr(.) %.>% 
    mutate(., vacc_schedule = "Without\nBooster")
    
  )
toc()

tic()
vaccine_eff_with_booster <- (
  pbmclapply(1:(nrow(dwan_grid)), 
             vacc_eff_grid, po_obj = po_vacc_eff_with_booster, mc.cores = 8) %.>% 
    form_dfr(.) %.>%
    mutate(., vacc_schedule = "With\nBooster")
)
toc()



vaccine_eff <- (
  vaccine_eff_without_booster %.>% 
    bind_rows(., vaccine_eff_with_booster) %.>% 
    gather(., key = "age_cohort", value = "hazard", -c(year, dwan, vacc_schedule)) %.>%
    mutate(., log_hazard = log(hazard, base = 10))
  )


ncolours <- seq(0, 1, length.out = 100)

modify_palette_by_target <- function(colours = seq_gradient_pal("#dd3e54", "#6be585")(ncolours), 
                                     target, 
                                     plot_var, 
                                     replace_colour = "#642B73") {
  
  # browser()
  # rank data to lineraize values
  if(target %in% plot_var) {
    message("target exists!!")
    int_target = target
  } else {
    # browser()
    message("value closest to the target used")
    int_target = plot_var[which.min(abs(target - plot_var))]
  }
  
  
  ranked_plot_var <- rank(plot_var)
  
  ranked_target <- ranked_plot_var[which(plot_var == int_target)]
  
  # normalize the target interval 
  ranked_range <- range(ranked_plot_var) 
  
  target_interval <- ranked_plot_var[which(plot_var == int_target)]*c(0.95, 1.05)
  
  normalized_target_interval <- (target_interval - ranked_range[1])/diff(ranked_range)
  
  # browser()
  
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    # Decide what values to replace
    # browser()
    ranked_x <- rank(x)
    
    normalized_ranked_x <- (ranked_x - range(ranked_x)[1])/diff(range(ranked_x))
    
    replace <- normalized_ranked_x > normalized_target_interval[1] & 
      normalized_ranked_x < normalized_target_interval[2]
    
    out <- ramp(x)
    # Actually replace values
    out[replace] <- replace_colour
    out
  }
}



vaccine_eff_sub <- (
  vaccine_eff %.>%
  filter(., 
         age_cohort == sprintf("Ef_I_%s", 1:2) & dwan == 111) %.>% 
  mutate(., age_cohort = ifelse(age_cohort == "Ef_I_1", age_names[1], age_names[2]))
)

vaccine_eff %.>% 
  filter(., age_cohort == sprintf("Ef_I_%s", 1:2)) %.>% 
  mutate(., age_cohort = ifelse(age_cohort == "Ef_I_1", age_names[1], age_names[2])) %.>% 
  ggplot(., aes(x = year, y = log_hazard)) +
  geom_line(size = 0.8, aes(colour = dwan, group = dwan)) +
  geom_line(data = vaccine_eff_sub,
              colour = "#642B73", 
            size = 0.8) +
  geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed", size = 1) +
  labs(x = "Time Since Vaccination (Years)", 
       y = expression(log[10](Infection~RR~Post-Vaccination)), 
       colour = "Immune\nDuration") +
  facet_grid(rows = vars(age_cohort), cols = vars(vacc_schedule), scales = "free_y") +
  # scale_color_gradient(low = "#dd3e54", high = "#6be585", breaks = ) +
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
  












