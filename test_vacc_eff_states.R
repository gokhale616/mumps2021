# requires 
source("./00/src.R", chdir = TRUE)
source("./plotting/prep_estm_tables.R", chdir = TRUE)

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
                           extra_end_t = 2100)
)

test_vacc_traj <- (
  po_vacc_eff_with_booster %.>%  
    trajectory(., 
               param = best_model_p_vec_vacc_eff, 
               method = "ode23", format = "d") %.>% 
    select(.,  
           `.id`, 
           year, 
           starts_with("Nv"), 
           starts_with("It"), 
           #starts_with("Iw"), 
           #starts_with("Is"), 
           starts_with("lambda"),
           starts_with("PR")
           ) %.>% 
    mutate(., 
           Itp_1 = It_1/Nv_1, 
           Itp_2 = It_2/Nv_2, 
           Itp_3 = It_3/Nv_3, 
           Itp_4 = It_4/Nv_4, 
           Itp_5 = It_5/Nv_5, 
           H_1 = lambdaw_1/lambdas_1, 
           H_2 = lambdaw_2/lambdas_2, 
           H_3 = lambdaw_3/lambdas_3, 
           H_4 = lambdaw_4/lambdas_4, 
           H_5 = lambdaw_5/lambdas_5) %.>% 
    prep_as_plot_trajectory(., init_year = 1967-1/1)
  )  
  

test_vacc_traj %.>% 
  mutate(., count = ifelse(count <= 0, NA, count)) %.>% 
  filter(., year > 2019) %.>% 
  ggplot(., aes(x = year, y = count)) +
  geom_line(aes(colour = age_cohort)) +
  #geom_bar(aes(fill = age_cohort), stat = "identity") +
  facet_grid(rows = vars(comp_exp), scale = "free") +
  scale_y_continuous(trans = "log10")


# mod_mumps_covariates_sigmoidal %.>%
#   select(., year, starts_with("N_")) %.>% 
#   mutate(., `.id` = 1) %.>% 
#   prep_as_plot_trajectory(., init_year = 1967-1/1) %.>% 
#   ggplot(., aes(x = year, y = count)) +
#   geom_line(aes(colour = age_cohort)) +
#   facet_grid(rows = vars(comp_exp), scale = "free") #+
#   #scale_y_continuous(trans = "log10")




  

