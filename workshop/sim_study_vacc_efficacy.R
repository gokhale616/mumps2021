source("../00/src.R", chdir = TRUE)

source("../00/cov_make_pomp_model_vacc_effectiveness.R", chdir = TRUE)

# load the parameter estimates from all other models 
source("../plotting/prep_estm_tables.R", chdir = TRUE)

best_model_vacc_eff <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec_vacc_eff <- (
  best_model_vacc_eff %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar, 
                 epsilon2, p_intro, t_intro, impact)) %.>% 
    unlist(.) %.>% sim_p_vals(., default_p_vals = param_vals_vacc_effective)
)


# simulate dynamics

# define a pomp obect to summarize the epidemiology
po_vacc_eff <- make_pomp_vacc_effective(., 
                                        covar = mod_mumps_covariates_sigmoidal, 
                                        temp_scale = 1/52)


#vacc_eff_frid <- function(c, dwan_data) {
  
  # Generate trajectories from all the compartments 
  comp_vacc_eff_traj <- (
    trajectory(po_vacc_eff, 
               param = best_model_p_vec_vacc_eff, 
               method = "ode23", format = "d") 
  ) %.>% 
    mutate(., 
           Is = Is_1 + Is_2 + Is_3 + Is_4 + Is_5, 
           Iw = Iw_1 + Iw_2 + Iw_3 + Iw_4 + Iw_5,
           It = It_1 + It_2 + It_3 + It_4 + It_5)
  
  
  # join the trajectory to the covariates 
  
  interpolated_covs %.>% 
    filter(., year > 1967-1/52) %.>% 
    mutate(., N = N_1 + N_2 + N_3 + N_4 + N_5) %.>% 
    right_join(.,
               comp_vacc_eff_traj, 
               join_by = "year") %.>% 
    select(., 
           year, 
           starts_with("N"), starts_with("Is"), 
           starts_with("Iw"), starts_with("It")
    ) %.>%
    transmute(., 
              year = year, 
              Ef_I_1 = 1 - (Iw_1/N_1)/(Is_1/N_1), 
              Ef_I_2 = 1 - (Iw_2/N_2)/(Is_2/N_2),
              Ef_I_3 = 1 - (Iw_3/N_3)/(Is_3/N_3),
              Ef_I_4 = 1 - (Iw_4/N_4)/(Is_4/N_4),
              Ef_I_5 = 1 - (Iw_5/N_5)/(Is_5/N_5),
              Ef_I = 1 - (Iw/N)/(Is/N)
    )
  
  
  
#}











