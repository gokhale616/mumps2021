# might have to change this to src_cl if cluster is needed
source("../00/src.R", chdir = TRUE)
source("../plotting/prep_estm_tables.R", chdir = TRUE)
source("../fit/treat_vacc_covar.R", chdir = TRUE)

# simulate dynamics to the check simulation match to theoretical results 
fin_year <- 2018

# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
    mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
    unlist(.) %.>% sim_p_vals(.)
)


# define a pomp obect to summarize the epidemiology
po_to_sim <- make_pomp(., 
                       covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1965, extra_end_t = fin_year, 
                       temp_scale = 1/1)


leaky_waning_effect <- function(c = 1, param_grid, 
                                po_obj = po_to_sim,
                                def_param = best_model_p_vec) {
  
  # browser()
  # select value of leakiness and waning to be substituted
  val_to_replace <- (
    param_grid %.>% 
      slice(., c) %.>% 
      unlist(.)  
    )
  
  
  # substitute defaults   
  def_param[c("dwan", "epsilon2", "t_intro", "p_intro")] <- c(val_to_replace[c("dwan", "epsilon2")], 2002, 3)
  
  # Generate trajectories from all the compartments 
  comp_traj <- (
    trajectory(po_to_sim, 
               param = def_param, 
               method = "ode23", format = "d") 
  )
  
  
  # combine with covariate data
  comp_traj_w_cov <- (
    comp_traj %.>% 
      right_join(., 
                 mod_mumps_covariates_sigmoidal %.>% 
                   filter(., year > 1965),
                 by = "year") %.>% 
      as_tibble(.) %.>% 
      mutate(., `.id` = as.numeric(`.id`)) %.>% 
      drop_na(.)
  )
  

  # summarize the trajectory 
  summarized_traj <- summarize_epidemiology(traj_cov_data = comp_traj_w_cov,
                                            p_vals = c(def_param,
                                                       params_for_Rp)) %.>% 
    select(., year, starts_with("I2p"), starts_with("Cs"), Reff)
  
  
  # collect stats
  reemergence_stats <- (
    summarized_traj %.>% 
      filter(., year > 2002 & Reff < 1) %.>%
      arrange(., year) %.>%   
      slice(., 1) %.>% 
      replace(., is.na(.), 0)
    )  
  
  # collect outcome
  outcome <- (
    reemergence_stats %.>%
      mutate(., 
             epsilon = param_grid$epsilon2[c],
             dwan = param_grid$dwan[c],
             time_of_reemergence = ifelse(year == 0, 2019, year)) %.>% 
      select(., -year))
  
  outcome
  
}


test_param_grid <- expand.grid(epsilon2 = seq(0, 1, by = 0.1), 
                               dwan = seq(5, 300, length.out = 10))

sim_study_outcome <- map_dfr(1:100, leaky_waning_effect, param_grid = test_param_grid)
