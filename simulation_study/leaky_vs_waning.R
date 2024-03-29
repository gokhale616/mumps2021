library(furrr)
library(rioja)
library(progressr)

# might have to change this to src_cl if cluster is needed
source("../00/src_cl.R", chdir = TRUE)
source("./load_pram_est.R", chdir = TRUE)
source("../fit/treat_vacc_covar.R", chdir = TRUE)

# simulate dynamics to the check simulation match to theoretical results 
fin_year <- 2020
fine_res <- 1/12

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
                       temp_scale = fine_res)


leaky_waning_effect <- function(x = 1, param_grid, 
                                po_obj = po_to_sim,
                                def_param = best_model_p_vec, 
                                p) {
  
  
  # progressor
  p()
  #browser()
  # select value of leakiness and waning to be substituted
  val_to_replace <- (
    param_grid %.>% 
      slice(., x) %.>% 
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
  
  
  if(fine_res < 1) {
    
    interp_range <- c(1965, fin_year)
    
    xnew <- seq(interp_range[1], interp_range[2], by = fine_res)
    
    interpolated_covs <- (
      mod_mumps_covariates_sigmoidal %.>% 
        bind_rows(., 
                  mod_mumps_covariates_sigmoidal %.>% 
                    slice(., rep(n(), interp_range[2]-2018))) %.>% 
        mutate(., 
               year = seq(1910, interp_range[2], by = 1)) %.>% 
        filter(., year > 1965) %.>%   
        select(., year, starts_with("N_")) %.>% 
        interp.dataset(y = ., x=.$year, 
                       xout = xnew, 
                       method = "linear") %.>% 
        as_tibble(.)
    )  
    
    
  } else {
    interpolated_covs <- mod_mumps_covariates_sigmoidal
  }
  
  
  # combine with covariate data
  comp_traj_w_cov <- (
    comp_traj %.>% 
      right_join(., 
                 interpolated_covs %.>% 
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
      filter(., year >= 2002 & Reff < 1) %.>%
      arrange(., year) %.>%   
      slice(., 1) %.>% 
      replace(., is.na(.), 0)
    )  
  
  #print((x/nrow(param_grid))*100)
  
  # collect outcome
  outcome <- (
    reemergence_stats %.>%
      mutate(., 
             epsilon = param_grid$epsilon2[x],
             dwan = param_grid$dwan[x],
             time_of_reemergence = ifelse(year == 0, 2019, year)) %.>% 
      select(., -year))
  
  outcome
  
}


test_param_grid <- expand.grid(epsilon2 = seq(0, 1, length.out = 10), 
                               dwan = seq(5, 300, length.out = 10))


if(FALSE) {
tic()
map(1:nrow(test_param_grid), leaky_waning_effect, param_grid = test_param_grid)
toc()
}



tic()
plan("multisession", workers = detectCores())

with_progress({
  
  p <- progressor(steps = nrow(test_param_grid))
  
  sim_study_res <- future_map(1:nrow(test_param_grid), 
                              leaky_waning_effect, 
                              param_grid = test_param_grid, 
                              p = p)
  })

toc()

#save(sim_study_res, file = "../result_data/sim_study/sim_study_res.rds")
