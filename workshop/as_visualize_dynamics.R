# source prerequisites for this project
source("../00/src.R", chdir = TRUE)
source("../fit/prep_inc_data.R", chdir = TRUE)

# source all the treated covariates
source("../fit/treat_vacc_covar.R", chdir = TRUE)
source("../simulation_study/load_pram_est.R", chdir = TRUE)


# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

leaky_model <- (
  all_result_df %.>%
    filter(., hypothesis == "leaky2", vacc_covariate == "constant") %.>%
    filter(., d_AIC == min(d_AIC))
  )
  


# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
    mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
    unlist(.) %.>% sim_p_vals(.)
)


leaky_model_p_vec <- (
  leaky_model %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>%   
    mutate(., dwan = Inf) %.>% 
    unlist(.) %.>% sim_p_vals(.)
  )


leaky_pvec <- mle_leaky2_3_constant_is$DEobj$optim$bestmem
leaky_eval <- mle_leaky2_3_constant_is$DEobj$optim$bestval

waning_pvec <- mle_waning_sigmoid_is$DEobj$optim$bestmem
waning_eval <- mle_waning_sigmoid_is$DEobj$optim$bestval


po_to_vis <- make_pomp(covar = mod_mumps_covariates_sigmoidal, 
                       extrapolate_simulation = TRUE, 
                       extra_start_t = 1967-1/52, extra_end_t = 2020, 
                       temp_scale = 1/52)


leaky_est_these <- names(leaky_pvec)

param_vals_est_l2_3 <- set_param(new_val = 2)


po_leaky_of <- (
  make_pomp(covar = mod_mumps_covariates_constant, 
            incidence_data = in_sample_mumps_case_reports, 
            params = param_vals_est_l2_3) %.>% 
    traj_objfun(., 
                est = leaky_est_these,
                ode_control = list(method = "ode23"), 
                fail.value = 1e10)
)
  

# this is the best fitting estimate for the leaky2_3 model with constant covariates
po_leaky_of(par = leaky_pvec)

leaky_pvec["epsilon2"]

leaky_pvec_mod <- leaky_pvec
leaky_pvec_mod["epsilon2"] <- 0.00

po_leaky_of(par = leaky_pvec_mod)


pvec_data <- tibble(p_vec = c(seq(0, 1, length.out = 20), leaky_pvec["epsilon2"] %.>% unname(.))) %.>% 
  arrange(., pvec) 


  rowwise(.)
  mutate(., of_eval = function() {
    leaky_pvec_mod <- leaky_pvec
    leaky_pvec_mod["epsilon2"] <- pvec
    
    po_leaky_of(par = leaky_pvec_mod)
    
  })


  
data <- map_dfr(1:nrow(pvec_data), function(i) {
  pvec_int <- pvec_data %.>% slice(., i) 
  
  pvec_val <- pvec_int %.>% unlist(.) %.>% unname(.)
  
  leaky_pvec_mod <- leaky_pvec
  leaky_pvec_mod["epsilon2"] <- pvec_val
  
  pvec_int %.>% 
    mutate(., 
           po_eval = po_leaky_of(par = leaky_pvec_mod)
           )
  
})
  
  





rp_vals_mod <-  leaky_pvec

# mod_these_params <- c("q", "sigma", "beta1", "dwan", "epsilon1", "t_intro", "p_intro", "iota", "epsilon2")
mod_these_params <- c(sprintf("q_age_%s", 1:5), "alpha", "dwan")

rp_vals_mod[mod_these_params] <- c(best_model_p_vec[sprintf("q_age_%s", 1:5)][1:4], 0.5, 0.00, Inf)


rp_vals_mod2 <- rp_vals_mod %.>% sim_p_vals(.)

rp_vals_mod2[c("p_intro", "epsilon2")] <- c(2, 0.1)

make_pomp(covar = mod_mumps_covariates_constant, 
          incidence_data = in_sample_mumps_case_reports, 
          params = param_vals_est_l2_3) %.>% 
  trajectory(., param = rp_vals_mod2, method = "ode23", format = "d") -> trajectory



# something fishy with the C compartment!! - Maybe a compiler issue! 

sample_plot_traj <- trajectory %.>% prep_as_plot_trajectory(., init_year = 1967-1/52)

sample_plot_traj %.>%
  filter(., year > 1978) %.>% 
  ggplot(., aes(x = year, y = count, colour = age_cohort)) +
  geom_line(size = 0.8) +
  facet_wrap(vars(comp_exp), scales = "free") +
  project_theme +
  scale_color_brewer(palette = "Blues", direction = -1) +
  cap_axes()



xnew <- seq(1950, 2018, by = 1/52)

system_covs <- (
  mod_mumps_covariates_sigmoidal %.>% 
    filter(., year > 1949) %.>%   
    select(., year, starts_with("N_")) %.>% 
    interp.dataset(y = ., x=.$year, 
                   xout = xnew, 
                   method = "linear") %.>% 
    as_tibble(.))  


traj_to_summarize <- (
  sample_traj %.>% 
    #mutate(., 
    #       Model= "Waning") %.>% 
  right_join(., 
             by = c("year"),  
             system_covs %.>% 
               replace_na(., 
                          system_covs %.>% 
                            filter(., year == 2018)) %.>% 
               mutate(., year = xnew)) 
  ) %.>% 
  mutate(., `.id` = 1) %.>% 
  drop_na(.)
  
  

append_rp_with_these <- c(N = 1e8, p = 0.86, nu = 1/80, ad = age_class_duration)

summarized_traj <- summarize_epidemiology(traj_cov_data = traj_to_summarize, 
                                          p_vals = c(rp_vals_mod, 
                                                     append_rp_with_these))




summarized_traj %.>% 
  select(., year, Reff, average_beta, mean_age_at_infection) %.>% 
  gather(., "summary", "value", -year) %.>% 
  ggplot(., aes(x = year, y = value)) +
  geom_line(size = 0.8) +
  facet_grid(rows = vars(summary), scales = "free_y") +
  project_theme


















  




