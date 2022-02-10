source("../00/src.R", chdir = TRUE)
source("../plotting/prep_estm_tables.R", chdir = TRUE)



# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
    select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
    mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
    unlist(.) %.>% sim_p_vals(.)
)


mai_po <- (
  mumps_covariates %.>% 
    mutate(., p1 = 0, p2 = 0) %.>% 
    make_pomp(extrapolate_simulation = TRUE, 
              extra_start_t = 2019, extra_end_t = 2020, covar = .) 
  )


pick_wts <- (
  trajectory(object = mai_po, param = best_model_p_vec, method = "ode23", format = "d") %.>% 
    select(., year, `.id`, starts_with("C_")) %.>% 
    slice(., n()) %.>% 
    prep_as_plot_trajectory(., init_year = 2019) %.>% 
    mutate(., wts = count) %.>% 
    select(., age_cohort, wts)
  ) 
  
age_data_with_wts <- (
  tibble(age = seq(0, 80, by = 1)) %.>% 
    mutate(., age_cohort = case_when(age < 5  ~ age_names[1], 
                                     age < 15 ~ age_names[2], 
                                     age < 25 ~ age_names[3], 
                                     age < 40 ~ age_names[4], 
                                     TRUE ~  age_names[5])
           ) %.>% 
    right_join(., pick_wts, by = "age_cohort") %>% 
    mutate(., 
           mean_age = sum(age*wts)/sum(wts))
  )
  
  
  
