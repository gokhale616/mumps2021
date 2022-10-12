# load all the prereq functions - make sure to change to src_cl.R when using it in cluster
source("../../00/src_cl.R", chdir = TRUE) 

# load the co-variate treatment script to make available for further functions  
source("../treat_vacc_covar.R", chdir = TRUE)

# load the in-fit incidence data to feed the times to the simulator
source("../prep_inc_data.R", chdir = TRUE)

# load the param constraints to extract the parameter namaes
source("../param_constraints.R", chdir = TRUE)

# some minor settings for the optimizer hyper parameter
bs_n_cores <- detectCores()
bs_nsim <- 5000
optim_maxit <- 1e5
optim_reltol <- 1e-5

# set seed for reproducibility
set.seed(986747881L)
RNGkind("L'Ecuyer-CMRG")

# load all of the MLE result objects in the global environment 
mle_result_path <- "../../result_data/mle"

list.files(path = mle_result_path, 
           full.names = TRUE)[-29] %.>% 
  lapply(., load,  envir = .GlobalEnv)


# One of the waning hypothesis was wrongly named to noloss
# correct the hypothesis name 
mle_waning_constant_is$Hypothesis <- "waning_constant"

# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_gamma_n_2_waning_slow_is, mle_gamma_n_2_waning_sigmoid_is, 
    mle_gamma_n_2_waning_rapid_is, mle_gamma_n_2_waning_constant_is,
    mle_gamma_waning_slow_is, mle_gamma_waning_sigmoid_is, 
    mle_gamma_waning_rapid_is, mle_gamma_waning_constant_is,
    mle_waning_slow_is, mle_waning_sigmoid_is, mle_waning_rapid_is, mle_waning_constant_is,
    mle_leaky2_2_slow_is, mle_leaky2_2_sigmoid_is, mle_leaky2_2_rapid_is, mle_leaky2_2_constant_is,
    mle_leaky2_3_slow_is, mle_leaky2_3_sigmoid_is, mle_leaky2_3_rapid_is, mle_leaky2_3_constant_is,
    mle_leaky2_4_slow_is, mle_leaky2_4_sigmoid_is, mle_leaky2_4_rapid_is, mle_leaky2_4_constant_is, 
    mle_slow_is, mle_sigmoid_is, mle_rapid_is, mle_constant_is
  )
)

# hard coding some parameter values for use downstream within estimation of R0s and Rps
# calculate the reproductive numbers 

mk_result_df <- function(c = 1, res = result_list) {
  # to look into the function at a specific iteration
  # if (c == 20) {browser()}
  
  # collect qualitative covariates
  extra_params <- res[[c]]$Hypothesis %.>% str_split(., pattern = "_") %.>% unlist(.)
  
  if(length(extra_params) == 3) {
    hypo_covar <- extra_params[1]
    p_intro    <- extra_params[2] %.>% as.numeric(.)
    vacc       <- extra_params[3]
  } else {
    hypo_covar <- extra_params[1]
    p_intro    <- NA
    vacc       <- extra_params[2]
  }
  
  # collate results in a dataframe
  res[[c]]$DEobj$optim$bestmem %.>% 
    as.list(.) %.>% 
    as_tibble(.) %>% 
    mutate(., 
           loglik = -res[[c]]$DEobj$optim$bestval, 
           npar = res[[c]]$DEobj$optim$bestmem %.>% length(.), 
           AIC = calculate_aic(loglik, npar), 
           hypothesis = hypo_covar, 
           p_intro    = p_intro,
           vacc_covariate = vacc) %.>% 
    select(., -c(loglik, npar))
  
}

# cllect all mles in a single data frame
all_mle_df <- map_dfr(1:length(result_list), mk_result_df)
  

# reset "NA" parameters to default and collect the four best models to simulate 
best_mle_df <- (
  all_mle_df %.>% 
    group_by(., hypothesis) %.>% 
    filter(., AIC == min(AIC)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC) %.>% 
    mutate(., 
           dwan = ifelse(is.na(dwan) == TRUE, Inf, dwan), 
           epsilon2 = ifelse(is.na(epsilon2) == TRUE, 0, epsilon2), 
           p_intro = ifelse(is.na(p_intro) == TRUE, 6, p_intro), 
           t_intro = ifelse(is.na(t_intro) == TRUE, 3e3, t_intro))
  )


# produce pomp pomp objects based on covar and type wanig distribution
make_appropriate_pomp <- function(covar_name, hypo_name, 
                                  incidence_data) {
  # generate the appropriate pomp  object 
  # generate the appropriate pomp  object 
  if(covar_name == "constant") {
    if (hypo_name == "gwaning") {
      mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_constant, 
                                incidence_data = incidence_data)
    } else if (hypo_name == "gwaning2") {
      mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_constant, 
                                    incidence_data = incidence_data) 
    } else {
      mod_po <- make_pomp(covar = mod_mumps_covariates_constant, 
                          incidence_data = incidence_data)
    }
  } else if(covar_name == "rapid") {
    if (hypo_name == "gwaning") {
      mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_rapid, 
                                incidence_data = incidence_data)
    } else if (hypo_name == "gwaning2") {
      mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_rapid, 
                                    incidence_data = incidence_data)
    } else {
      mod_po <- make_pomp(covar = mod_mumps_covariates_rapid, 
                          incidence_data = incidence_data)
    }
  } else if(covar_name == "sigmoid") {
    if (hypo_name == "gwaning") {
      mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                incidence_data = incidence_data)
    } else if (hypo_name == "gwaning2") {
      mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_sigmoidal, 
                                    incidence_data = incidence_data)  
    } else {
      mod_po <- make_pomp(covar = mod_mumps_covariates_sigmoidal, 
                          incidence_data = incidence_data)
    }
  } else if(covar_name == "slow") {
    if (hypo_name == "gwaning") {
      mod_po <- make_gamma_pomp(covar = mod_mumps_covariates_slow, 
                                incidence_data = incidence_data)
    } else if (hypo_name == "gwaning2") {
      mod_po <- make_gamma_n_2_pomp(covar = mod_mumps_covariates_slow, 
                                    incidence_data = incidence_data)  
    } else {
      mod_po <- make_pomp(covar = mod_mumps_covariates_slow, 
                          incidence_data = incidence_data)
    }
  } else {
    stop("invalid covar name! choose from {'slow', 'sigmoid', 'rapid', 'constant'}")
  }
  
  mod_po
}


# prepare simulated trajectories
sim_traj_at_mle <- (
  map_dfr(1:nrow(best_mle_df), function(c, res = best_mle_df) {
    
    # browser()
    sliced_res <- res %.>% slice(., c)
    
    # collect some information for later  
    hypo_name <-  sliced_res %.>% select(., hypothesis) %.>% unlist(.)
    covar_name <- sliced_res %.>% select(., vacc_covariate) %.>% unlist(.)
    
    
    mod_po <- make_appropriate_pomp(covar_name = covar_name, hypo_name = hypo_name, 
                                    incidence_data = in_sample_mumps_case_reports)
    
    
    # produce a vector of parameters to simulate from
    sim_params <- (
      sliced_res %.>% 
        select(., -c(hypothesis, vacc_covariate)) %.>% 
        unlist(.) %.>% 
        sim_p_vals(.)
    )
    
    # simulate the observation model and add the hypothesis name 
    sim_params %.>% 
      sim_obs_model(mod_po, params = ., times = time(mod_po), nsim = bs_nsim, 
                    summarise_obs_process = FALSE, 
                    root_transform = FALSE) %.>% 
      mutate(., 
             hypothesis = hypo_name, 
             covar = covar_name)
    }) %.>% 
    # format correctly for later use - when being used in the optimizer i.e.
    spread(., key = age_class, value = cases) %.>% 
    mutate(., total = NA_real_)
  )
  
  
# sim_traj_at_mle %.>%
#   filter(., hypothesis == "gwaning2") %.>%
#   gather(., key = "age", value = "cases", -c(year, `.id`, hypothesis, covar)) %.>%
#   ggplot(., aes(x = year, y = cases, group = `.id`)) +
#   geom_line()+
#   facet_grid(rows = vars(age), scales = "free_y")
  


# setting the parameter names to estimate   
param_to_est <- (
  list(
    leaky2 = param_range_leaky2$lower %.>% names(.),
    waning = param_range_waning$lower %.>% names(.),
    gwaning = param_range_waning$lower %.>% names(.),
    gwaning2 = param_range_waning$lower %.>% names(.),
    noloss = hypo_ind_params$lower %.>% names(.))
  )
  


# transform the initial guesses 
trans_best_mle_df <- (
  best_mle_df %.>% 
    mutate_at(., .vars = vars(starts_with("q_age_"), starts_with("rho_age_"), 
                              beta1, epsilon2), .funs = logit) %.>% 
    mutate_at(., .vars = vars(starts_with("psi_"), 
                              sigma, dwan, t_intro), .funs = log)
)

# this function is an optim wrapr to reduce code repetition
# optimization method has has been set to Nelder-Mead
optim_wrpr <- function(par, fn) {
  optim(par = par, fn = fn, method = "Nelder-Mead", 
        control = list(maxit = optim_maxit, trace = 1, reltol = optim_reltol))
  
}

##############################################################################################################
################ Following are series functions to execute trajectory matching protocols #####################
##############################################################################################################

# hypothesis independent parameters to be transformed
indi_log_trans_param <- c(sprintf("psi_%s", c(1:5, "u")), "sigma")
indi_logit_trans_param <-  c(sprintf("rho_age_%s", c(1:5, "u")), sprintf("q_age_%s", c(1:5)), "beta1")


###### write a traj match protocol function for the noloss model ######
trajmatch_protocol <- function(x, hypo_name, covar_name, 
                               params = param_vals_est, 
                               log_vars, logit_vars) {
  # browser()
  # pick the i^{th} trajectory
  simulated_inc_data <- (
    sim_traj_at_mle %.>% 
      filter(., `.id` == x & hypothesis == hypo_name) %.>% 
      select(., -c(`.id`, hypothesis, covar))
  )
  
  param_to_est_int <- param_to_est[[hypo_name]]
  
  obj_fun <- (
    make_appropriate_pomp(covar_name = covar_name, 
                          hypo_name = hypo_name, 
                          incidence_data = simulated_inc_data) %.>% 
      traj_objfun(., 
                  est = param_to_est_int, 
                  paramnames = param_names_est,
                  params = param_vals_est, 
                  partrans = parameter_trans(log = c(indi_log_trans_param, log_vars), 
                                             logit = c(indi_logit_trans_param, logit_vars)), 
                  fail.value = 1e20, 
                  ode_control = list(method = "ode23")
      )
  )
  
  # set up an initial guess
  init_guess <- (
    trans_best_mle_df %.>% 
      filter(., hypothesis == hypo_name) %.>% 
      select(., -c(hypothesis, vacc_covariate)) %.>% 
      unlist(.)[param_to_est_int]
  )
  
  
  res <- optim_wrpr(par = init_guess, fn = obj_fun)
  
  res$par %.>% 
    as.list(.) %.>% 
    as_tibble(.) %.>% 
    mutate(., 
           logLik = -res$value, 
           `.id` = x, 
           hypo_name = hypo_name, 
           covar_name = covar_name)
  

}

##############################################################################################################
##############################################################################################################
##############################################################################################################

CI_res_data_path <- "../../result_data/CI/"
CI_res_path_fn <- function(fn) { paste0(CI_res_data_path, fn)}



# new seed used for leaky2
#set.seed(886747981)




