# this script has utility functions for generating pomp style objective functions to be used in 
# the optimization routines 
# these are mostly wrappers

# this functions generates a utility function to combine pomp with deoptim

# wrapper - this function produces pomp style objective function for DEoptim()
DE_traj_objfun <- function(x, objfun, est, 
                           seed = 986747881L) {
  
  # produce an empty vector
  par_v <- rep(NA, times = length(est))
  # name the vector
  names(par_v) <- est
  
  # assign the guess to a new internal vector
  par_v[est] <- un_list_name(x)
  
  # supply to the pomp objective function
  objfun(par = par_v)
  }

# DE_traj_objfun <- function(x, objfun, est) {
#   # produce an empty vector 
#   par_v <- rep(NA, times = length(est))
#   # name the vector 
#   names(par_v) <- est
#   
#   # assign the guess to the new internal vector
#   par_v[est] <- unlist(x) %.>% unname(.)
#   
#   # supply to the pomp objective function
#   objfun(par = par_v)
#}

#  this is a wrapper function that carries out parameter estimation using DEoptim()
DE_traj_match <- function(param_constraints,
                          waning_distribution = "exp",
                          params = p_vals,
                          ninit = np_val, ode_control = NULL, 
                          hypo_name, best_past_est = NULL,
                          seed = 986747881L,
                          other_DE_controls = my_controls,
                          ...) {
  
  message(cat(c("Est: ", names(param_constraints$lower))))
  # browser()
  # generate a pomp objective function - this step also includes defining the pomp object
  # NOTE: depending upon the argument, "waning_distribution" -- exponential or gamma distributed waning rate 
  # will be chosen.
  # NOTE: names of the parameters estimated are taken from the lower constraint vector   
  
  if(waning_distribution == "exp") {
    pomp_objfun <- (
      make_pomp(...) %.>%
        # define the objective function 
        traj_objfun(., 
                    est = names(param_constraints$lower), 
                    params = params, fail.value = 1e20, 
                    ode_control = ode_control)
    ) 
  } else if(waning_distribution == "gamma") {
    pomp_objfun <- (
      make_gamma_pomp(...) %.>%
        # define the objective function 
        traj_objfun(., 
                    est = names(param_constraints$lower), 
                    params = params, fail.value = 1e20, 
                    ode_control = ode_control)
    ) 
  } else if(waning_distribution == "gamma_n2") {
    pomp_objfun <- (
      make_gamma_n_2_pomp(...) %.>% 
        traj_objfun(.,
                    est = names(param_constraints$lower),
                    params = params, fail.value = 1e20,
                    ode_control = ode_control)
      
    )  
  } else {
      stop("Invalid waning rate distribution: use 'exp' or 'gamma'")
    }
  
  # browser()
  # generate a grid of initial guesses
  if(is.null(best_past_est)) {
    
    init_guess_grid <- sobol_design(lower = param_constraints$lower, 
                                    upper = param_constraints$upper, 
                                    nseq = ninit)   
  } else {
    
    init_guess_grid <- sobol_design(lower = param_constraints$lower, 
                                    upper = param_constraints$upper, 
                                    nseq = ninit) %.>% 
      bind_rows(., best_past_est) %.>% 
      replace(., is.na(.), 0) 
    
  }
  
  # browser()
  
  # set seed for reproducible parallel computation
  set.seed(986474881L)
  
  RNGkind("L'Ecuyer-CMRG")
  
  # set multi-core cluster for parallel computation optimal solution
  no_cores <- detectCores()  
  
  registerDoParallel(cores = no_cores)  
  
  cl <- makeCluster(no_cores, type="FORK")
  
  
  # feed all this info to the evolutionary optimizer
  DEobj <- DEoptim(fn = DE_traj_objfun, 
                   est = names(param_constraints$lower), 
                   objfun = pomp_objfun,
                   seed = seed,
                   lower = param_constraints$lower, 
                   upper = param_constraints$upper, 
                   control = c(other_DE_controls, 
                               list(cluster = cl, 
                                    NP = nrow(init_guess_grid), 
                                    initialpop = init_guess_grid %.>% 
                                      as.matrix(.))
                   )
  ) 
  
  stopCluster(cl)
  
  # collect results here
  result <- list(initial_pop = init_guess_grid, 
                 DEobj = DEobj, 
                 Hypothesis = hypo_name)
  
  # written the result
  result  
  
}

# this convenience function allows to change one or more values of a named vector  
# useful for leaky2 'p_intro' parameter 

set_param <- function(default_param = param_vals_est, 
                      param_to_change = "p_intro", new_val) {
  
  new_param_v <- default_param
  new_param_v[param_to_change] <- new_val
  
  new_param_v
  
}




