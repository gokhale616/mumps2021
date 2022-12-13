source("./SVEIR.R")
# View the pomp object
#spy(sveir_po)       

rm(filtered_vacc)
rm(mumps_case_reports)
rm(mumps_demog_geog_annual)
rm(mumps_covariates)
rm(mumps_weekly_case_reports)
rm(C)
rm(contact_matrix)


# load the mle values 
load("./o_noisy_optim.rds")
load("./op_noisy_optim.rds")


est_these <- c("sigma", "R0", "beta1", "dwan", "rho", "psi")

n_sim <- 5000


# extract parameter values 
o_noisy_mle <- c(rp_vals, ip_vals)
o_noisy_mle[est_these] <- o_noisy_optim$par

op_noisy_mle <- c(rp_vals, ip_vals)
op_noisy_mle[est_these] <- op_noisy_optim$par


# Simulating synthetic data for o_noisy mle
o_noisy_traj <- (
  sveir_det_po %.>% 
    trajectory(., 
               params = o_noisy_mle, 
               method = "ode23")
)


o_noisy_p_mat <- o_noisy_mle %.>% param_replicate_matrix(., n = n_sim)

o_noisy_mle_ts <- (
  rmeasure(sveir_po, 
           x = o_noisy_traj, 
           times = time(sveir_po), 
           params = o_noisy_p_mat) %.>% 
    do.call(rbind, lapply(1:n_sim, function(x) {.[,x,]})) %.>% 
    t(.) %.>% 
    as_tibble(.) %>% 
    mutate(., Year = seq(0,45)) %.>% 
    gather(., key = `.id`, value = Cases, -Year) %.>% 
    mutate(., 
           #`.id` = as_factor(str_split(`.id`, "")[[1]][2]),
           Cases = ifelse(Year == 0, NA, Cases)
    )
) 

#save(o_noisy_mle_ts, file = "./o_noisy_mle_ts.rds")

# Simulating synthetic data for op_noisy mle
op_noisy_traj <- (
  sveir_det_po %.>% 
    trajectory(., 
               params = op_noisy_mle, 
               method = "ode23")
)


op_noisy_p_mat <- op_noisy_mle %.>% param_replicate_matrix(., n = n_sim)

op_noisy_mle_ts <- (
  rmeasure(sveir_po, 
           x = op_noisy_traj, 
           times = time(sveir_po), 
           params = op_noisy_p_mat) %.>% 
    do.call(rbind, lapply(1:n_sim, function(x) {.[,x,]})) %.>% 
    t(.) %.>% 
    as_tibble(.) %>% 
    mutate(., Year = seq(0,45)) %.>% 
    gather(., key = `.id`, value = Cases, -Year) %.>% 
    mutate(., 
           #`.id` = as_factor(str_split(`.id`, "")[[1]][2]), 
           Cases = ifelse(Year == 0, NA, Cases)
    )
) 




# a function to parallelize parameteric bootstrap

optim_protocol <- function(counter, ts_data, mle) {
  # browser()
  # extract the right time series
  ts_int <- (
    ts_data %.>% 
      filter(., `.id` == sprintf("V%d", counter)) %.>% 
      select(., -`.id`) 
  )
  
  # define the pomp object
  po_int <- (
    ts_int %.>% 
      pomp(data = ., 
           t0 = -1000, 
           times = "Year",
           skeleton = vectorfield(Csnippet(sveir_skel)),
           rinit = Csnippet(sveir_rinit), 
           dmeasure = Csnippet(sveir_dmeas), 
           rmeasure = Csnippet(sveir_rmeas), 
           accumvars = zero_names, 
           statenames = state_names, 
           paramnames = c(rp_names, ip_names),
           params = c(rp_vals, ip_vals)
      )
  )
  
  # define the objecive function
  po_objfun_int <- (
    traj_objfun(po_int,
                params = c(rp_vals, ip_vals),
                est = est_these, 
                fail.value = 1e20, 
                ode_control = list(method = "ode23")
    )
  )
  
  
  # run Nelder-Mead with mle as the initial guess
  optim_int <- (
    optim(par = mle[est_these], 
          fn = po_objfun_int,
          control = list(trace = 10, 
                         maxit = 1e5)
    )
  )
  
  #extract the converged parameter value
  optim_int$par %.>% 
    as.list(.) %.>% 
    as_tibble(.) %.>% 
    mutate(., 
           `.id` = counter, 
           iter = optim_int$counts["function"] %.>% unname(.))
  
}



# testing the optim protocol function
# optim_protocol(counter = 500, ts_data = o_noisy_mle_ts, mle = o_noisy_mle)
# optim_protocol(counter = 517, ts_data = op_noisy_mle_ts, mle = op_noisy_mle)






