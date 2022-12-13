# loads the skin for the pomp objects 
source("./SVEIR.R")

# simulates and saves synthetic data for the study
source("./sim_SVEIR.R")

# defines pomp objects for the study 
source("./define_po.R")


#o_noisy_po %.>% spy(.)
#op_noisy_po %.>% spy(.)

est_these <- c("sigma", "R0", "beta1", "dwan", "rho", "psi")
trans_est_these <- parameter_trans(log = c("sigma", "R0", "dwan", "psi"), 
                                   logit = c("beta1", "rho"))


# use true values as initial guesses for faster convergence
# transform initialguesses these for unconstrained optimization
init_t <- rp_vals[est_these]
init_t[names(init_t) %in% c("beta1", "rho")] <- logit(init_t[names(init_t) %in% c("beta1", "rho")]) 
init_t[names(init_t) %nin% c("beta1", "rho")] <- log(init_t[names(init_t) %nin% c("beta1", "rho")]) 


# form objective functions 
# objfun with just observation noise 
o_noisy_objfun <- (
  traj_objfun(o_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              partrans = trans_est_these,
              paramnames = c(rp_names, ip_names),
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
)


# objfun with observation and process noise 
op_noisy_objfun <- (
  traj_objfun(op_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              partrans = trans_est_these,
              paramnames = c(rp_names, ip_names),
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
)


# objective function evaluation at the MLE 
# o_noisy_objfun(rp_vals[est_these])
# op_noisy_objfun(rp_vals[est_these])

# running a local search
o_noisy_optim <- (
  optim(par = init_t, 
        fn = o_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
        )
  )


op_noisy_optim <- (
  optim(par = init_t, 
        fn = op_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
  )
)

#if(FALSE){
# back transform parameters and save

# o_noisy
bck_t_invlogit_o_noisy <- invlogit(o_noisy_optim$par[names(o_noisy_optim$par) %in% c("beta1", "rho")]) 
bck_t_invlog_o_noisy <- exp(o_noisy_optim$par[names(o_noisy_optim$par) %nin% c("beta1", "rho")]) 

o_noisy_par <- o_noisy_optim$par

o_noisy_par[names(o_noisy_optim$par) %in% c("beta1", "rho")] <- bck_t_invlogit_o_noisy
o_noisy_par[names(o_noisy_optim$par) %nin% c("beta1", "rho")] <- bck_t_invlog_o_noisy

#save(o_noisy_par, file = "./o_noisy_par.rds")

# op_noisy
bck_t_invlogit_op_noisy <- invlogit(op_noisy_optim$par[names(op_noisy_optim$par) %in% c("beta1", "rho")]) 
bck_t_invlog_op_noisy <- exp(op_noisy_optim$par[names(op_noisy_optim$par) %nin% c("beta1", "rho")]) 

op_noisy_par <- op_noisy_optim$par

op_noisy_par[names(op_noisy_par) %in% c("beta1", "rho")] <- bck_t_invlogit_op_noisy
op_noisy_par[names(op_noisy_par) %nin% c("beta1", "rho")] <- bck_t_invlog_op_noisy



#save(op_noisy_par, file = "./op_noisy_par.rds")
}
