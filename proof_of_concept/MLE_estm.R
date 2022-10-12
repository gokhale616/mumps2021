# loads the skin for the pomp objects 
source("./SVEIR.R")

# simulates and saves synthetic data for the study
source("./sim_SVEIR.R")

# defines pomp objects for the study 
source("./define_po.R")


#o_noisy_po %.>% spy(.)
#op_noisy_po %.>% spy(.)

est_these <- c("sigma", "R0", "beta1", "dwan", "rho", "psi")
# form objective functions 
# objfun with just observation noise 
o_noisy_objfun <- (
  traj_objfun(o_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
  )

# objfun with observation and process noise 
op_noisy_objfun <- (
  traj_objfun(op_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
)


# objective function evaluation at the MLE 
# o_noisy_objfun(rp_vals[est_these])
# op_noisy_objfun(rp_vals[est_these])

# running a local search
o_noisy_optim <- (
  optim(par = rp_vals[est_these], 
        fn = o_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
        )
  )

save(o_noisy_optim, file = "./o_noisy_optim.rds")
op_noisy_optim <- (
  optim(par = rp_vals[est_these], 
        fn = op_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
  )
)
save(op_noisy_optim, file = "./op_noisy_optim.rds")

