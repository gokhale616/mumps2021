# depends on ./SVEIR.R
load("./synthetic_data.rds")

# prepare the data 
o_noisy_data <- (
  synthetic_data %.>%
    filter(., model == "Observation\nNoise Only" & `.id` == 5) %.>% 
    select(., Year, Cases) %.>% 
    mutate(., Year = as.numeric(Year)) %.>% 
    bind_rows(., 
              tibble(Year = 0, Cases = NA)) %.>% 
    arrange(., Year)
)
  

op_noisy_data <- (
  synthetic_data %.>%
    filter(., model != "Observation\nNoise Only" & `.id` == 5) %.>% 
    select(., Year, Cases) %.>% 
    mutate(., Year = as.numeric(Year)) %.>% 
    bind_rows(., 
              tibble(Year = 0, Cases = NA)) %.>% 
    arrange(., Year)
)


# prepare the pomp objects
o_noisy_po <- (
  o_noisy_data %.>% 
    pomp(., 
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

op_noisy_po <- (
  op_noisy_data %.>% 
    pomp(., 
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


#o_noisy_po %.>% spy(.)
#op_noisy_po %.>% spy(.)
est_these <- c("sigma", "R0", "beta1", "dwan", "rho", "psi")
# form objective functions 
o_noisy_objfun <- (
  traj_objfun(o_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
  )

op_noisy_objfun <- (
  traj_objfun(op_noisy_po,
              params = c(rp_vals, ip_vals),
              est = est_these, 
              fail.value = 1e20, 
              ode_control = list(method = "ode23"))
)


# objective function evaluation at the MLE 
o_noisy_objfun(rp_vals[est_these])
op_noisy_objfun(rp_vals[est_these])

# running a local search
o_noisy_optim <- (
  optim(par = rp_vals[est_these], 
        fn = o_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
        )
  )

save(o_noisy_optim, file = "o_noisy_optim.rds")
op_noisy_optim <- (
  optim(par = rp_vals[est_these], 
        fn = op_noisy_objfun,
        control = list(trace = 10, 
                       maxit = 1e4)
  )
)
save(op_noisy_optim, file = "op_noisy_optim.rds")

load("o_noisy_optim.rds")
load("op_noisy_optim.rds")
nsim = 1000
# simulate from the observation model - only observation noise
o_est <- sim_p_vals(o_noisy_optim$par, c(rp_vals, ip_vals))

o_traj_array <- o_noisy_po %.>% trajectory(., params = o_est, method = "ode23")

all_o_vals_mat <- o_est %.>% param_replicate_matrix(., n = nsim)

o_raw_rmeasure_sim <- (
  rmeasure(o_noisy_po, 
           x = o_traj_array, 
           times = time(o_noisy_po), 
           params = all_o_vals_mat) %.>% 
    do.call(rbind, lapply(1:nsim, function(x) {.[,x,]})) %.>% 
    t(.) %.>% 
    as_tibble(.)
) 

# set names and bring the data t the right format
colnames(o_raw_rmeasure_sim) <- 1:nsim
o_det_rmeasure_sum <- (
  o_raw_rmeasure_sim %.>% 
    mutate(., Year = seq(0,45)) %.>% 
    gather(., key = `.id`, value = Cases, -Year) %.>% 
    mutate(., 
           Cases = ifelse(Year == 0, NA, Cases),
    ) %.>% 
    select(., -`.id`) %.>% 
    group_by(., Year) %.>% 
    dplyr::summarise(., 
                     qs = quantile(Cases, c(0.025, 0.5, 0.975), na.rm = TRUE), 
                     prob = c("0.2", "0.5", "0.8"), 
                     .groups = 'drop') %.>%
    spread(., key = prob, value = qs) %.>% 
    ungroup(.) %.>% 
    mutate(., model = "Observation\nNoise Only")
)


# simulate from the observation model - observation and process noise
op_est <- sim_p_vals(op_noisy_optim$par, c(rp_vals, ip_vals))

op_traj_array <- op_noisy_po %.>% trajectory(., params = op_est, method = "ode23")

all_op_vals_mat <- op_est %.>% param_replicate_matrix(., n = nsim)

op_raw_rmeasure_sim <- (
  rmeasure(op_noisy_po, 
           x = op_traj_array, 
           times = time(op_noisy_po), 
           params = all_op_vals_mat) %.>% 
    do.call(rbind, lapply(1:nsim, function(x) {.[,x,]})) %.>% 
    t(.) %.>% 
    as_tibble(.)
) 

# set names and bring the data t the right format
colnames(op_raw_rmeasure_sim) <- 1:nsim

op_det_rmeasure_sum <- (
  op_raw_rmeasure_sim %.>% 
    mutate(., Year = seq(0,45)) %.>% 
    gather(., key = `.id`, value = Cases, -Year) %.>% 
    mutate(., 
           Cases = ifelse(Year == 0, NA, Cases),
    ) %.>% 
    select(., -`.id`) %.>% 
    group_by(., Year) %.>% 
    dplyr::summarise(., 
              qs = quantile(Cases, c(0.025, 0.5, 0.975), na.rm = TRUE), 
              prob = c("0.2", "0.5", "0.8"), 
              .groups = 'drop') %.>%
    spread(., key = prob, value = qs) %.>% 
    ungroup(.) %.>% 
    mutate(., model = "Process\n And Observation\nNoise")
)


# combine the synthetic data for plotting
synthetic_data_for_plot <- (
  synthetic_data %.>% 
    filter(., `.id` == 5) %.>% 
    select(., -.id)
  )

# plot for the fit 

proof_of_concept_plt <- (
  op_det_rmeasure_sum %.>% 
    bind_rows(., 
              o_det_rmeasure_sum) %.>% 
    ggplot(., aes(x = Year))+
    geom_line(aes(y = `0.5`, colour = "Model\nSimulation")) +
    geom_ribbon(aes(ymin = `0.2`, ymax = `0.8`, fill = "95% Prediction\nInterval"), alpha = 0.5)+
    geom_line(data = synthetic_data_for_plot, aes(y = Cases, colour = "Synthetic\nData"), size = 0.8)+
    labs(y = "Cases")+
    facet_grid()+
    facet_grid(rows = vars(model)) +
    scale_x_continuous(breaks = c(1,seq(9, 45, by = 9))) + 
    scale_colour_manual(name = "Trajectory",
                        values = c("Median" = "Black", "Synthetic\nData" = "Red"))+
    scale_fill_manual(name = "", values = c("95% Prediction\nInterval" = "Black"))+
    project_theme+
    cap_axes() +
    guides(colour = guide_legend(nrow = 2))
  )




















