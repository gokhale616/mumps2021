pof_dir <- "../proof_of_concept/"

load(paste0(pof_dir,"synthetic_data.rds"))
load(paste0(pof_dir, "o_noisy_optim.rds"))
load(paste0(pof_dir,"op_noisy_optim.rds"))
load(paste0(pof_dir,"o_noisy_bootstrap_res.rds"))
load(paste0(pof_dir, "op_noisy_bootstrap_res.rds"))


source(paste0(pof_dir, "SVEIR.R"))
source(paste0(pof_dir, "define_po.R"))

# df of mMLES
pof_param_ests <- (
  c(sigma = 365.25/25, R0 = 13, beta1 = 0.11, dwan = 100, rho = 1e-3, psi = 0.1) %.>% 
    rbind(., 
              op_noisy_optim$par 
              ) %.>%  
    rbind(., o_noisy_optim$par) %.>% 
    as_tibble(.) %.>% 
    mutate(., sigma = 365.25/sigma) %.>% 
    gather(., key = "Parameter", value = "Estimate") %.>% 
    mutate(., 
           model = rep(c("True value", "Observation\nNoise Only", 
                         "Process\n And Observation\nNoise"), 6)) %.>% 
    spread(., key = model, value = Estimate) %.>% 
    slice(., 4, 6, 1, 2, 5, 3) %.>% 
    mutate(., 
           Parameter = case_when(Parameter == "R0"~"$R_0$",
                                 Parameter == "beta1"~"$\\beta_1$", 
                                 Parameter == "sigma"~"$\\sigma^{-1}$ (Days)", 
                                 Parameter == "dwan"~ "$\\delta^{-1}$ (Years)",  
                                 Parameter == "rho"~ "$\\rho$", 
                                 Parameter == "psi"~ "$\\psi$" 
                                 )
           ) %>% 
    gather(., key = "model", value = "mle", -Parameter) %>% 
    mutate(mle = round(mle,5))
  )


pof_param_ci <- (
  o_noisy_bootstrap_res %.>% 
  mutate(., model = "Observation\nNoise Only") %.>% 
  rbind(., 
        op_noisy_bootstrap_res %.>% 
          mutate(., model = "Process\n And Observation\nNoise")
        ) %.>% 
  select(., 
         -c(`.id`, iter)
         ) %>% 
  mutate(., sigma = 365.25/sigma) %.>% 
  gather(., key = "Parameter", value = "Estimate", -model) %.>% 
  filter(., Estimate >= 0 ) %.>%   
  group_by(., model, Parameter) %>% 
  summarise(., 
            qs = quantile(Estimate, c(0.025, 0.975), na.rm = TRUE), 
            prob = c("0.025", "0.975"), 
            .groups = 'drop') %.>%
  ungroup(.) %>% 
  mutate(., 
         Parameter = case_when(Parameter == "R0"~"$R_0$",
                               Parameter == "beta1"~"$\\beta_1$", 
                               Parameter == "sigma"~"$\\sigma^{-1}$ (Days)", 
                               Parameter == "dwan"~ "$\\delta^{-1}$ (Years)",  
                               Parameter == "rho"~ "$\\rho$", 
                               Parameter == "psi"~ "$\\psi$" 
         ),
         qs = round(qs,5),
         ) %.>%
    spread(., key = prob, value= qs) %.>% 
    mutate(., ci = paste0("(",`0.025`,", ", `0.975`,")")) %>% 
    select(., 1, 2, 5)  
    )
  
  
pof_param_tbl_data <- (
  pof_param_ests %.>% 
    full_join(.,
              pof_param_ci,
              by = c("Parameter", "model")) %.>% 
    mutate(., 
           mle_ci = ifelse(is.na(ci) == TRUE, mle, paste(mle,ci))
           ) %.>% 
    select(., 1, 2, 5) %.>% 
    spread(., key = model, value = mle_ci) %.>% 
    select(., c(1,4,2,3)) %.>% 
    slice(., c(6, 1, 5, 2, 4, 3))
)


pof_param_ests_kbl <- (
  pof_param_tbl_data %.>% 
    kbl(., 
        align = "c", 
        digits = 4, 
        linesep = "",
        booktabs = T, 
        format = "latex", 
        caption = "Model specific parameter estimates were obtained by maximizing the likelihood function", 
        escape = FALSE) %.>% 
    kable_styling(.,
                  position='left', full_width = F,
                  latex_options=c('striped', 'HOLD_position', "scale_down")) %.>%
    add_header_above(., 
                     bold = TRUE, 
                     c(" " = 2,  "Model" = 2))
  )




nsim = 1000
set.seed(986747881L)
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
    geom_ribbon(aes(ymin = `0.2`, ymax = `0.8`, fill = "Prediction\nInterval"), alpha = 0.5)+
    geom_line(data = synthetic_data_for_plot, aes(y = Cases, colour = "Synthetic\nData"), size = 0.8)+
    labs(y = "Cases")+
    facet_grid()+
    facet_grid(rows = vars(model)) +
    scale_x_continuous(breaks = c(1,seq(9, 45, by = 9))) + 
    scale_colour_manual(name = "Trajectory",
                        values = c("Median" = "Black", "Synthetic\nData" = "Red"))+
    scale_fill_manual(name = "", values = c("Prediction\nInterval" = "Black"))+
    project_theme+
    cap_axes() +
    guides(colour = guide_legend(nrow = 2))
)




















