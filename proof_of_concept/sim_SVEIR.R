source("./SVEIR.R")
# View the pomp object
#spy(sveir_po)       


# Simulating synthetic data 
sveir_sim <- (
  sveir_po %>% 
    simulate(nsim = 10, format="d", include.data=FALSE, seed = 986746881L)
) 



sveir_traj <- sveir_po %.>% trajectory(., method = "ode23", format ="d")

all_comp_o_plt <- (
  sveir_traj %.>%
    as_tibble(.) %.>%
    mutate(., N = S + E + I + R + V) %.>%
    select(., -`.id`) %.>%
    filter(., Year > 0) %.>%
    gather(., key = Compartment, value = Comp_count, -c(Year)) %>%
    ggplot(aes(x = Year, y = Comp_count)) +
    geom_line() +
    facet_wrap(.~Compartment, scales = "free") +
    project_theme+
    cap_axes()
  )
  

sveir_sim_C <- (
  sveir_sim %.>% 
    select(., Year, C, `.id`)
)

sveir_traj_C <- (
  sveir_traj %.>% 
    select(., Year, C)
)

sveir_sim_C_plt <- (
  sveir_sim_C %.>%
    filter(., Year > 0) %.>%
    ggplot(.) +
    geom_line(aes(x = Year, y = C, alpha = `.id`)) +
    geom_line(data = sveir_traj_C %.>% filter(., Year > 0),
              aes(x = Year, y = C), colour = "red", size = 0.8) +
    labs(y = "True Incidence")


)
  

# simulate only using the observation model 
sveir_traj_array <- sveir_det_po %.>% trajectory(., method = "ode23")

all_p_vals_mat <- c(rp_vals, ip_vals) %.>% param_replicate_matrix(., n = 10)

raw_rmeasure_sim <- (
  rmeasure(sveir_det_po, 
         x = sveir_traj_array, 
         times = time(sveir_det_po), 
         params = all_p_vals_mat) %.>% 
  do.call(rbind, lapply(1:10, function(x) {.[,x,]})) %.>% 
  t(.) %.>% 
  as_tibble(.)
  ) 

# set names and bring the data t the right format
colnames(raw_rmeasure_sim) <- 1:10
det_rmeasure_sim_long <- (
  raw_rmeasure_sim %.>% 
    mutate(., Year = seq(0,45)) %.>% 
    gather(., key = `.id`, value = Cases, -Year) %.>% 
    mutate(., 
           `.id` = as_factor(`.id`), 
           model = "Observation\nNoise Only"
           )
)

synthetic_data <-(
  sveir_sim %.>% 
    select(., Year, `.id`, Cases) %.>% 
    mutate(., 
           Year = as.integer(Year), 
           `.id`= as.integer(`.id`),
           model = "Process\n And Observation\nNoise"
           ) %.>%
    as_tibble(.) %.>% 
    mutate(., `.id`= as_factor(`.id`)) %.>% 
    bind_rows(., det_rmeasure_sim_long) %.>% 
    arrange(., `.id`) %.>% 
    filter(., Year > 0)
      
)

# sveir_traj_Cs <- (
#   sveir_traj_C %.>% 
#     filter(., Year > 0) %.>% 
#     mutate(., 
#            Cs = C*1e-3)
# )


# synthetic_data %.>%
#   ggplot(.)+
#   geom_line(aes(x = Year, y = Cases, group = `.id`, alpha = `.id`))+
#   geom_line(data = sveir_traj_Cs,
#             aes(x = Year, y = Cs), colour = "red", size = 0.3)+
#   facet_grid(rows = vars(model)) +
#   scale_x_continuous(breaks = c(1,seq(9, 45, by = 9))) +
#   project_theme+
#   cap_axes()

# save data  in the prcessed data folder 
save(synthetic_data, file = "./synthetic_data.rds")
















