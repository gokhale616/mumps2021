source("../../00/src.R", chdir = TRUE)

library(zoo)

bb_stoc_stop <- Csnippet(
  "
  // make normal random draws for the Brownian bridge
  double dw = rnorm(0, dt*sigma_wn);
  
  // define a temporal boundary to make sure that bridge state doen not update before t = 0. 
  // this will avoid eccentric evaluations
  
  double r_treat;

  if (t < 0) {
    B += 0;
    r_treat = 0;
  } else {
      // update the brownian bridge here   
      B += ((-B)/(T-t) + dw); 
  
      // add the curvature 
      r_treat = (pow((t/T), k) + B); 
      
  }
  
  
  // add a couple of conditions to make sure that the treated bb stays within the [0,b] bounds
  // maybe a logit max function can be used here??!
  
  if((T-t) < hack || r_treat > b) {
    r =  b;
  } else if (r_treat < 0) {
      r =  0;
  } else {
      r = r_treat;  
  }
  
   
"
)



bb_rmeasure  <- Csnippet( 
"
p = r;
"
)

bb_init <- Csnippet( 
"
B = 0;
r = 0;

"
)


state_names <- c("B", "r")

rp_names <- c("k", "T", "b", "sigma_wn", "hack")

rp_vals <- c(k = 0.2, `T` = 20, b = 0.85, sigma_wn = 365.25/30, hack = 0.45)


bb_po <- (
  tibble(year = seq(-10, 20, by = 1/52), p = NA) %.>% 
    pomp(.,
         t0 = -50,
         times = "year",
         rinit    = bb_init, 
         rprocess = euler(bb_stoc_stop, delta.t = 1/365.25),
         rmeasure = bb_rmeasure,
         accumvars = "B", 
         statenames = state_names, 
         paramnames = rp_names,
         params = rp_vals)
  )


loop_over <- c(0.15, 0.35, 0.65, 1, 2.0, 3, 6)#seq(0.00, 5, length.out = 4)

k_vals <- map_dfr(1:length(loop_over), 
                  function(c, k_vec = loop_over){
  
  rp_vals_int <- rp_vals
  
  rp_vals_int["k"] <- k_vec[c] 
  
  bb_po %.>% 
    simulate(., nsim = 1, param = rp_vals_int, 
             format="d", include.data=FALSE) %.>% 
    mutate(., 
           B = ifelse(B > 0.1|B < -0.1, 0, B)) %.>% 
    gather(., 
           key = "state_var", value = "state_value", 
           factor_key = TRUE, -c(year, .id)) %.>% 
    mutate(., 
           k = k_vec[c])
  
  
})


k_vals_c <- k_vals %.>% 
  mutate(., bridge_curvature = case_when(
    k < 1 ~ "Convex", 
    k > 1 ~ "Concave", 
    TRUE ~ "Linear")
    )





bb_traj <- k_vals_c %.>% 
  filter(., state_var %in%c("year", "p"))



bb_roll_median <- (
  bb_traj %.>% 
    select(., year, k, state_value) %.>%   
    group_by(., k) %.>%   
    summarize(.,
              med_year_value = rollmedian(year, k = 53),
              med_state_value = rollmedian(state_value, k = 53)*100) %.>% 
    ungroup(.) %.>% 
    filter(., 
           med_year_value >= 0)
    )
  

bb_traj %.>% 
  filter(., year >= 0) %.>% 
  ggplot(.) +
  geom_line(aes(x = year, y = state_value*100, 
                colour = bridge_curvature, group = k), size = 0.5) +
  geom_line(data = bb_roll_median, 
            aes(x = med_year_value, y = med_state_value, group = k, 
                linetype = "Annual"),
            colour = "#71B280", 
            size = 0.8)+
  labs(y = "Vaccine coverage (%)", 
       x = "Year", 
       colour = "Bridge\nCurvature", 
       linetype = "Simulated\nTrajectory") +
  #facet_wrap(.~state_var, scales = "free") +
  scale_y_continuous(breaks = seq(0, 85, by = 17), limits = c(0, 85))+
  scale_color_viridis_d(option = "F") +
  scale_linetype_manual(values = c("Weekly" = "solid", "Annual" = "dashed")) +
  project_theme +
  cap_axes() +
  guides(colour = guide_legend(title.position = "left", 
                               nrow = 2), 
         linetype = guide_legend(title.position = "left", 
                                 nrow = 2, 
                                 override.aes = list(colour = "black", 
                                                     size = 0.5)))





