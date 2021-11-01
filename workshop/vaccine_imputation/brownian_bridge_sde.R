source("../../00/src.R", chdir = TRUE)

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

rp_vals <- c(k = 0.2, `T` = 20, b = 0.85, sigma_wn = 365.25/100, hack = 0.45)


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


loop_over <- seq(0, 5, by = 0.1)

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
    k < 1 ~ "convex", 
    k > 1 ~ "concave", 
    TRUE ~ "linear")
    )





bb_traj <- k_vals_c %.>% 
  filter(., state_var %in%c("year", "p"))


bb_traj %.>% 
  ggplot(.,
         aes(x = year, y = state_value, 
             colour = bridge_curvature, group = k)) +
  geom_line(size = 0.3) +
  labs(x = "Year", 
       y = "State value", 
       colour = "Bridge\nCurvature") +
  facet_wrap(.~state_var, scales = "free") +
  scale_y_continuous(breaks = seq(0, 0.85, by = .17))+
  scale_color_viridis_d() +
  project_theme +
  cap_axes





