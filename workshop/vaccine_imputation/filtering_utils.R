# make a one step ahead simulator 
source("./VSEI_1_step.R")

# this is a pomp object that simulates 1 step in the future and it also preserves the preceding time-step
po_1_step <- (
  tibble(year = seq(0, 1/52, 1/52), 
         cases = NA) %.>% 
    make_pomp_vsei(data = .)
  )

ivp_names <- c("S_0", "E_0", "I_0", "B_0", "p_0")

fugazi_state0 <- rp_vals[ivp_names]



# this function converts updated states into parameter vectors and appends it with the randomly generated the 
# curvature parameter
# 
gen_param_state_vals <- function(updated_states = fugazi_state0, 
                                  n_particles = 10) {
  
  param_replicate_matrix(updated, n = n_particles) %.>% 
    t(.) %.>% 
    as_tibble(.) %.>% 
    mutate(., k = rtrunc(n(), spec = "norm", a = 0, b = 15, mean = 4, sd = 1)) 
  
}


# this is function simulates particles one step forward 
sim_particles_one_step <- function(state_param_array, 
                                   default_p_vals = rp_vals) {
  
  map_dfr(1:nrow(state_param_array), function(c) {
    
    # pulls and sets the the param vectory to simulate from 
    particle_configuration <- (
      state_param_array %.>% 
        slice(., c) %.>% 
        unlist(.) %.>% 
        sim_p_vals(., default_p_vals = default_p_vals)
    )
    
    # simulates an updated state
    particle_configuration %.>% 
      simulate(object = po_1_step, 
               param = particle_configuration, 
               format = "d") %.>% 
      slice(., n()) %.>% 
      mutate(., 
             particle_id = c, 
             k = particle_configuration['k'])
      
     
    
    
  })
  
}

#  a simple for loop to simulate until
for (i in 1:15) {
  
  
  
}































