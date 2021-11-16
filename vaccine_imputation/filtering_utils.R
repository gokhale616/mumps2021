set.seed(986747881L)

# this is a pomp object that simulates 1 step in the future and it also preserves the preceding time-step
define_time_po_1_step <- function(from = 0, to = 1/52) {

  tibble(year = c(from, to), 
         cases = NA) %.>% 
    make_pomp_vsei(data = .)
  

}
  
po_1_step <- define_time_po_1_step()


# this function converts updated states into parameter vectors and appends it with the randomly generated the 
# curvature parameter
gen_param_state_vals <- function(updated_states = fugazi_state0, 
                                 particles = 10) {

  updated_states %.>% 
    # mutate(., 
    #        k              = rtrunc(n(), spec = "norm", a = 0, b = 15, mean = 4, sd = 1), 
    #        sigma_wn_scale = rtrunc(n(), spec = "norm", a = 1, b = 30, mean = 17, sd = 1)
    #        ) 
    mutate(.,
           k              = runif(n(), min = 0, max = 5),
           sigma_wn_scale = runif(n(), min = 1, max = 50)
           )
}


# this is function simulates particles one step forward 
sim_particles_one_step <- function(state_param_matrix, 
                                   default_p_vals = rp_vals) {
  
  
  
  map_dfr(1:nrow(state_param_matrix),
          function(c) {
    # browser()
    # pulls and sets the the param vector to simulate from 
    particle_configuration <- (
      state_param_matrix %.>% 
        slice(., c) %.>% 
        unlist(.) %.>% 
        sim_p_vals(., default_p_vals = default_p_vals)
    )
    # browser()
    # simulates an updated state
    particle_configuration %.>% 
      simulate(object = po_1_step, 
               param = particle_configuration, 
               format = "d") %.>% 
      slice(., n()) %.>% 
      mutate(., 
             `.id` = c,
             obs_t = NA) %.>% 
      select(., -c(Reff, cases)) 
     
    
    
  }) 
  
}


# filtering distribution
over_dispersed_normal <- function(obs_cases, sim_true_cases, rho = 0.06, psi = 1.2) {
  
   # define mean and variance of the   standard normal distribution
   m <- sim_true_cases*rho 
   v <- m*(1 - rho + m*psi^2) 
  
  # a bit more book-keeping 
  tol <- 1.0e-18 
  
  
  dnorm(obs_cases, m, sqrt(v)+tol, log = TRUE)
  }


# this function normalizes the weights
normalize_particle_weights <- function(weights) {
  weights/sum(weights)
}



# this function re-samples from the particles using the filtering distribution
resample_particle_indices <- function(j_norm_weights) {
  
  # number of particles
  j_particles <- length(j_norm_weights)
  
  
  # cumulative sum over normalized weights (wn)
  j_cumsum_wn <- cumsum(j_norm_weights)
  
  # re-assignment of the index
  # initial draw from a uniform distribution
  u <- c(runif(1, max = 1/j_particles), rep(NA, times = (j_particles-1)))
  
  # generate evenly spaced sampling points 
  u <- map_dbl(1:j_particles, function(c) {ifelse(c == 1, u[1], u[1]+(c-1)*1/j_particles)})
  
  # produce a vector of resampling indices to populate
  k_j <- rep(NA, j_particles)
  # initialize the set of raw indices at one 
  p = 1  
  
  for(j in 1:j_particles) {
    
    while(u[j] > j_cumsum_wn[p]) {
      p <- p + 1
    }
    # assign the re-sampling index
    k_j[j] <- p
  }
  
  # return indices used to sample from the particle population
  k_j
  
  
}


# load initial conditions
load("./init_for_vsei.rds")

pfilter_once <- function(case_data = mumps_weekly_case_reports %.>% 
                           select(., year, cases), 
                         j_particles = 2, 
                         init_states = init_for_vsei %<>% 
                           unlist(.)) {
  
  
  # Generate a data frame of current states.
  # the first row of this data frame is always going to be initial conditions
  # we always start with 1 - for all particles
  current_states <- (
    tibble(year  = rep(case_data$year[1], times = j_particles),
           `.id` = 1:j_particles, 
           V_0   = init_states["V"],
           S_0   = init_states["S"],
           E_0   = init_states["E"], 
           I_0   = init_states["I"], 
           C_0   = 0, 
           B_0   = 0, 
           p_0   = 0, 
           obs_t = case_data$year[1])
  )
  
  
  for(n in 2:nrow(case_data)) {
      
    # add the parameters to simulate from 
    param_state <- (
      current_states  %.>%  
        filter(., year == case_data$year[(n-1)]) %.>% 
        gen_param_state_vals(., particles = j_particles)
    )
    
    
    # using these param-state values, simulate one step forward
    updated_states <- (
      param_state %.>% 
      sim_particles_one_step(.) %.>%  
      mutate(., 
             year = case_data$year[n], 
             obs_t = case_data$year[n]
      )
      )
    
    # add data and and generate weights
    weights <- (
      updated_states %.>%
      right_join(.,
                 case_data %.>% 
                   filter(., year == case_data$year[n]), 
                 by = "year"
                 ) %.>% 
      mutate(., 
             weights = over_dispersed_normal(cases, C)) %.>% 
      select(., weights) %.>% 
      unlist(.) %.>%   
      unname(.) 
      ) 

    # normalize weights
    norm_weigths <- normalize_particle_weights(weights) 
    
    # identify new indices
    filtered_indices <- resample_particle_indices(norm_weigths)
    
    # sampled using filtered indices
    resampled_updated_states <- (
      updated_states %.>% 
        slice(., 
              filtered_indices)
    ) 
    
    # change the colnames to match the colnames of 
    colnames(resampled_updated_states) <- colnames(current_states)
    
    # update current state 
    current_states <- (
      current_states %.>% 
        bind_rows(., 
                  resampled_updated_states)
    )
  }
  
  # reset col-names
  colnames(current_states) <- c("year", ".id", "V", "S", "E", "I", "C", "B", "p", "obs_t")
  
  current_states
  
}





