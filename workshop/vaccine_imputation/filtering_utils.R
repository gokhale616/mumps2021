source("./VSEI_1_step.R", chdir = TRUE)
library("lubridate")
library("future.apply")
# load and process data for vaccine imputation - 
# this was done here because we want to use the time variable will be 
# used in few of the function definitions below
mumps_weekly_case_reports <- (
  read_csv("../../raw_data/tycho_20201118-150400.csv") %.>%
    # select relevant variables for forming the time series
    select(.,
           PeriodStartDate, 
           PeriodEndDate,
           CountValue) %.>% 
    # find week midpoints for plotting convenience
    mutate(.,  
           PeriodMidDate = PeriodStartDate + (PeriodEndDate-PeriodStartDate)/2
    ) %.>% 
    select(., 
           -c(PeriodStartDate, PeriodEndDate)) %.>% 
    # summarize case report volume over all states
    group_by(., 
             PeriodMidDate) %.>% 
    summarise(., cases = sum(CountValue, na.rm = TRUE)) %.>% 
    ungroup(.) %.>% 
    arrange(., 
            PeriodMidDate) %.>% 
    # define a period to run the analysis
    filter(., 
           PeriodMidDate > as.Date("1967-12-01") & PeriodMidDate < as.Date("1985-01-01")) %.>% 
    # define a year variable for pomp
    mutate(., 
           year_val = year(PeriodMidDate), 
           week_val = week(PeriodMidDate), 
           year = (year_val-1968)+week_val/52
    ) %.>%  
    select(., PeriodMidDate, year, cases) %.>% 
    # add a missing row for pomp
    bind_rows(., 
              tibble(PeriodMidDate = NA, 
                     year = 0, 
                     cases = NA)
    )  %.>% 
    arrange(., year) %.>% 
    mutate(., 
           year = seq(0, 17, length.out = nrow(.)))
)



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
                                  n_particles = 10) {
  
  param_replicate_matrix(param_v = updated_states, n = n_particles) %.>% 
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


simultate_n_one_step <- function(n_particles) {
  
  # use this to generate the time variable to provide to the one step simulator
  year_val <- mumps_weekly_case_reports$year
  
  na_string <- rep(NA, times = (length(year_val)-1))
  
  # Generate a data frame of current states.
  # it is required to have a row of initial conditions 
  current_states <- (
    tibble(year = year_val,
           V_0 = c(0,                         na_string),
           S_0 = c(round(1/10*219e6),         na_string),
           E_0 = c(round(0.0004003011*219e6), na_string), 
           I_0 = c(round(0.0001539356*219e6), na_string), 
           B_0 = c(0,                         na_string), 
           p_0 = c(0,                         na_string), 
           obs_t = year_val, 
           `.id` = 0)
  )
  
  # names used for the outcome data-frame
  final_statenames <- (
    colnames(current_states)[2:(ncol(current_states)-1)] %.>% 
      str_remove(., pattern = "_0")
    )

  # across the time, check the integrate the following steps 
  for(i in 2:(length(na_string)+1)) {
    
    browser()
    
    # extract previousl states 
    states_vector <- (
      current_states %.>% 
        slice(., (i-1)) %.>% 
        unlist(.)
      )
    
    # generate the state-param array given a state vector 
    param_states <- gen_param_state_vals(updated_states = states_vector, 
                                         n_particles = n_particles)
    
    
    
    
    simulated_step <- (
      current_states %.>% 
        slice(., (i-1)) %.>% 
        unlist(.) %.>% 
        sim_p_vals(., default_p_vals = rp_vals) %.>% 
        simulate(po_1_step, param = ., include.data = FALSE, format = "d") %.>% 
        slice(., n()) %.>% 
        select(., 
               final_statenames %.>% 
                 all_of(.)
        )
    )
    
    current_states[i,2:(ncol(current_states)-1)] <- simulated_step
    
    
    
  } 
  
  # reset colnames
  colnames(current_states) <- c("year", final_statenames, "obs_t")
  
  current_states
}








test_one_step_sim <- (
  simultate_n_one_step() %.>% 
  gather(., key = "comp", value = "value", -year) 
)  
  


if(FALSE) {
tic()
  test_one_step_sim %.>% 
    mutate(., 
          value =  ifelse(.$comp %in% c("V", "S", "E", "I"), value/rp_vals["pop"], value)
          ) %.>% 
    ggplot(., 
           aes(x = year, y = value)) +
    geom_line() +
    facet_wrap(.~comp, scales = "free") +
    scale_x_continuous(breaks = c(0 ,6, 12, 18), limits = c(0, 18)) +
    project_theme +
    cap_axes

toc()
}







