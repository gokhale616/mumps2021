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

# ivp_names <- c("S_0", "E_0", "I_0", "B_0", "p_0")
# 
# fugazi_state0 <- rp_vals[ivp_names]



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





simultate_n_one_step <- function() {
  
  year_val <- mumps_weekly_case_reports$year
  
  na_string <- rep(NA, times = (length(year_val)-1))
  # browser()
  current_states <- (
    tibble(year = year_val,
           V_0 = 0,
           S_0 = c(round(1/10*219e6),         na_string),
           E_0 = c(round(0.0004003011*219e6), na_string), 
           I_0 = c(round(0.0001539356*219e6), na_string), 
           B_0 = 0, 
           p_0 = 0, 
           my_t = year_val)
  )
  
  final_statenames <- (
    colnames(current_states)[2:(ncol(current_states)-1)] %.>% 
      str_remove(., pattern = "_0")
    )
  
  
  for(i in 2:(length(na_string)+1)) {
    # browser()
    #po_1_step <- define_time_po_1_step(from = current_states$year[i-1], 
    #                                   to =   current_states$year[i])
    
    # browser()
    simulated_step <- (
      current_states %.>% 
        slice(., (i-1)) %.>% 
        unlist(.) %.>% 
        sim_p_vals(., default_p_vals = rp_vals) %.>% 
        simulate(po_1_step, param = ., include.data = FALSE, format = "d")
    )
    
    current_states[i,2:(ncol(current_states)-1)] <- (
      simulated_step %.>% 
        slice(., n()) %.>% 
        select(., final_statenames)
    )
    
  } 
  
  # browser()
  
  # reset colnames
  colnames(current_states) <- c("year", final_statenames, "my_t")
  
  current_states
}

# tic()
# simultate_n_one_step()
# toc()



test_one_step_sim <- (
  simultate_n_one_step() %.>% 
    gather(., key = "comp", value = "value", -year) %.>% 
    ggplot(., 
           aes(x = year, y = value)) +
    geom_line() +
    facet_wrap(.~comp, scales = "free")
  )









