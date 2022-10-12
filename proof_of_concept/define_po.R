# prepare the data 
# time series with only observation noise
o_noisy_data <- (
  synthetic_data %.>%
    filter(., model == "Observation\nNoise Only" & `.id` == 5) %.>% 
    select(., Year, Cases) %.>% 
    mutate(., Year = as.numeric(Year)) %.>% 
    bind_rows(., 
              tibble(Year = 0, Cases = NA)) %.>% 
    arrange(., Year)
)

# time series with observation and process noise
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
# pomp object with just observation noise
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

# pomp object with observation and process noise
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
