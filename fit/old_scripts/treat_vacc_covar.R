# prepare the vacc-coverage for p1
imputed_vacc_coverage <- (
  filtered_vacc %.>% 
    right_join(., 
               by = "year",
               mumps_weekly_case_reports %.>% 
                 select(., -cases))  %.>% 
    select(., PeriodMidDate, p) %.>% 
    group_by(., PeriodMidDate)  %.>% 
    summarize(., 
              median = median(p)) %.>% 
    ungroup(.)  %.>% 
    mutate(., 
           year = year(PeriodMidDate)
    ) %.>%
    group_by(., year) %.>% 
    summarise(., 
              p1 = mean(median)) %.>% 
    drop_na(.)  
)



set.seed(986747881)

# interpolate the covars such that vaccine uptake is ....

# ... slow ...
log_par_defaults_slow <- c(t0_p1 = 1979, b_p1 = 1, c_p1 = 0.5, d_p1 = 0.8, 
                           t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)

mod_mumps_covariates_slow <- treat_covar_data(par = log_par_defaults_slow) 

# plot_covars(covar_data = mod_mumps_covariates_slow, filter_from = 1965)


# ... sigmoidal ...
log_par_defaults_sigmoidal <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
                                t0_p2 = 1994, b_p2 = 1, c_p2 = 1, d_p2 = 0.8)

mod_mumps_covariates_sigmoidal <- treat_covar_data(par = log_par_defaults_sigmoidal) 

# plot_covars(covar_data = mod_mumps_covariates_sigmoidal, filter_from = 1965)


# ... rapid ...
log_par_defaults_rapid <- c(t0_p1 = 1968, b_p1 = 1, c_p1 = 0.5, d_p1 = 1.5, 
                            t0_p2 = 1989, b_p2 = 1, c_p2 = 0.5, d_p2 = 1.5)

mod_mumps_covariates_rapid <- treat_covar_data(par = log_par_defaults_rapid) 



# .... constant ......
mod_mumps_covariates_constant <- (
  mod_mumps_covariates_slow %.>% 
    select(., -p2) %.>% 
    right_join(., 
               mumps_covariates %.>% 
                 select(., year, p2), 
               by = "year")
)



# replace the p1 with the imputed data for all of the treated covars
add_imputed_p1 <- function(data) {rows_update(data, imputed_vacc_coverage, by = "year")}

mod_mumps_covariates_slow %<>% add_imputed_p1()
mod_mumps_covariates_sigmoidal %<>% add_imputed_p1()
mod_mumps_covariates_rapid %<>% add_imputed_p1()
mod_mumps_covariates_constant %<>% add_imputed_p1()












