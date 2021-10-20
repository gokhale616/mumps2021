# soucrce the pre requisites 
source("../00/src.R", chdir = TRUE)

set.seed(986747881)

# interpolate the covars such that vaccine uptake is ....

# ... slow ...
log_par_defaults <- c(t0_p1 = 1979, b_p1 = 1, c_p1 = 0.5, d_p1 = 0.8, 
                      t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)

mod_mumps_covariates_slow <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_slow, filter_from = 1965)


# ... sigmoidal ...
log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
                      t0_p2 = 1994, b_p2 = 1, c_p2 = 1, d_p2 = 0.8)

mod_mumps_covariates_sigmoidal <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_sigmoidal, filter_from = 1965)


# ... rapid ...
log_par_defaults <- c(t0_p1 = 1968, b_p1 = 1, c_p1 = 0.5, d_p1 = 1.5, 
                      t0_p2 = 1989, b_p2 = 1, c_p2 = 0.5, d_p2 = 1.5)

mod_mumps_covariates_rapid <- treat_covar_data(par = log_par_defaults) 


# plot_covars(covar_data = mod_mumps_covariates_rapid, filter_from = 1965)


#  ... slow rapid ... 

log_par_defaults <- c(t0_p1 = 1979, b_p1 = 1, c_p1 = 0.5, d_p1 = 0.8, 
                      t0_p2 = 1989, b_p2 = 1, c_p2 = 0.5, d_p2 = 1.5)

mod_mumps_covariates_slow_rapid <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_slow_rapid, filter_from = 1965)

# ... slow sigmoid ...
log_par_defaults <- c(t0_p1 = 1979, b_p1 = 1, c_p1 = 0.5, d_p1 = 0.8, 
                      t0_p2 = 1994, b_p2 = 1, c_p2 = 1, d_p2 = 0.8)

mod_mumps_covariates_slow_sigmoid <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_slow_sigmoid, filter_from = 1965)

# .... slow constant ....

mod_mumps_covariates_slow_constant <- (
  mod_mumps_covariates_slow %.>% 
    select(., -p2) %.>% 
    right_join(., 
               mumps_covariates %.>% 
                 select(., year, p2), 
               by = "year")
)


# ... rapid slow ...
log_par_defaults <- c(t0_p1 = 1968, b_p1 = 1, c_p1 = 0.5, d_p1 = 1.5, 
                      t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)

mod_mumps_covariates_rapid_slow <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_rapid_slow, filter_from = 1965)


# ... rapid sigmoid ...
log_par_defaults <- c(t0_p1 = 1968, b_p1 = 1, c_p1 = 0.5, d_p1 = 1.5, 
                      t0_p2 = 1994, b_p2 = 1, c_p2 = 1, d_p2 = 0.8)

mod_mumps_covariates_rapid_sigmoid <- treat_covar_data(par = log_par_defaults)  


# plot_covars(covar_data = mod_mumps_covariates_rapid_sigmoid, filter_from = 1965)


# ... sigmoid slow ...
log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
                      t0_p2 = 1999, b_p2 = 1, c_p2 = 0.5, d_p2 = 0.8)

mod_mumps_covariates_sigmoidal_slow <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_sigmoidal_slow, filter_from = 1965)


# ... sigmoid rapid ...
log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.8, 
                      t0_p2 = 1989, b_p2 = 1, c_p2 = 0.5, d_p2 = 1.5)

mod_mumps_covariates_sigmoidal_rapid <- treat_covar_data(par = log_par_defaults) 

# plot_covars(covar_data = mod_mumps_covariates_sigmoidal_rapid, filter_from = 1965)







