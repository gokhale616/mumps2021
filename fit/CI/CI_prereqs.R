# load all the prereq functions 
source("../../00/src.R", chdir = TRUE) 

# load all MLES coresponding to all the models 
source("../../plotting/prep_estm_tables.R", chdir = TRUE)

# select the right parameter estimates
all_result_df %.>% 
  select(., -c(R0, Rp, impact, best_fit_covar)) %.>% 
  group_by(., hypothesis) %.>% 
  filter(., d_AIC == min(d_AIC)) %.>% 
  ungroup(.)


# prepare simulated trajectories





