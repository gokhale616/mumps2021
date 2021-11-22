# source prerequisites for this project
source("../00/src.R", chdir = TRUE)

# source all the tretaed covariates
source("../fit/treat_vacc_covar.R", chdir = TRUE)

po_to_vis <- make_pomp(., covar = mod_mumps_covariates_constant)

# spy(po_to_vis)

rp_vals_mod <- param_vals_est



mod_these_params <- c("q", "sigma", "beta1", "dwan", "epsilon1", "t_intro", "p_intro", "iota", "epsilon2")

rp_vals_mod[mod_these_params] <- c(0.22, 365.25/13, 0.11, Inf, 0, 2000, 2, 1, 0.5)

sample_traj <- trajectory(po_to_vis, param = rp_vals_mod, method = "ode23", format = "d")


prep_as_plot_trajectory <- function(traj_data, init_year = 1976) {

  sample_traj %.>% 
    filter(., year > init_year) %.>% 
    select(., -`.id`) %.>% 
    gather(., key = "comp", value = "count", -year) %.>% 
    mutate(., 
           comp_exp = str_split_fixed(comp, pattern = "_", n = 2)[,1], 
           age_number = str_split_fixed(comp, pattern = "_", n = 2)[,2], 
    ) %.>% 
    mutate(., 
           age_cohort = case_when(age_number == 1 ~ age_names[1], 
                                  age_number == 2 ~ age_names[2], 
                                  age_number == 3 ~ age_names[3], 
                                  age_number == 4 ~ age_names[4], 
                                  age_number == 5 ~ age_names[5], 
                                  TRUE            ~ age_number) %.>% as_factor(.)
    ) %.>% 
    select(., year, comp_exp, age_cohort, count)
    
}



sample_plot_traj <- sample_traj %.>% prep_as_plot_trajectory(.)

sample_plot_traj %.>% 
  ggplot(., aes(x = year, y = count, colour = age_cohort)) +
  geom_line(size = 0.8) +
  facet_wrap(vars(comp_exp), scales = "free") +
  project_theme +
  scale_color_brewer(palette = "Blues", direction = -1) +
  cap_axes 






  