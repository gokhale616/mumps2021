# # contact matrix
contact_plt <- (
  contact_matrix %.>%
    plot_contact_matrix(contact_matrix = .) +
    project_theme +
    theme(legend.position = "bottom", 
          text = element_text(size = 20)) +
    cap_axes +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top"))
)

# load all the treated covariates 
source("../fit/treat_vacc_covar.R", chdir = TRUE)

# anno layer to mark the where the case records begin
#"#B2022F"

anno_case_layer <- (
  annotate(geom = "rect", 
           xmin = 1977, xmax = 2018, 
           ymin = 0, ymax = 1, 
           fill = "#b91d73", alpha = 0.2)
  )

x_lims <- c(1950, 2020)
x_breaks <- seq(1950,2020, by = 10)

y_lims <- c(0, 1)
y_breaks <- seq(0, 1, by = 0.25)


interpolate_colour <- "#FF8C00"

default_colour <- "grey30"

# plot for the probability of reporting age
prob_report_plt <- (
  mod_mumps_covariates_slow %.>%
    filter(., year > 1949) %.>% 
    select(., year, eta_a) %.>% 
    ggplot(., aes(x = year, y = eta_a)) +
    geom_line(size = 0.8, colour = default_colour) +
    geom_point(pch = 21, fill = "white", size = 2, colour = default_colour) +
    labs(x = "", y = "Age-stratified\n Case Records  ") + #P(Age reported) Case reports\n  age recorded
    annotate(geom = "text", label = "Case Data", 
             angle = 90, x = 1974, y = 0.5) +
    anno_case_layer +
    scale_x_continuous(limits = x_lims, breaks = x_breaks) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks, 
                       labels = scales::percent) +
    project_theme +
    cap_axes
  )


# plot for normalized births
norm_births_plt <- (
  mod_mumps_covariates_slow %.>% 
    filter(., year > 1949) %.>% 
    select(., year, Births) %.>% 
    mutate(., 
           normalized_Births = Births/max(Births)) %.>% 
    ggplot(., aes(x = year, y = normalized_Births)) +
    geom_line(size = 0.8, colour = default_colour) +
    geom_point(pch = 21, fill = "white", size = 2, colour = default_colour) +
    labs(x = "", y = "Normalized Births    ") +
    anno_case_layer +
    scale_x_continuous(limits = x_lims, breaks = x_breaks) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks, 
                       labels = scales::percent) +
    project_theme +
    cap_axes
  )


# plot  neonatal cover
p1_cover_plt <- (
  mod_mumps_covariates_slow %.>% 
  filter(., year > 1949) %.>%   
  select(., year, p1) %.>% 
  mutate(., vacc_cover  = ifelse(year > 1967 & year < 1985, "Interpolated", "Observed")) %.>% 
  ggplot(., aes(x = year, y = p1, colour = vacc_cover)) +
  geom_line(aes(group = 1),  size = 0.8) +
  geom_point(pch = 21, fill = "white", size = 2)  +
  labs(x = "", y = "Neonatal Dose", colour = "Vaccine\nCoverage") +
  anno_case_layer +
  scale_x_continuous(limits = x_lims, breaks = x_breaks) +
  scale_y_continuous(limits = y_lims, breaks = y_breaks, 
                     labels = scales::percent) +
  scale_colour_manual(values = c(interpolate_colour, default_colour)) +
  project_theme +
  cap_axes +
  theme(legend.position = "none")  
  )




# plot booster cover
p2_plt_data <- (
  mod_mumps_covariates_slow %.>% 
  filter(., year > 1949) %.>% 
  mutate(., shape = "Slow") %.>% 
  bind_rows(., 
            mod_mumps_covariates_sigmoidal %.>% 
              mutate(., shape = "Sigmoidal")) %.>% 
  bind_rows(., 
            mod_mumps_covariates_rapid %.>% 
              mutate(., shape = "Rapid")) %.>%
  bind_rows(., 
            mod_mumps_covariates_constant %.>% 
              mutate(., shape = "Constant")) %.>%
  select(., year, p2, shape) %.>% 
  mutate(., 
         vacc_cover = ifelse(year > 1987 & year < 2000, "Interpolated", "Observed"))
  )


p2_cover_plt <- (
  p2_plt_data %.>% 
    ggplot(., aes(x = year, y = p2, colour = vacc_cover)) +
    geom_segment(aes(xend = if_else(lead(shape) == shape, lead(year), NA_integer_), 
                     yend = lead(p2)), size = 0.8) +
    geom_point(shape = 21, size = 2, fill = "white") +
    labs(x = "", y = "Booster Dose", color = "Vaccine\nCoverage") +
    anno_case_layer +
    scale_x_continuous(limits = x_lims, breaks = x_breaks) +
    scale_y_continuous(limits = y_lims, breaks = y_breaks, 
                       labels = scales::percent) +
    scale_colour_manual(values = c(interpolate_colour, default_colour)) +
    project_theme +
    cap_axes +
    guides(colour = guide_legend(direction = "horizontal", 
                                 nrow = 2))
  )
  
  

vacc_plt_legend <- get_legend(p2_cover_plt)  

p2_cover_plt <- p2_cover_plt + theme(legend.position = "none")

# age independent covariates

ai_grid_plot <- (
  plot_grid(prob_report_plt, norm_births_plt, 
            p1_cover_plt, p2_cover_plt, nrow = 2, 
            labels = c("A", "B", "C", "D"), align = "hv")
  )

ai_grid_plot_w_leg <- (
  plot_grid(ai_grid_plot, vacc_plt_legend, nrow = 2, rel_heights = c(1, 0.1))
)

# plot for age stratified covariates - population
pop_plt_data <- (
  mumps_covariates %.>% 
  select(., year, starts_with("N_")) 
  )

colnames(pop_plt_data) <- c("year", age_names)

pop_plt <- (
  pop_plt_data %.>% 
    filter(., year > 1949) %.>% 
    gather(., key = "age_cohort", value = "pop", -year, factor_key = TRUE) %.>% 
    ggplot(., aes(x = year, y = pop, fill = age_cohort)) +
    geom_area() +
    labs(x = "Year", y = "Population Size") +
    annotate(geom = "rect", 
             xmin = 1977, xmax = 2018, 
             ymin = 0, ymax = 4e8, 
             fill = "#b91d73", alpha = 0.2) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks) +
    scale_y_continuous(limits = c(0, 4e8), breaks = seq(0, 4e8, by = 1e8), 
                       labels = scales::scientific) +
    scale_fill_brewer(palette = "Greens", direction = -1) +
    project_theme +
    cap_axes +
    theme(legend.position = "none")  
    ) 
    
# plot for age stratified covariates - migartion  

mig_plt_data <- (
  mumps_covariates %.>% 
    select(., year, starts_with("MU_")) 
)

colnames(mig_plt_data) <- c("year", age_names)

mig_plt <- (
  mig_plt_data %.>% 
    filter(., year > 1949) %.>% 
    gather(., key = "age_cohort", value = "mig", -year, factor_key = TRUE) %.>% 
    ggplot(., aes(x = year, y = mig, fill = age_cohort)) +
    geom_area() +
    labs(x = "Year", y = "Migration Rate", fill = "Age\ncohort") +
    annotate(geom = "rect", 
             xmin = 1977, xmax = 2018, 
             ymin = -0.06, ymax = 0.06, 
             fill = "#b91d73", alpha = 0.2) +
    scale_x_continuous(limits = x_lims, breaks = x_breaks) +
    scale_y_continuous(limits = c(-0.06, 0.06), 
                       breaks = c(-0.06, -0.03, 0, 0.03, 0.06), 
                       labels = scales::scientific) +
    scale_fill_brewer(palette = "Greens", direction = -1) +
    project_theme +
    cap_axes +
    guides(fill = guide_legend(direction = "horizontal", nrow = 2))
)

as_plt_legend <- get_legend(mig_plt)

mig_plt <- mig_plt + theme(legend.position = "none")  

#age dependent covariates     

as_grid_plot <- plot_grid(pop_plt, mig_plt, nrow = 1, labels = c("E", "F"))

as_grid_plot_w_leg <- plot_grid(as_grid_plot, as_plt_legend, nrow = 2, rel_heights = c(1, 0.2))

# cobine_all_plots
all_cov_plt <- plot_grid(ai_grid_plot_w_leg, as_grid_plot_w_leg, nrow = 2, rel_heights = c(1, 0.6))

  

