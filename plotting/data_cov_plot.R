# plot incidence data 
mumps_case_reports_l <- (
  mumps_case_reports %.>%
    gather(., key = "age_cohort", value = "cases", -c(year, total)) %.>% 
    mutate(., 
           sqrt_cases = sqrt(cases),
           age_cohort = as_factor(age_cohort))
    
) 

abs_inc_plt <- (
  mumps_case_reports_l %.>% 
  ggplot(., 
         aes(x = year, y = sqrt_cases, fill = age_cohort)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = expression(sqrt(Cases)),
       x = "") +
  scale_fill_manual(values = age_cols) +
  scale_x_continuous(expand = c(0.005, 0), 
                     breaks = seq(1977, 2018, by = 8)) +  
  scale_y_continuous(expand = c(0.005, 0), 
                     lim = c(0, 250), 
                     breaks = seq(0, 250, 125)) +
  project_theme + 
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()
        )
  )


prop_inc_plt <- (
  mumps_case_reports_l %.>% 
  ggplot(., 
         aes(x = year, y = cases, fill = age_cohort)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Year", 
       y = "Percent of total cases", 
       fill = "Age\nCohort") +
  scale_fill_manual(values = age_cols) +
  scale_x_continuous(expand = c(0.005, 0), breaks = seq(1977, 2018, by = 8)) +
  scale_y_continuous(expand = c(0.005, 0), 
                     breaks = c(0, 0.5, 1),
                     labels = scales::percent_format(accuracy = 2)) +  
  project_theme
  )


inc_legend <- prop_inc_plt %.>% get_legend(.) 

prop_inc_plt <- prop_inc_plt + theme(legend.position = "none") 

incidence_plt <- (
  plot_grid(abs_inc_plt, prop_inc_plt,
            nrow = 2, labels = c("A", "B"),
            align = "hv", axis = "b")
  ) 

# final incidence plot
incidence_w_legend_plt <- (
  plot_grid(incidence_plt, inc_legend, 
            rel_heights = c(1, 0.1), nrow = 2)
  ) 


incidence_line_plt <- mumps_case_reports_l %.>% 
  ggplot(., aes(x = year, y = sqrt_cases)) +
  geom_line(size = 0.8) +
  labs(y = expression(sqrt(Cases)), 
       x = "Year") +
  facet_grid(rows = vars(age_cohort), scales = "free_y") +
  project_theme
  



# contact matrix
contact_plt <- (
  contact_matrix %.>% 
  plot_contact_matrix(contact_matrix = .) +
  scale_x_discrete(expand = c(0, 0.005)) +
  scale_y_discrete(expand = c(0, 0.005)) +
  project_theme +
  theme(legend.position = "bottom")
  ) 
  

# co-variate data panel
mumps_covariates_plt_data <- plot_covars(mumps_covariates)








  
  