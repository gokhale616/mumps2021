# treat incidence data for plotting

mumps_case_reports_l <- (
  mumps_case_reports %.>%
    gather(., key = "age_cohort", value = "cases", -c(year, total), factor_key = TRUE) %.>% 
    mutate(., 
           sqrt_cases = sqrt(cases), 
           sqrt_total = sqrt(total)
           )
    
) 

covars_for_case_conversion <- (
  mumps_covariates %.>% 
    mutate(., total_pop = N_1 + N_2 + N_3 + N_4 + N_5) %.>% 
    filter(., year > 1975 & year < 2019) %.>% 
    select(., year, total_pop)
  )

 

mumps_incidence_rate_l <- (
  mumps_case_reports_l %.>%
  filter(., age_cohort == "unknown") %.>%   
  select(., year, total) %.>% 
  right_join(., 
             covars_for_case_conversion, 
             by = c("year")) %.>% 
  mutate(., ac_inc_rate = (total/total_pop)*1e5) %.>%   
  select(., -c(total, total_pop))
  )
  



# y_max_abs <- max(mumps_case_reports_l$sqrt_total, na.rm = TRUE)
# 
# abs_inc_plt <- (
#   mumps_case_reports_l %.>% 
#     filter(., age_cohort == "unknown") %.>% 
#     ggplot(., 
#            aes(x = year, y = sqrt_total)) +
#     geom_area(stat = "identity", alpha = 0.9, fill = "#b21f1f") +
#     labs(y = expression(sqrt(Cases)),
#          x = "") +
#     scale_x_continuous(#expand = c(0.000, 0.000), 
#                        breaks = gen_x_breaks) + 
#     scale_y_continuous(#expand = c(0.0, 0), 
#                        breaks = (y_max_abs*c(0, 0.25, 0.5, 0.75, 1)) %.>% floor(.)) +
#     project_theme +
#     theme(axis.text.x = element_blank()) +
#     cap_axes +
#     theme(plot.margin = unit(rep(mar_val, 4), "cm")) 
#   )


# y_max_inc <- max(mumps_incidence_rate_l$ac_inc_rate, na.rm = TRUE)


inc_rate_plt <- (
  mumps_incidence_rate_l %.>% 
    ggplot(., 
           aes(x = year, y = ac_inc_rate)) +
    geom_area(stat = "identity", alpha = 0.9, fill = "#fcb045") +
    labs(y = expression(Cases~per~10^5),
         x = "") +
    scale_x_continuous(#expand = c(0.000, 0.000), 
      breaks = gen_x_breaks) + 
    scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), limits = c(0, 10)) +
    project_theme +
    theme(axis.text.x = element_blank()) +
    cap_axes +
    theme(plot.margin = unit(rep(mar_val, 4), "cm")) 
)


prop_inc_plt <- (
  mumps_case_reports_l %.>% 
  ggplot(., 
         aes(x = year, y = cases, fill = age_cohort)) +
  geom_area(stat = "identity", position = "fill") +
  labs(x = "Year", 
       y = "Percent of total cases", 
       fill = "Age\nCohort") +
  scale_fill_manual(values = age_cols) +
  scale_x_continuous(#expand = c(0, 0), 
                     breaks = gen_x_breaks) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent_format(accuracy = 1)) +  
  project_theme +
  cap_axes + 
  theme(plot.margin = unit(rep(mar_val, 4), "cm")) 
  )


inc_legend <- prop_inc_plt %.>% get_legend(.) 

prop_inc_plt <- prop_inc_plt + theme(legend.position = "none") 

incidence_plt <- (
  plot_grid(inc_rate_plt,
            prop_inc_plt,
            nrow = 2, labels = c("A", "B"),
            rel_heights = c(0.5, 1), 
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



