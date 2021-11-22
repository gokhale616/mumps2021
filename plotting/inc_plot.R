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
  



y_max_abs <- max(mumps_case_reports_l$sqrt_total, na.rm = TRUE)

abs_inc_plt <- (
  mumps_case_reports_l %.>%
    drop_na(.) %.>% 
    ggplot(.,
           aes(x = year, y = sqrt_cases)) +
    geom_line(aes(colour = age_cohort), size = 0.8, alpha = 0.95) +
    # geom_area(aes(colour = age_cohort, fill = age_cohort), 
    #           position = position_dodge(), alpha = 0.3, size = 0.9) +
    labs(y = expression(paste(sqrt(Cases), phantom(100))),
         x = "") +
    scale_x_continuous(#expand = c(0.000, 0.000),
                       breaks = gen_x_breaks) +
    scale_y_continuous(#expand = c(0.0, 0),
                       breaks = (y_max_abs*c(0, 0.25, 0.5, 0.75, 1)) %.>% floor(.)) +
    scale_colour_manual(values = age_cols) +
    scale_fill_manual(values = age_cols) +
    project_theme +
    # theme(axis.text.x = element_blank()) +
    cap_axes +
    theme(axis.title = element_text(size = 8.5), 
          plot.margin = unit(rep(mar_val, 4), "cm"), 
          legend.position = "none"
          )
  )


# y_max_inc <- max(mumps_incidence_rate_l$ac_inc_rate, na.rm = TRUE)


anno_data <- tibble(label = c("Initial\ndecline", "First\nresurgence", "Second\nresurgence", "", ""),
                    label_2 = "Sustained\ntransmission",
                    x_2 = 2013,
                    x = c(1981, 1987, 2006, 2009, 2016), 
                    y = 8.5, 
                    yend = 6)


anno_data_2 <- tibble(label = "Sustained\ntransmission",
                      x = 2012, 
                      y = 8.5, 
                      yend = 6)


#"#EF629F"
inc_rate_plt <- (
  mumps_incidence_rate_l %.>% 
    ggplot(., 
           aes(x = year, y = ac_inc_rate)) +
    geom_area(stat = "identity", alpha = 0.9, fill = "#FF4E50", colour = "grey30") +
    geom_segment(data = anno_data, 
                 aes(x = x, xend = x, y = y, yend = yend), 
                 arrow = arrow(length = unit(0.09, "npc")),
                 colour = "grey30") +
    geom_segment(aes(x = 2008.95, xend = 2016.05, y = 8.49, yend = 8.49), colour = "grey30") +
    geom_text(data = anno_data, 
              aes(label = label, x = x, y = 12), colour = "grey30", size = 3) +
    geom_text(data = anno_data_2, 
              aes(label = label, x = x, y = 12), colour = "grey30", size = 3) +
    labs(y = expression(paste(Cases~per~10^5, phantom(1000000))),
         x = "") +
    scale_x_continuous(#expand = c(0.000, 0.000), 
      breaks = gen_x_breaks) + 
    scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0, 15)) +
    project_theme +
    # theme(axis.text.x = element_blank()) +
    cap_axes +
    theme(axis.title = element_text(size = 8.5),
          plot.margin = unit(rep(mar_val, 4), "cm")) 
)


prop_inc_plt <- (
  mumps_case_reports_l %.>% 
  ggplot(., 
         aes(x = year, y = cases, fill = age_cohort)) +
  geom_area(stat = "identity", position = "fill") +
  labs(x = "Year", 
       y = expression(paste(Percent~total, phantom(10000000))),
       fill = "Age\nCohort") +
  scale_fill_manual(values = age_cols) +
  scale_x_continuous(#expand = c(0, 0), 
                     breaks = gen_x_breaks) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent_format(accuracy = 1)) +  
  project_theme +
  cap_axes + 
  theme(axis.title = element_text(size = 8.5), 
        plot.margin = unit(rep(mar_val, 4), "cm")) 
  )


inc_legend <- prop_inc_plt %.>% get_legend(.) 

prop_inc_plt <- prop_inc_plt + theme(legend.position = "none") 

incidence_plt <- (
  plot_grid(inc_rate_plt,
            abs_inc_plt, 
            prop_inc_plt,
            nrow = 3, labels = c("A", "B", "C"),
            # rel_heights = c(1, 1), 
            align = "hv", axis = "b")
  ) 

# final incidence plot
incidence_w_legend_plt <- (
  plot_grid(incidence_plt, inc_legend, 
            rel_heights = c(1, 0.1), nrow = 2)
  ) 


# maps go here
map_mumps <- function(mumps_geog, fill_var = Incidence, 
                      breaks) {
  
  # browser()
  # generate a quosure type variables
  enquo_fill_var <- enquo(fill_var)
  
  
  
  #browser()
  # depends on the preprocess_tycho_data.R
  mumps_geog %.>% 
    ggplot(., 
           aes(x = X, y = Y, group = Group)) +
    geom_polygon(aes(fill = !!enquo_fill_var)) +
    geom_path(size = 0.06, colour = "grey30")+
    labs(x = "", y = "", 
         fill = expression(Cases~per~10^5)) +
    annotate(geom = "text", x = 0.4e6, y = 0.65e6, 
             label = mumps_geog$Year %.>% unique(.), 
             parse = TRUE, size = 3)+
    scale_fill_gradient(low = "#3a7bd5", high = "#FF4E50", na.value = NA,
                        breaks = breaks, 
                        limits = c(0, 10))+
    project_theme +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank(),
          legend.position = c(0.72, 0.1), 
          legend.title = element_text(size=9),
          legend.text = element_text(size=8), 
          panel.grid = element_blank(), 
          plot.margin = unit(rep(mar_val, 4), "cm")) +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal", 
                                 barheight = 0.3)
    )
  
}

#"#F9D423"

set_map_data <- function(year_range) {
  
  # browser()
  
  len_v <- length(year_range)
  
  mumps_demog_geog_annual %.>% 
    filter(., Year %in% year_range) %.>% 
    mutate(., Year = ifelse(len_v == 1, 
                            year_range[1] %.>% as.character(.), 
                            paste0(year_range[1], "-", year_range[len_v])
                            )
           ) %.>% 
    group_by(., Year, State) %.>% 
    mutate(., Incidence = mean(Incidence, na.rm = TRUE)) %.>% 
    ungroup(.)
}


mumps_demog_geog_annual_1985_91 <- (
  set_map_data(year_range = c(1985:1989))  
) 
  
# mumps_demog_geog_annual_1985_89$Incidence %.>% max(., na.rm = TRUE)
# mumps_demog_geog_annual_1985_89$Incidence %.>% min(., na.rm = TRUE)

map_1985_91_plt <-  mumps_demog_geog_annual_1985_91 %.>% 
  map_mumps(., breaks = c(2, 4, 6, 8))



mumps_demog_geog_annual_2006_12 <- (
  set_map_data(year_range = c(2006:2012))  
) 

# mumps_demog_geog_annual_2006_18$Incidence %.>% max(., na.rm = TRUE)
# mumps_demog_geog_annual_2006_18$Incidence %.>% min(., na.rm = TRUE)

map_2006_12_plt <-  mumps_demog_geog_annual_2006_12 %.>% 
  map_mumps(., breaks = c(2, 4, 6, 8))



map_grid_plt <- plot_grid(map_1985_91_plt, map_2006_12_plt, 
                          labels = c("D", "E"))



 # put all of the figure 1
incidence_age_geog <- plot_grid(incidence_w_legend_plt, map_grid_plt, 
                                nrow = 2, rel_heights = c(1, 0.5))




incidence_line_plt <- mumps_case_reports_l %.>% 
  ggplot(., aes(x = year, y = sqrt_cases)) +
  geom_line(size = 0.8) +
  labs(y = expression(sqrt(Cases)), 
       x = "Year") +
  facet_grid(rows = vars(age_cohort), scales = "free_y") +
  project_theme
  






