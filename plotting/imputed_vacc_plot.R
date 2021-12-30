# extract mean and median of the filtered distribution


imputed_vacc_coverage_all <- (
  filtered_vacc %.>% 
    right_join(., 
               by = "year",
               mumps_weekly_case_reports))
    
    
imputed_vacc_coverage <- (
  imputed_vacc_coverage_all %.>% 
    select(., PeriodMidDate, p) %.>% 
    group_by(., PeriodMidDate)  %.>% 
    summarize(., 
              median = median(p), 
              mean = mean(p)) %.>% 
    ungroup(.)  %.>% 
    gather(., "statistic", "trace", -PeriodMidDate)  
)


# defining an annual time of vaccine coverage
imputed_vacc_coverage_fin <- (
  imputed_vacc_coverage %.>%
    mutate(.,
           year = year(PeriodMidDate)) %.>%
    group_by(., year, statistic) %.>%
    mutate(.,
           Weekly = trace,
           Annual = mean(trace)) %.>%
    ungroup(.) %.>% 
    select(., -c(trace, year)) %.>%
    gather(.,
           key = "imputation_res", value = "trace", -c(PeriodMidDate, statistic)) %.>% 
    mutate(., imputation_res = factor(imputation_res, levels = c("Weekly", "Annual")))
)

# anno_case_data
anno_case_data <- (
  imputed_vacc_coverage_all %.>% 
    drop_na(.) %.>% 
    select(., PeriodMidDate, cases) %.>% 
    distinct(.) 
  )


# show the missing data plot for neonatal dose
neonatal_dose_plot <- (
  mumps_covariates %.>% 
    filter(., year>1960) %.>% 
    mutate(., p1 = ifelse(year > 1967 & year < 1986, NA, p1)) %.>% 
    select(., year, p1) %.>% 
    ggplot(., aes(x = year, y = p1)) +
    geom_line(colour = "grey30", size = 1) +
    geom_point(colour = "grey30", shape = 21, size = 2, fill = "white") +  
    annotate(geom = "rect", xmin = 1968, xmax = 1985, ymin = 0, ymax = 1, 
             fill = "deepskyblue1", alpha = 0.4) +
    labs(x = "Year", y = "Neonatal Vaccine Coverage           ")+
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), 
                       labels = scales::percent) + 
    project_theme +
    cap_axes()#+
    #theme(text = element_text(size = 15))
    )

    
  



# plot on imputed vaccine values
imputed_vacc_coverage_plot <- (
  imputed_vacc_coverage_fin %.>%
    ggplot(.) +
    geom_line(aes(x = PeriodMidDate, 
                  y = trace, 
                  colour = imputation_res,
                  linetype = statistic), size = 1.2) +
    labs(x = "Year",
         y = "Missing Vaccine Coverage",
         linetype = "Filtered\nStatistic")  +
    facet_grid(rows = vars(imputation_res)) +  
    project_theme +
    scale_x_date(breaks = as.Date(c("1968-01-03", "1974-01-01", "1980-01-01", "1985-01-01")), 
                 labels = date_format("%Y")) +
    scale_y_continuous(limits = c(0, 0.86),
                       breaks = c(0, 0.29, 0.57, 0.86), 
                       labels = scales::percent) +
    scale_linetype_manual(values = c(4, 1), 
                          labels = c("Mean", "Median")) +
    scale_colour_manual(values = c("Weekly" = "grey30", 
                                   "Annual" = "deepskyblue1")) +
    cap_axes() +
    theme(legend.position = c(0.25, 0.8)#, 
          #text = element_text(size = 15)
          ) +
    guides(linetype = guide_legend(direction = "horizontal", 
                                 nrow = 2), 
           colour = "none")
) 



case_data_for_imptn_plot <- (
  anno_case_data %.>% 
    ggplot(.) + 
    geom_area(aes(x = PeriodMidDate, y = cases), fill = "grey30") +
    labs(y = "Weekly Cases            ", x = "") +
    scale_x_date(breaks = as.Date(c("1968-01-03", "1974-01-01", "1980-01-01", "1985-01-01")), 
                 labels = date_format("%Y")) +
    scale_y_continuous(labels = scales::scientific) +
    project_theme+
    cap_axes() #+
    #theme(text = element_text(size = 15))
    )



# 
imputed_vacc_coverage_plt <- (
  plot_grid(case_data_for_imptn_plot, imputed_vacc_coverage_plot, nrow = 2, rel_heights = c(0.5,1), 
          align = "hv", axis = "lr", labels = c("B", "C"))
  )



imputed_vacc_coverage_grid_plt <- (
  plot_grid(neonatal_dose_plot, imputed_vacc_coverage_plt, nrow = 2, rel_heights = c(0.5,1), 
            align = "h", axis = "bl", labels = c("A", ""))
)
