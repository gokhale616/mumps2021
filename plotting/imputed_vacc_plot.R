# extract mean and median of the filtered distribution
imputed_vacc_coverage <- (
  filtered_vacc %.>% 
    right_join(., 
               by = "year",
               mumps_weekly_case_reports %.>% 
                 select(., -cases))  %.>% 
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
           key = "imputation_res", value = "trace", -c(PeriodMidDate, statistic))
)


# plot on imputed vaccine values
imputed_vacc_coverage <- (
  imputed_vacc_coverage_fin %.>%
    ggplot(., aes(x = PeriodMidDate, y = trace)) +
    geom_line(aes(linetype = statistic), size = 1.2) +
    labs(x = "Year",
         y = "Vaccine Coverage",
         linetype = "Filtered\nStatistic")  +
    facet_grid(rows = vars(imputation_res)) +  
    project_theme +
    scale_x_date(breaks = as.Date(c("1968-01-03", "1974-01-01", "1980-01-01", "1985-01-01")), 
                 labels = date_format("%Y")) +
    scale_y_continuous(limits = c(0, 0.86),
                       breaks = c(0, 0.29, 0.57, 0.86), 
                       labels = scales::percent) +
    scale_linetype_manual(values = c(4, 1), 
                          labels = c("Mean", "Median"))+
    cap_axes +
    theme(legend.position = c(0.25, 0.8), 
          text = element_text(size = 20)) +
    guides(colour = guide_legend(direction = "horizontal", 
                                 nrow = 2))
) 


# 


