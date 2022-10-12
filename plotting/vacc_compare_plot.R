# load the cdc estimates 
mmwr_mmr_ests <- (
  read_csv("../raw_data/MMWR_MMR_ests.csv") %.>% 
    select(., Year, `Estimate-CDC`, CI) %.>%  
    transmute(., 
              year = Year, 
              mean = `Estimate-CDC`, 
              ci_lower = mean-CI, 
              ci_upper = mean+CI)
    
    
)

# use one of the co-variate to extraxt recent vaccine covareges
who_mmr_ests <- (
  mod_mumps_covariates_sigmoidal %.>% 
    filter(., year %in% mmwr_mmr_ests$year) %.>% 
    select(., year, p1, p2)
  )



# combine the two data sets

mmr_ests_combined <- (
  mmwr_mmr_ests %.>% 
    right_join(., 
               who_mmr_ests, 
               by = "year") %.>% 
    gather(., key = "dose", value = "est", c(p1, p2)) %.>% 
    mutate(.,
           dose = ifelse(dose == "p1", "Neonatal", "Booster"),
           est = est*100)
)


# comparison plot for the available data 
vacc_ests_compare_plt <- (
  mmr_ests_combined %.>% 
    ggplot(.) +
    geom_linerange(aes(x = est, ymin = ci_lower, ymax = ci_upper, colour = dose))+
    geom_point(aes(x = est, y = mean, fill = year, colour = dose), pch = 21, size = 2.5) +
    geom_smooth(aes(x = est, y = mean, fill = year, colour = dose), 
                method='lm', se = FALSE)+
    labs(x = "WHO reported\nvaccine estimates (%)",
         y = "MMWR reported\nvaccine estimates (%, >=1 Dose)",
         fill = "Year",
         colour = "Dose"
         )+
    project_theme+
    scale_fill_gradient(low = "white", high = "grey30", 
                        limits = c(2000, 2017), breaks = c(2000, 2008, 2017)) +
    scale_colour_manual(values = c("#F16529", "#0575E6"), 
                        breaks = c("Neonatal", "Booster")) +
    cap_axes() +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal"), 
           colour = guide_legend(nrow = 2)
           )
  )
  





