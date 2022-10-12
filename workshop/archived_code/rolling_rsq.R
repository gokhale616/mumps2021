source("./00/src.R", chdir = TRUE)
test_data <- read_csv("./mumps_data_for_Rsq.csv")





test_data_2 <- (
  test_data %.>% 
    group_by(., age_class) %.>%   
    mutate(., 
           roll_or = slider::slide_dbl(.x = cur_data(), .f = ~calculate_Rsq_temp(.x), 
                                       .before = 5,
                                       .complete = TRUE))
  )
       


test_data_2 %.>% 
  ggplot(.)+
  geom_line(aes(x = year, y = roll_or, colour = "Rolling~R^2~(6~year~window)"))+
  facet_grid(rows= vars(age_class), scales = "free_x")+
  scale_colour_manual(name = "", 
                      values = c("Rolling~R^2~(6~year~window)" = "#1f77b4"), labels = parse_format())+
  project_theme+
  cap_axes()







