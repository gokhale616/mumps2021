source("../00/src.R", chdir = TRUE)

# test how the pdf looks for different values 
time_tibble <- tibble(time = seq(0, 100, by = 1/10))


prob_exp <- function(time, rate) {
  exp(-rate*time)
}



time_tibble %<>% mutate(., p_immunity_lost = prob_exp(time, rate = (1/111+1/80)))


anno_segment <- (
  time_tibble %.>% 
    filter(., time %in% c(5, 10, 20, 40, 80)) %.>% 
    transmute(., x_c = time, y_c = p_immunity_lost, zero = 0)
  )


time_tibble %.>% 
  ggplot(., aes(x = time, y = p_immunity_lost)) +
  geom_area(alpha = 0.6, fill = "#ff7e5f") +
  geom_segment(data = anno_segment, aes(x = x_c, xend = x_c, y = 0, yend = y_c), linetype = "dotdash") +
  geom_segment(data = anno_segment, aes(x = 0, xend = x_c, y = y_c, yend = y_c), linetype = "dotdash") +
  geom_point(data = anno_segment, aes(x = x_c, y = y_c), shape = 21, colour = "#4286f4", fill = "white", 
             size = 2) +
  labs(x = "Time since immunization (Years)", 
       y = "Percent Immune") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), labels = scales::percent) +
  project_theme +
  cap_axes

