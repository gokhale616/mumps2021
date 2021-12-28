source("../00/src.R", chdir = TRUE)
source("../plotting/prep_estm_tables.R", chdir = TRUE)

# get some values ready for this plot 
# mle_values
best_model <- all_result_df %.>% filter(., best_fit_covar == 1) 

# convert tibble to vector of param values for simulations 
best_model_p_vec <- (
  best_model %.>% 
  select(., -c(R0, Rp, hypothesis, vacc_covariate, d_AIC, best_fit_covar)) %.>% 
  mutate(., epsilon2 = 0, p_intro = 6, t_intro = 3000) %.>% 
  unlist(.) %.>% sim_p_vals(.)
  )


# test how the pdf looks for different values 
time_tibble <- tibble(time = seq(0, 100, by = 1/1e3))


prob_exp <- function(time, rate) {
  exp(-rate*time)
}

log_t <- function(prob,rate) {
  -log(prob)/rate
}

time_tibble %<>% mutate(., p_immunity_lost = prob_exp(time, rate = (1/best_model$dwan+1/80)))

# base on the the estimated R0 what is the critical level of vaccination?
critic_lev_vacc <- 1-1/best_model$R0

t_critic_lev_vacc_reached <- log_t(prob = critic_lev_vacc, rate = (1/best_model$dwan+1/80))

anno_segment <- (
  time_tibble %.>% 
    filter(., time %in% c(5, 10, 20, 40, 80)) %.>% 
    transmute(., 
              x_c = time, 
              y_c = p_immunity_lost, zero = 0) %.>% 
    bind_rows(., 
              tibble(x_c = t_critic_lev_vacc_reached, 
                     y_c = critic_lev_vacc, 
                     zero = 0)) %.>% 
    arrange(., x_c)
  )


y_low_lim <- prob_exp(100, rate = (1/best_model$dwan+1/80))

immune_distbn_plot <- (
  time_tibble %.>% 
    ggplot(.) +
    geom_segment(aes(x = 0, xend = time, y = p_immunity_lost, yend = p_immunity_lost, 
                     colour = p_immunity_lost)) +
    geom_point(data = anno_segment, aes(x = x_c, y = y_c, colour = y_c), shape = 21, fill = "white", 
               size = 4) +
    geom_label(data = anno_segment, 
              aes(label = paste0("(", round(x_c, 1), ", " , round(y_c, 2)*100, "%)"), 
                  x = x_c, y = y_c, fill = y_c), 
              nudge_x = 5.5, nudge_y = 0.015, colour = "white") +
    annotate(geom = "segment", y = 0.945, yend = 0.945, x = 20, xend = 15, 
             arrow = arrow(length = unit(2, "mm")), colour = "grey30") +
    annotate(geom = "text", y = 0.945, x = 30, label= "Critical vaccination level", 
             arrow = arrow(length = unit(2, "mm")), colour = "grey30") +
    labs(x = "Time Since Immunization (Years)", 
         y = "Percent Immune\nPost Vaccination") +
    scale_y_continuous(limits = c(y_low_lim, 1), breaks = seq(round(y_low_lim, 1), 1, by = 0.15), 
                       labels = scales::percent) +
    scale_colour_gradient(low = "grey30", high = "#6be585", 
                          limits = c(y_low_lim, 1), breaks = seq(round(y_low_lim, 1), 1, by = 0.15), 
                          guide = FALSE) +
    scale_fill_gradient(low = "grey30", high = "#6be585",
                          limits = c(y_low_lim, 1), breaks = seq(round(y_low_lim, 1), 1, by = 0.15), 
                          guide = FALSE) +
    project_theme +
    cap_axes
  )


# lets see how the next generation matrix looks

eigen_analysis <- calculate_R0_mq(p_vals = c(best_model_p_vec, params_for_R0))


ngm <- eigen_analysis$K_mat[1:5, 6:10] 
dimnames(ngm) <- list(age_names, age_names)

ngm_plot <- (
  ngm %.>% 
    plot_contact_matrix(., col_max = "#ec2F4B",
                        limits = c(0, 15), 
                        breaks = seq(0, 15, 5), 
                        fill_lab = "Average Number\nOf New Cases") +
    project_theme + 
    cap_axes +
    guides(fill = guide_colorbar(frame.colour = "black", 
                                 ticks.colour = "black", 
                                 title.position = "top", 
                                 direction = "horizontal", 
                                 barheight = 0.5))
  )

eigen_vector_tibble <- (
  eigen_analysis$other_output$vectors[1:5, 1] %.>% 
    tibble(.) %.>%
    transmute(., age_class = age_names, reproductive_contri = .)
  )
  
eigen_vector_tibble %.>%   
  ggplot(., aes(x = age_class, y = reproductive_contri)) +
  geom_bar(stat = "identity", fill = "#ec2F4B" , width = 0.45) +
  labs(x = "Age Cohort", 
       y = "Proportion of Cases") +
  project_theme +
  cap_axes










