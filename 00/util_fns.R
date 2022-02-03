# group of functions useful for hacking the colour bar if you require to add on a target
ncolours <- seq(0, 1, length.out = 100)

modify_palette_by_target <- function(colours = seq_gradient_pal("#dd3e54", "#6be585")(ncolours), 
                                     target, 
                                     plot_var, 
                                     replace_colour = "#642B73") {
  
  # browser()
  # rank data to lineraize values
  if(target %in% plot_var) {
    message("target exists!!")
    int_target = target
  } else {
    # browser()
    message("value closest to the target used")
    int_target = plot_var[which.min(abs(target - plot_var))]
  }
  
  
  ranked_plot_var <- rank(plot_var)
  
  ranked_target <- ranked_plot_var[which(plot_var == int_target)]
  
  # normalize the target interval 
  ranked_range <- range(ranked_plot_var) 
  
  target_interval <- ranked_plot_var[which(plot_var == int_target)]*c(0.99, 1.01)
  
  normalized_target_interval <- (target_interval - ranked_range[1])/diff(ranked_range)
  
  # browser()
  
  ramp <- scales::colour_ramp(colours)
  
  function(x) {
    # Decide what values to replace
    # browser()
    ranked_x <- rank(x)
    
    normalized_ranked_x <- (ranked_x - range(ranked_x)[1])/diff(range(ranked_x))
    
    replace <- normalized_ranked_x > normalized_target_interval[1] & 
      normalized_ranked_x < normalized_target_interval[2]
    
    out <- ramp(x)
    # Actually replace values
    out[replace] <- replace_colour
    out
  }
}




# a wraper function to recursively bind rows of a list 
form_dfr <- function(list){do.call(rbind, lapply(1:length(list), function(x) {list[[x]]} ))}

# this 'convenience' function modifies vector of default params to set hypothesis specific param estimates
sim_p_vals <- function(estm_vect, default_p_vals = param_vals_est) {
  
  replace_these <- names(default_p_vals)[names(default_p_vals) %in% names(estm_vect)] 
  
  tmp_p_vals <- default_p_vals
  tmp_p_vals[replace_these] <- estm_vect[replace_these]
  
  tmp_p_vals
  
}


# a wrapper to unlist and unname a list object.
un_list_name <- function(x) {
  x %.>% 
    unlist(.) %.>% 
    unname(.)
}  


# operator for easy analysis
`%nin%` <- Negate(`%in%`)

# this function calculates quntiles for a grouped data frame 
calculate_quantile <- function(grouped_data, var, p = c(0.025, 0.5, 0.975)) {
  # browser()
  enquo_var <- enquo(var)
  
  summarise(., 
            qs = quantile(!!enquo_var, p, na.rm = TRUE), 
            prob = as.character(p), 
            .groups = 'drop') %.>%
    spread(., 
           key = prob, value = qs) %.>% 
    ungroup(.)
  
}


# this is a convenience function to gather model compartments
gather_comp_1_sim <- function(sim_data) {
  sim_data %.>% 
    select(., -c(`.id`, C)) %.>% 
    gather(., key = "Comp", value = "value", -year) 
}


# function to calculate the AICc for a given loglikilhood
calculate_aic <- function(loglik, npar) {
  return(2*npar - 2*loglik)
}


# function: produces n replicates to implement R measure 
param_replicate_matrix <- function(param_v, n = 1000) {
  t(t(param_v)) %*% rep(1, n)
}   

# function: converts a factor to numeric useful for preserving the Year vector 
as_numeric_factor <- function(x) {as.numeric(levels(x))[x]}

# function: simulates from the observation model specified in pomp
sim_obs_model <- function(po_obj, params, times, 
                          nsim, method = "ode23", 
                          summarise_obs_process = TRUE, 
                          root_transform = TRUE, 
                          quantile_prob = c(0.025, 0.5, 0.975), 
                          quantile_prob_names = c("0.025", "0.5", "0.975")) {
    
    # browser()
    # apply the observation process and generate simulations   
    # solve the ode system
    soln_array <- po_obj %.>% 
      trajectory(., params = params, method = method, format = "array") 
    
    # simulate the obs process - this is an 3D array 
    obs_soln_array <- po_obj %.>%
      rmeasure(., x = soln_array, times = times, 
               params = param_replicate_matrix(params, n = nsim)) 
    
    # dress the result to produce a neat data frame that has the simulated obs process
    result_df <- do.call(rbind, lapply(1:nsim, function(x) {t(obs_soln_array[,x,])})) %.>% 
      as_tibble(.) %.>% 
      transmute(., 
                D_1 = D_1, D_2 = D_2, D_3 = D_3, 
                D_4 = D_4, D_5 = D_5, D_u = D_u,
                year = rep(times, times = nsim), 
                `.id` = rep(1:nsim, each = length(times))
                ) %.>% 
      select(., year, `.id`, starts_with("D_"))
  
  
    # feed the names of the age classes
    change_these_colnames <- which(colnames(result_df) %in% c(sprintf("D_%d", 1:5), "D_u"))
    
    colnames(result_df)[change_these_colnames] <- age_names_u
  
    # preset in a long form
    result_df_l <- result_df %.>% 
      gather(., 
             key = "age_class", 
             value = "cases", -c(year, `.id`), factor_key = TRUE) 
      
    
    if(root_transform == TRUE) {
      
      message("simulations are root-transformed!!")
      result_df_l_trans <- (
        result_df_l %.>% 
          mutate(., cases = sqrt(cases))
      )
    
    } else {
      
        message("simulations are in natural scale!!")
        result_df_l_trans <- result_df_l
    }
    
    
    
    if(summarise_obs_process == TRUE) {
      
      fin_result <- (
        result_df_l_trans %.>%
        select(., -`.id`) %.>% 
        group_by(., year, age_class) %.>%
        summarise(., 
                  qs = quantile(cases, quantile_prob, na.rm = TRUE), 
                  prob = quantile_prob_names, 
                  .groups = 'drop') %.>%
        spread(., key = prob, value = qs) %.>% 
        ungroup(.)
        )
        
      
    } else {
      
        fin_result <- result_df_l
      
    }
      
    fin_result
  
}
 


# my if_else for NA to rplace na values for various parameter columns 

na_if_else <- function(x, replace_na_with) {if_else(is.na(x) == TRUE, replace_na_with, x)}


################################ Reporductive numbers go here ################################################

# functions to transfrom back and forth between q and R0 

calculate_q <- function(R0) {
  return(R0/45.39024)
}

calculate_R0 <- function(q) {
  return(q*45.39024)
}


# plot the contact matrix
plot_contact_matrix <- function(contact_matrix = contact_sym_Lm10, 
                                plt_title = NULL, plt_subtitle = NULL, 
                                col_min = "#ffd89b", col_max = "#480048", limits = c(0, 8), 
                                breaks = seq(0, 8, by = 2), fill_lab = "Mean Daily\nContacts") {
  
  sapply(c("reshape2", "ggthemes"), library, character.only = TRUE)
  
  ggplot(data = melt(contact_matrix), mapping = aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_continuous(low = col_min, high = col_max, limits = limits,
                          breaks =  breaks) +
    coord_flip() + labs(title = plt_title, subtitle = plt_subtitle,  
                        x = "Contact Age", y = "Reporter Age", fill = fill_lab) + 
    theme(aspect.ratio = 1, 
          legend.position="bottom",
          axis.text.x = element_text(angle = 90)) -> contact_matrix_plt  
  
  return(contact_matrix_plt) 
  
}

############################################################################################################## 
#################################### functions for preprocessing ############################################# 
######################################## demography data #####################################################
##############################################################################################################
# break the migration term into two values 
break_mu_into_2 <- function(x, c = 0.1) {
  #browser()
  abs_x <- abs(x)
  
  mu_into_2 <- case_when(
    x > 0 ~ sort(c((abs_x+c), -c)),
    x < 0 ~ sort(c(-(abs_x+c), c)),
    x == 0 ~ c(0, 0))
  
  names(mu_into_2) <- c("efflux", "influx")
  
  as_tibble(t(mu_into_2))
  
}

#break_mu_into_2(2)


# function: to make the spline for the mu
mu_splinefun <- function(t) {
  return(nu_splinefun(t) - (change_pop_splinefun(t)/total_pop_splinefun(t)))
}

# function: estimator coefficient of variation
cVar <- function(x) {return((sd(x)*100/mean(x)))}


# function: to interpolate the vaccine uptake
logistic <- function(t, t0 = 1970, a = 0.86, b = 1, 
                     c = 1, d = 2, 
                     ts = 1968, te = 1979) {
  
  if(ts > te) {
    stop("ts < te is not satified")
  }
  
  case_when(t < ts ~ 0,
            t > (ts-1) & t < (te+1) ~ a*((b + c*exp(-d*(t-t0)))^-1),
            t > te ~ a) 
  
  
}

# function: to treat the covarariates for the modified objective function 

# log_par_defaults <- c(t0_p1 = 1972, b_p1 = 1, c_p1 = 1, d_p1 = 0.5, 
#                       t0_p2 = 1994, b_p2 = 1, c_p2 = 1, d_p2 = 1) 

treat_covar_data <- function(par = log_par_defaults) {
  # browser()
  # Extract p1 
  p1 <- mumps_covariates %.>% 
    select(., p1, year) %.>% 
    filter(., !(year > 1967 & year < 1980)) 
  # Treat p1
  p1_log_in <- tibble(year = seq(1968, 1979)) %>% 
    mutate(p1 = logistic(t = year, t0 = par["t0_p1"],
                         a = 0.86, b = par["b_p1"], 
                         c = par["c_p1"], d = par["d_p1"],
                         ts = 1968, te = 1979))
  # Reattach p1
  p1 %>% 
    bind_rows(p1_log_in) %>% 
    arrange(year) -> p1_in
  
  
  # Extract p2
  mumps_covariates %>% 
    select(p2, year) %>% 
    filter(!(year > 1988 & year < 2000)) ->  p2
  # Treat p2
  p2_log_in <- tibble(year = seq(1989, 1999)) %>% 
    mutate(p2 = logistic(t = year, t0 = par["t0_p2"],
                         a = 0.85, b = par["b_p2"], 
                         c = par["c_p2"], d = par["d_p2"],  
                         ts = 1989, te = 1999))
  # Reattach p2
  p2 %>% 
    bind_rows(p2_log_in) %>% 
    arrange(year) -> p2_in
  
  # get rid of right constant interpolation of p1 an p2 
  mumps_covariates %>% 
    select(-starts_with("p")) -> mumps_demog_data
  
  # Reattach p1 and p2 
  p1_in %>%
    full_join(p2_in, by = "year") %>% 
    full_join(mumps_demog_data, by = "year") -> treated_mumps_covariates  
  
  treated_mumps_covariates
  
  
}


# function: plot all the co-variates 
plot_covars <- function(covar_data, filter_from = 1950) {
  # browser()
  # levels of the age specific (as) data 
  cov_levels_as <- c("Population", "Migration (Year^-1)")
  
  # levels of the age independent (ai) data 
  cov_levels_ai <- c("Births", "Age Record Probabilty", "Neonatal Dose", "Booster Dose")
  
  covar_data %>% nrow() -> cov_nrow
  
  
  # generating 
  
  anno_data1 <- data.frame(Covariate = cov_levels_as, 
                           y = c(0.3e8, -4.7e-2), 
                           x = 2015, 
                           labs = c("bold(e)", "bold(f)"))
  
  
  anno_data2 <- data.frame(Covariate = cov_levels_ai, 
                           ycord = c(0.2, 0.50, 0.2, 0.2), 
                           lab = c("","bold(Case~Data)","",""), 
                           y = c(0.15, 0.15, 0.15, 0.15), 
                           x = 2015, 
                           labs = c("bold(a)", "bold(b)", "bold(c)", "bold(d)")) 
  
  
  
  # plot age specific covariates
  covar_data %>% 
    select(Year, starts_with("N_"), starts_with("MU_")) %>% 
    gather(key = "AgeCohort", value = "Value", -Year) %>% 
    mutate(AgeCohort = rep(rep(age_names, each = cov_nrow), times = length(cov_levels_as)), 
           Covariate = factor(rep(cov_levels_as, each = cov_nrow*length(age_names)), 
                              levels = cov_levels_as)) %>% 
    filter(Year > filter_from) %>% 
    ggplot(aes(x = Year, y = Value)) +
    geom_area(aes(fill = AgeCohort)) +
    annotate(geom = "rect", xmin = 1977, xmax = Inf, 
             ymin = -Inf, ymax = Inf, 
             fill = "#B2022F", colour = "#B2022F", 
             linetype = "dotdash", alpha = 0.2) +
    geom_text(data = anno_data1, aes(x = x, y = y, label = labs), 
              parse = TRUE) +
    facet_wrap(.~Covariate, nrow = 1, scales = "free") +
    labs(y = "", 
         fill = "Age\nCohort") +
    scale_fill_brewer(palette = "Greens", direction = -1) +
    scale_x_continuous(expand = c(0,0.005)) +
    scale_y_continuous(expand = c(0,0.005), labels = scientific) +
    project_theme + 
    cap_axes +
    theme(strip.text.x = element_blank()) +
    guides(linetype = "none", 
           fill=guide_legend(nrow=2)) -> age_stratified_cov_plt
  
  age_class_legend <- get_legend(age_stratified_cov_plt)
  
  age_stratified_cov_plt <- age_stratified_cov_plt +
    theme(legend.position = "none")
  
  # browser()
  # plot age independent covariates
  covar_data %>% 
    select(Year, Births, p1, p2, eta_a) %>%
    mutate(Births = Births/max(Births)) %>% 
    gather(key = "Covariate", value = "Value", -Year) %>% 
    mutate(Covariate = factor(case_when(Covariate == "Births" ~ cov_levels_ai[1],
                                        Covariate == "eta_a" ~ cov_levels_ai[2],
                                        Covariate == "p1" ~ cov_levels_ai[3],
                                        Covariate == "p2" ~ cov_levels_ai[4]), 
                              levels = cov_levels_ai), 
           cond_colour_p = case_when(Covariate == cov_levels_ai[3] & Year > 1967 & Year < 1985 ~ TRUE,  
                                     Covariate == cov_levels_ai[4] & Year > 1988 & Year < 2000 ~ TRUE, 
                                     TRUE ~ FALSE), 
           cond_colour_l = case_when(Covariate == cov_levels_ai[3] & Year > 1966 & Year < 1985 ~ TRUE,  
                                     Covariate == cov_levels_ai[4] & Year > 1987 & Year < 2000 ~ TRUE, 
                                     TRUE ~ FALSE)) %>%
    filter(Year > filter_from) %>%
    ggplot(aes(x = Year, y = Value)) +
    geom_line(aes(colour = cond_colour_l, group = 1), size = 0.8) +
    geom_point(aes(colour = cond_colour_p), pch = 21, fill = "white", size = 2) +
    annotate(geom = "rect", xmin = 1977, xmax = Inf, 
             ymin = -Inf, ymax = Inf, fill = "#B2022F", 
             colour = "#B2022F", linetype = "dotdash", alpha = 0.2) +
    geom_text(data = anno_data2, aes(x = 1974.5, y = ycord, label = lab, angle = 90), 
              inherit.aes = FALSE, parse = TRUE) +
    geom_text(data = anno_data2, aes(x = x, y = y, label = labs), 
              inherit.aes = FALSE, parse = TRUE) +
    labs(y = "", 
         x = "") +
    facet_wrap(.~Covariate, nrow = 2, scales = "free") +
    scale_colour_manual(name = "Vaccine\nCoverage", labels = c("Observed","Interpolated"), 
                        values = c("grey30", "#FF4B2B")) +
    scale_x_continuous() +
    scale_y_continuous(limits = c(0, 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = scales::percent
                       ) +
    project_theme +
    cap_axes +
    theme(legend.spacing.y = unit(0., "cm"),
          panel.spacing.y = unit(1.75, "lines"), 
          strip.text.x = element_blank()) -> age_independent_cov_plt
  
  plot_grid(age_independent_cov_plt, age_stratified_cov_plt, 
            align = "hv", axis = "lr", nrow = 2, 
            rel_heights = c(1, 0.43)) -> cov_plt  
  
  # add the legend for age specific plots 
  
  plot_grid(cov_plt, age_class_legend, 
            nrow = 2, rel_heights = c(1, 0.09)) -> cov_plt_with_legend
  
  plot(cov_plt_with_legend)
}

# plot_covars(covar = mumps_covariates)

plot_sweep_dynamics <- function(data, 
                                filter_from = 1976,
                                only_cases = TRUE, 
                                colour_lab,
                                y_lab = NULL,
                                sim_cases = FALSE, 
                                ...) {
  # browser()
  if(sim_cases == TRUE) {
    
    plt <- (
      data %.>% 
      ggplot(., aes(x = year, y = `0.5`, group = Value, colour = Value)) +
      geom_line(aes(colour = Value), size = 0.8) +
      geom_linerange(aes(ymin=`0.025`, ymax=`0.975`))+  
      labs(colour = colour_lab, 
           x = "Year", 
           y = y_lab) +
      facet_grid(rows = vars(age_class), scales="free_y") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(labels = scales::scientific_format(digits = 2))+
      scale_colour_continuous(low = "#2C5364", high = "#FF416C", breaks = c(50, 150, 250, 350)) +
      project_theme +
      guides(colour = guide_colorbar(frame.colour = "black", 
                                     ticks.colour = "black", 
                                     title.position = "top")
             )
      )
  
  } else {
    
      if(only_cases == TRUE) {
        data_int <- data %.>% 
          select(., year, starts_with("C_"), Value)
      } else {
        data_int <- data 
      }
      
      colnames(data_int)[which(colnames(data_int) %in% sprintf("C_%d", 1:5))] <- age_names
      
      # treat data_int for plotting 
      plt_data <- data_int %.>%
        filter(., 
               year > filter_from) %.>% 
        gather(., 
               key = "compartment", value = "cases", 
               -c(year, Value), factor_key = TRUE) 
      
      
      
      # browser()
      plt <- (
        plt_data %.>% 
          ggplot(., aes(x = year, y = cases)) +
          geom_line(aes(colour = Value, group = Value), size = 0.8) +
          labs(colour = colour_lab, 
               x = "Year", 
               y = y_lab) +
          facet_grid(rows = vars(compartment), scales="free_y") +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(labels = scales::scientific_format(digits = 2))+
          scale_colour_continuous(low = "#2C5364", high = "#FF416C", breaks = c(50, 150, 250, 350)) +
          project_theme +
          guides(colour = guide_colorbar(frame.colour = "black", 
                                         ticks.colour = "black", 
                                         title.position = "top")
          )
        
      )  
  }
  
  
  
  
  plt
  
}





















