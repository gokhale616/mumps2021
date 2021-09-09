# operator for easy analysis
`%nin%` <- Negate(`%in%`)

# function to calculate the AICc for a given loglikilhood
calculate_aic <- function(loglik, npar) {
  return(2*npar - 2*loglik)
}


# function: calculate quantiles for sample of a distribution 
give_quantile_list <- function(quantiles = c(0.025, 0.975), denominator = 100) {
  # browser()
  purrr::map(quantiles, ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = map_chr(quantiles, ~paste0(.x*denominator, "%")))
}

# function: produces n replicates to implement R measure 
param_replicate_matrix <- function(param_v, n = 1000) {
  t(t(param_v)) %*% rep(1, n)
}   

# function: converts a factor to numeric useful for preserving the Year vector 
as_numeric_factor <- function(x) {as.numeric(levels(x))[x]}

# function: simulates from the observation model specified in pomp
sim_obs_model <- function(po_obj, params, times, 
                          nsim, method = "ode45", 
                          long_form = FALSE, soln_only = FALSE) {
  
  # produce solution to the sys. of ode
  if(soln_only == FALSE) {
    
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
      as.data.frame(.) %.>% 
      transmute(., 
                D_1 = D_1, D_2 = D_2, D_3 = D_3, 
                D_4 = D_4, D_5 = D_5, Du = Du,
                Year = rep(times, times = nsim), 
                `.id` = rep(1:nsim, each = length(times))
                ) %.>% 
      select(., Year, `.id`, starts_with("D")) %.>% 
      mutate_at(., 
                .vars = vars(starts_with("D")), 
                .funs = function(x) {ifelse(.$Year == min(.$Year), NA, x)}
                ) 
  } else {
      
      # only solve the ODE system
      result_df <- po_obj %.>% 
        trajectory(., params = params, method = method, format = "data.frame")
  
  }
    
  
  # if the long form is required for plotting purposes
  if(long_form == TRUE) {
    
    fin_result <- result_df %.>% 
      gather(., key = "AgeClass", value = "SimObs", -c(Year, `.id`)) 
    
  } else { 
    
    fin_result <- result_df
  
  }
  
  
    fin_result
  
}
 

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

# convert  the list of profile likelihood results into a dataframe
list_to_tibble <- function(counter = 1, res_obj,
                           model_name = model_nms,
                           give_everything_GA = FALSE) {
  
  # browser()
  est_object <- res_obj[[counter]]
  model_nm <- model_name[counter]
  
  
  if(give_everything_GA == TRUE) {
    return(est_object)
    
  } else {
    # browser()
    est_object$GAobj@solution[1,] %.>%
      t(.) %.>% 
      as_tibble(.) %.>% 
      mutate(., 
             logLik =  est_object$GAobj@fitnessValue[1], 
             npar = length(est_object$GAobj@solution[1,]), 
             model = model_nm) -> param_Est     
    
    return(param_Est)
  }
  
}  

# my if_else for NA to rplace na values for various parameter columns 

na_if_else <- function(x, replace_na_with) {if_else(is.na(x) == TRUE, replace_na_with, x)}

# find where two curves intersect - used to calculate confidence intervals for the likelihood surfaces
# the function is sourced from andrewheiss's package reconPlots

curve_intersect <- function(curve1, curve2) {
    
    # browser()
    
    # Approximate the functional form of both curves
    curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
    curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)
    
    # Calculate the intersection of curve 1 and curve 2 along the x-axis
    point_x <- uniroot.all(function(x) curve1_f(x) - curve2_f(x), 
                           c(min(curve1$x), max(curve1$x)))
    
    # Find where point_x is in curve 2
    point_y <- c(min(curve2_f(point_x)), max(curve2_f(point_x)))
  
  return(list(x = point_x, y = point_y))
}

# utility functions around the contact matrix - Know that this bit of analysis is going to useful for quickly 
# accessing the next genertaion matrix for use contact matrix, calculate R0 based on the NextGen method and 
# so on, but very basic...

# Dirac delta function # still working on this function --  
delta <- function(i, j) return(ifelse(i == j, 1, 0)) 

calculate_R0_mq <- function(params, no_vaccine = TRUE) {
  # R effective calculations are still pending  
  N = params$pop                      # vector of total population
  
  q = params$q_mult                   # vector of susceptibility probabilities
  sigma = params$sigma                # rate of becoming infectious 
  gamma = params$gamma                # recovery rate
  
  if(no_vaccine  == TRUE) {
    epsilon = 0              # leakiness  
    alpha = 0                # primary vaccine failure
    delta = 0                # rate of waning failure  
  } else{
      epsilon = params$epsilon            # leakiness  
      alpha = params$alpha                # primary vaccine failure
      delta = params$delta                # rate of waning failure
  }
  
  C = params$contact_matrix           # matrix of contacts
  mu = params$age_durations           # vector of aging and mortality 
  nu = params$nu                      # average birth rate
  # browser()
  # calculate the necessary intermediate compartments   
  # define age specific population sizes   
  Ns = nu*N*mu^-1  
  # Si 
  # Vi
  
  # Define emplty matrices of F and V to calculate the the NextGen Matrix  
  NA_mat = matrix(NA, nrow = nrow(C), ncol = ncol(C)) 
  zero_mat = matrix(0, nrow = nrow(C), ncol = ncol(C))
  
  
  # Intialize 
  F_mat_nz = NA_mat
  
  if(no_vaccine == TRUE) {
    for(i in seq_along(q)) {
      for(j in seq_along(q)) {
        F_mat_nz[i,j] = q[i]*C[i,j]*Ns[i]/Ns[j]
      }
    }  
  } else {
      for(i in seq_along(q)) {
        for(j in seq_along(q)) {
          F_mat_nz[i,j] = q[i]*C[i,j]*(S[i] + epsilon*V[i])/N[j]}
    }
  }
  
  # Calculate the non-zero block of the F_mat matrix
  # Define the F matrix
  F_mat = rbind(cbind(zero_mat, F_mat_nz), cbind(zero_mat, zero_mat))
  
  # Initialize the block matrices of the V matrix, (explain what hey are later) 
  V1 = NA_mat 
  V2 = zero_mat
  V3 = NA_mat
  V4 = NA_mat
  
  for(i in seq_along(q)) {
    for(j in seq_along(q)) {
      if(i == 1) {
        V1[i,j] =  (mu[i] + sigma)*delta(i,j)
        V3[i,j] = -sigma*delta(i,j)
        V4[i,j] =  (mu[i] + gamma)*delta(i,j)  
      } else {
          V1[i,j] =  (mu[i] - mu[i-1] + sigma)*delta(i,j)
          V3[i,j] = -sigma*delta(i,j)
          V4[i,j] =  (mu[i] - mu[i-1] + gamma)*delta(i,j)  
      }
    }
  }
  
  V_mat = rbind(cbind(V1, V2), cbind(V3, V4))
  # browser()
  # calculate the next generation matrix
  K_mat = F_mat%*%solve(V_mat)
  
  R0 = max(eigen(K_mat)$values)
  
  other_output = eigen(K_mat) 
  
  final_output <- list(F_mat = F_mat, V_mat = V_mat, K_mat = K_mat, 
                       R0 = R0, other_output = other_output)
  
  return(final_output)
  
}

# utility function to calculate multiple R0s based on a data of q parameters 

give_R0_estimate <- function(counter, q = q_data) {
  
  q_vec <- unlist(q[counter,])
  
  params <- list(pop = 3e8, nu = 1/80, 
                 age_durations = age_class_duration^-1, 
                 gamma = param_vals["gamma"],
                 sigma = param_vals["sigma"], 
                 q_mult = q_vec,
                 contact_matrix = C)
  
  
  calculate_R0_mq(params = params)$R0
  
}


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
                                col_min = "#ffd89b", col_max = "#480048") {
  
  sapply(c("reshape2", "ggthemes"), library, character.only = TRUE)
  
  ggplot(data = melt(contact_matrix), mapping = aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_continuous(low = col_min, high = col_max, limits = c(0, max(contact_matrix))) +
    coord_flip() + labs(title = plt_title, subtitle = plt_subtitle,  
                        x = "Contact age", y = "Reporter age", fill = "Mean Annual\nContacts") + 
    theme(aspect.ratio = 1, 
          legend.position="bottom",
          axis.text.x = element_text(angle = 90)) -> contact_matrix_plt  
  
  return(contact_matrix_plt) 
  
}

############################################################################################################## 
#################################### functions for preprocessing ############################################# 
######################################## demography data #####################################################
##############################################################################################################

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
  
  # Extract p1 
  mumps_covar_data %>% 
    select(p1, Year) %>% 
    filter(!(Year > 1967 & Year < 1980)) ->  p1 
  # Treat p1
  p1_log_in <- tibble(Year = seq(1968, 1979)) %>% 
    mutate(p1 = logistic(t = Year, t0 = par["t0_p1"],
                         a = 0.86, b = par["b_p1"], 
                         c = par["c_p1"], d = par["d_p1"],
                         ts = 1968, te = 1979))
  # Reattach p1
  p1 %>% 
    bind_rows(p1_log_in) %>% 
    arrange(Year) -> p1_in
  
  
  # Extract p2
  mumps_covar_data %>% 
    select(p2, Year) %>% 
    filter(!(Year > 1988 & Year < 2000)) ->  p2
  # Treat p2
  p2_log_in <- tibble(Year = seq(1989, 1999)) %>% 
    mutate(p2 = logistic(t = Year, t0 = par["t0_p2"],
                         a = 0.85, b = par["b_p2"], 
                         c = par["c_p2"], d = par["d_p2"],  
                         ts = 1989, te = 1999))
  # Reattach p2
  p2 %>% 
    bind_rows(p2_log_in) %>% 
    arrange(Year) -> p2_in
  
  # get rid of right constant interpolation of p1 an p2 
  mumps_covar_data %>% 
    select(-starts_with("p")) -> mumps_demog_data
  
  # Reattach p1 and p2 
  p1_in %>%
    full_join(p2_in, by = "Year") %>% 
    full_join(mumps_demog_data, by = "Year") -> treated_mumps_covar_data  
  
  return(treated_mumps_covar_data)
  
  
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
                           ycord = c(3.2e6, 0.65, 0.1, 0.1), 
                           lab = c("","bold(Case~Data)","",""), 
                           y = c(3250000, 0.38, 0.08, 0.08), 
                           x = 2016, 
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
    theme(strip.text.x = element_blank()) +
    guides(linetype = "none", 
           fill=guide_legend(nrow=2)) -> age_stratified_cov_plt
  
  age_class_legend <- get_legend(age_stratified_cov_plt)
  
  age_stratified_cov_plt <- age_stratified_cov_plt +
    theme(legend.position = "none")
  
  # plot age independent covariates
  covar_data %>% 
    select(Year, Births, p1, p2, eta_a) %>% 
    gather(key = "Covariate", value = "Value", -Year) %>% 
    mutate(Covariate = factor(case_when(Covariate == "Births" ~ cov_levels_ai[1],
                                        Covariate == "eta_a" ~ cov_levels_ai[2],
                                        Covariate == "p1" ~ cov_levels_ai[3],
                                        Covariate == "p2" ~ cov_levels_ai[4]), 
                              levels = cov_levels_ai), 
           cond_colour_p = case_when(Covariate == cov_levels_ai[3] & Year > 1967 & Year < 1980 ~ TRUE,  
                                     Covariate == cov_levels_ai[4] & Year > 1988 & Year < 2000 ~ TRUE, 
                                     TRUE ~ FALSE), 
           cond_colour_l = case_when(Covariate == cov_levels_ai[3] & Year > 1966 & Year < 1980 ~ TRUE,  
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
                        values = c("grey30", "#FFD92F")) +
    scale_x_continuous(expand = c(0,0.005)) +
    scale_y_continuous(expand = c(0,0.005)) +
    project_theme +
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

# plot_covars()

# workshop related functions 
plot_dynamics <- function(data, filter_from = 1976, dynamics, ...) {
  
  if(dynamics == "stoc") {
    
    comp_levels <- c(sprintf("S_%d",1:5),  
                     sprintf("E1_%d",1:5), sprintf("I1_%d",1:5),
                     sprintf("E2_%d",1:5), sprintf("I2_%d",1:5),
                     sprintf("R_%d",1:5), sprintf("V_%d",1:5), 
                     sprintf("A_%d",1:5), sprintf("G_%d",1:5), 
                     sprintf("D_%d",1:5), "Du",
                     sprintf("Nsim_%d",1:5), 
                     "Nsim", "Rcal","Ncal")
  } else if (dynamics == "det") {
  
      comp_levels <- c(sprintf("S_%d",1:5),  
                       sprintf("E1_%d",1:5), sprintf("I1_%d",1:5),
                       sprintf("E2_%d",1:5), sprintf("I2_%d",1:5),
                       sprintf("R_%d",1:5), sprintf("V_%d",1:5), 
                       sprintf("A_%d",1:5), sprintf("G_%d",1:5), 
                       sprintf("Nsim_%d",1:5), 
                       "Nsim", "Rcal","Ncal")
  } else {
      
      stop("'plot' argument does not match stoc/det")
  }
  
  
  data %>%
    filter(Year > filter_from) %>% 
    mutate(Nsim = S_1 + S_2 + S_3 + S_4 + S_5 + 
             E1_1 + E1_2 + E1_3 + E1_4 + E1_5 +
             I1_1 + I1_2 + I1_3 + I1_4 + I1_5 +
             E2_1 + E2_2 + E2_3 + E2_4 + E2_5 +
             I2_1 + I2_2 + I2_3 + I2_4 + I2_5 +
             V_1 + V_2 + V_3 + V_4 + V_5 + 
             R_1 + R_2 + R_3 + R_4 + R_5,
           Nsim_1 = S_1 + E1_1 + I1_1 + E2_1 + I2_1 + R_1 + V_1,
           Nsim_2 = S_1 + E1_1 + I1_1 + E2_1 + I2_1 + R_1 + V_1,
           Nsim_3 = S_1 + E1_1 + I1_1 + E2_1 + I2_1 + R_1 + V_1,
           Nsim_4 = S_1 + E1_1 + I1_1 + E2_1 + I2_1 + R_1 + V_1,
           Nsim_5 = S_1 + E1_1 + I1_1 + E2_1 + I2_1 + R_1 + V_1, 
           Rcal = 3.3e8 - (S_1 + S_2 + S_3 + S_4 + S_5 + 
                             E1_1 + E1_2 + E1_3 + E1_4 + E1_5 +
                             I1_1 + I1_2 + I1_3 + I1_4 + I1_5 +
                             E2_1 + E2_2 + E2_3 + E2_4 + E2_5 +
                             I2_1 + I2_2 + I2_3 + I2_4 + I2_5 +
                             V_1 + V_2 + V_3 + V_4 + V_5),
           Ncal = S_1 + S_2 + S_3 + S_4 + S_5 + 
                  E1_1 + E1_2 + E1_3 + E1_4 + E1_5 +
                  I1_1 + I1_2 + I1_3 + I1_4 + I1_5 +
                  E2_1 + E2_2 + E2_3 + E2_4 + E2_5 +
                  I2_1 + I2_2 + I2_3 + I2_4 + I2_5 +
                  V_1 + V_2 + V_3 + V_4 + V_5 + 
             Rcal) %>% 
    gather(key = "Compartment", value = "Cases", -c(".id", "Year")) %>%
    ggplot(aes(x = Year, group = .id, color = .id, y = Cases))+
    geom_line()+
    facet_wrap(~factor(Compartment, levels = comp_levels), scales="free_y", ncol = 5) +
    scale_y_continuous(labels = scales::scientific_format(digits = 2))+
    theme(...)
  
}























