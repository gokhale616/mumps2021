# load all of the result objects in the global environment 
mle_result_path <- "../result_data/mle"

list.files(path = mle_result_path, 
           full.names = TRUE)[-25] %.>% 
  lapply(., load,  envir = .GlobalEnv)


# One of the waning hypothesis was wrongly named to noloss
# correct the hypothesis name 
mle_waning_constant_is$Hypothesis <- "waning_constant"

# form a single list of all of the result objects to loop over later
result_list <- (
  list(
    mle_gamma_waning_slow_is, mle_gamma_waning_sigmoid_is, 
    mle_gamma_waning_rapid_is, mle_gamma_waning_constant_is,
    mle_waning_slow_is, mle_waning_sigmoid_is, mle_waning_rapid_is, mle_waning_constant_is,
    mle_leaky2_2_slow_is, mle_leaky2_2_sigmoid_is, mle_leaky2_2_rapid_is, mle_leaky2_2_constant_is,
    mle_leaky2_3_slow_is, mle_leaky2_3_sigmoid_is, mle_leaky2_3_rapid_is, mle_leaky2_3_constant_is,
    mle_leaky2_4_slow_is, mle_leaky2_4_sigmoid_is, mle_leaky2_4_rapid_is, mle_leaky2_4_constant_is, 
    mle_slow_is, mle_sigmoid_is, mle_rapid_is, mle_constant_is
  )
)


mk_result_df <- function(c = 1, res = result_list) {
  # to look into the function at a specific iteration
  # if (c == 20) {browser()}
  
  # collect qualitative covariates
  extra_params <- res[[c]]$Hypothesis %.>% str_split(., pattern = "_") %.>% unlist(.)
  
  if(length(extra_params) == 3) {
    hypo_covar <- extra_params[1]
    p_intro    <- extra_params[2] %.>% as.numeric(.)
    vacc       <- extra_params[3]
  } else {
    hypo_covar <- extra_params[1]
    p_intro    <- NA
    vacc       <- extra_params[2]
  }
  
  # calculate the reproductive numbers 
  params_for_R0 <- c(N = 100e6, nu = 1/80, p = 0, ad = age_class_duration)
  params_for_Rp <- c(N = 100e6, nu = 1/80, p = 1, ad = age_class_duration)
  
  if(hypo_covar == "gwaning") {
    R0 <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_R0) %.>% 
      calculate_R0_mq(.)$reprodutive_number 
    
    Rp <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_Rp) %.>% 
      calculate_R0_mq(.)$reprodutive_number   
  } else {
      R0 <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_R0) %.>% 
        calculate_R0_mq(.)$reprodutive_number 
      
      Rp <- c(res[[c]]$DEobj$optim$bestmem %.>% sim_p_vals(.), params_for_Rp) %.>% 
        calculate_R0_mq(.)$reprodutive_number   
  }
  
  
  # collate results in a dataframe
  res[[c]]$DEobj$optim$bestmem %.>% 
    as.list(.) %.>% 
    as_tibble(.) %>% 
    mutate(., 
           R0 = R0,
           Rp = Rp,
           impact = 1-Rp/R0,
           loglik = -res[[c]]$DEobj$optim$bestval, 
           npar = res[[c]]$DEobj$optim$bestmem %.>% length(.), 
           AIC = calculate_aic(loglik, npar), 
           hypothesis = hypo_covar, 
           p_intro    = p_intro,
           vacc_covariate = vacc) %.>% 
    select(., -c(loglik, npar))
  
}


all_result_df <- (
  map_dfr(1:length(result_list), mk_result_df) %.>% 
    mutate(., 
           d_AIC = AIC - min(AIC)) %.>% 
    mutate(., 
           best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) %.>% 
    select(., -AIC)
)


# table of estimates - process model 
table_hypo_compare <- (
  all_result_df %.>% 
    select(., -c(best_fit_covar, starts_with("q_age"), 
                 starts_with("rho_age"), starts_with("psi"))) %.>% 
    arrange(., d_AIC) %.>% 
    group_by(., hypothesis) %.>% 
    filter(., d_AIC == min(d_AIC)) %.>% 
    ungroup(.) %.>%
    mutate(., 
           sigma = 365.25/sigma, 
           hypothesis = case_when(hypothesis == "waning" ~  "Waning (Exponential)", 
                                  hypothesis == "gwaning" ~  "Waning (Erlang, n = 3)", 
                                  hypothesis == "leaky2" ~  "Leaky", 
                                  hypothesis == "noloss" ~  "No Loss")) %.>% 
    mutate_if(., is.numeric, function(x){round(x, digits = 3)}) %.>%   
    mutate(., 
           vacc_covariate = str_to_title(vacc_covariate),
           hypothesis = str_to_title(hypothesis)) %.>% 
    gather(., 
           key = "Quantity", value = "Estimate", 
           -c(hypothesis)) %.>% 
    spread(., key = hypothesis, value = Estimate) %.>% 
     mutate(., 
            Quantity = case_when(Quantity == "d_AIC"~"$\\Delta AIC$",
                                 Quantity == "R0"~"$R_0$",
                                 Quantity == "Rp"~"$R_p$",
                                 Quantity == "impact"~"Vaccine impact ($\\psi$)",
                                 Quantity == "beta1"~"$\\beta_1$", 
                                 Quantity == "sigma"~"$\\sigma^{-1}$ (days)", 
                                 Quantity == "dwan"~ "$\\delta^{-1}$ (years)",  
                                 Quantity == "epsilon2"~ "$\\epsilon$", 
                                 Quantity == "p_intro"~ "$p_{intro}$", 
                                 Quantity == "t_intro"~ "$t_{intro}$", 
                                 Quantity == "vacc_covariate"~ "Booster shape"
                                 )) %.>% 
      mutate(., `Parameter/Quantity` = Quantity) %.>% 
      select(., -Quantity) %.>% 
      select(., c(5, 2, 4, 3, 1)) %.>% 
      slice(., c(2, 7, 8, 5, 1, 9, 3, 4, 6, 10, 11))
)


make_table <- function(table_data) {
  browser()
  table_data %.>% 
    kbl(., 
        align = "c", 
        digits = 4, 
        booktabs = T, 
        format = "latex", 
        #caption = "Model specific parameter estimates and derived quantitities were obtained by maximizing the likelihood function", 
        caption = "blah",
        escape = FALSE) %.>% 
    kable_styling(.,
                  position='left', full_width = F,
                  #latex_options=c('striped', 'hold_position', "scale_down"))%.>%
                  latex_options=c('hold_position'))%.>%
    add_header_above(., 
                     bold = TRUE, 
                     c(" " = 1, "Model" = 4)) #%.>% 
    #landscape(.)
}


est_mle_kbl <- (
  table_hypo_compare %.>% 
  kbl(., 
      align = "c", 
      digits = 4, 
      linesep = "",
      booktabs = T, 
      format = "latex", 
      caption = "Model specific parameter estimates and derived quantitities were obtained by maximizing the likelihood function", 
      escape = FALSE) %.>% 
   kable_styling(.,
                 position='left', full_width = F,
                 latex_options=c('striped', 'HOLD_position', "scale_down")) %.>%
  add_header_above(., 
                   bold = TRUE, 
                   c(" " = 1, "Model" = 4))
  )


# table for the main manuscript
# est_mle_kbl <- make_table(table_data =  table_hypo_compare)




# at the end load the treated covars so that this can they can be referenced while making plots 
# make this happen nicely  
source("../fit/treat_vacc_covar.R", chdir = TRUE)


