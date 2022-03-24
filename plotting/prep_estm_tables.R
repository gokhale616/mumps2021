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

# hard coding some parameter values for use downstream within estimation of R0s and Rps
# calculate the reproductive numbers 
params_for_R0 <- c(N = 100e6, nu = 1/80, p = 0, ad = age_class_duration)
params_for_Rp <- c(N = 100e6, nu = 1/80, p = 1, ad = age_class_duration)


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
    select(., -npar)
  
}


all_result_df <- (
  map_dfr(1:length(result_list), mk_result_df) %.>% 
    mutate(., 
           d_AIC = AIC - min(AIC)) %.>% 
    mutate(., 
           best_fit_covar = ifelse(AIC == min(AIC), 1, 0)) %.>% 
    ungroup(.) #%.>% 
    #select(., -AIC)
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
           p_intro = ifelse(p_intro == 3, "[15, 25)", p_intro),
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
                                 Quantity == "AIC"~"$AIC$",
                                 Quantity == "loglik"~"$log\\big(\\mathcal{L}(\\Theta)\\big)$",
                                 Quantity == "R0"~"$R_0$",
                                 Quantity == "Rp"~"$R_p$",
                                 Quantity == "impact"~"$\\xi$",
                                 Quantity == "beta1"~"$\\beta_1$", 
                                 Quantity == "sigma"~"$\\sigma^{-1}$ (Days)", 
                                 Quantity == "dwan"~ "$\\delta^{-1}$ (Years)",  
                                 Quantity == "epsilon2"~ "$\\epsilon$", 
                                 Quantity == "p_intro"~ "$p_{intro}$ (Age cohort)", 
                                 Quantity == "t_intro"~ "$t_{intro}$", 
                                 Quantity == "vacc_covariate"~ "Booster shape"
                                 )) %.>% 
      mutate(., `Parameter/Quantity` = Quantity) %.>% 
      select(., -Quantity) %.>% 
      select(., c(5, 2, 4, 3, 1)) %.>% 
      slice(., c(3, 1, 7, 9, 10, 6, 2, 11, 4, 5, 12, 13))
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



#load confidence intervals
CI_result_path <- "../result_data/CI"

list.files(path = CI_result_path, 
           full.names = TRUE) %.>% 
  lapply(., load,  envir = .GlobalEnv)



# make a single data frame for CIs and back transfrom
all_hypo_bootstrap_res <- (
  leaky2_bootstrap_res %.>% 
    bind_rows(., 
              waning_bootstrap_res) %.>% 
    bind_rows(., gwaning_bootstrap_res) %.>% 
    bind_rows(., noloss_bootstrap_res) %.>% 
    mutate_at(., 
              .vars = vars(starts_with("q_age_"), starts_with("rho_age_"), beta1, epsilon2), .funs = invlogit) %.>% 
    mutate_at(., .vars = vars(starts_with("psi_"), 
                              sigma, dwan, t_intro), .funs = exp) %.>% 
    # some modifications to avoid default conflicts downstream
    mutate(., 
           dwan     = ifelse(hypo_name %in% c("waning", "gwaning"), dwan, Inf), 
           epsilon2  = ifelse(hypo_name == "leaky2", epsilon2, 0), 
           t_intro  = ifelse(hypo_name == "leaky2", t_intro, 3000), 
           p_intro  = ifelse(hypo_name == "leaky2", 3, 6),
           ) 
  )


# calculate additional parameters before summarizing

if(paste0(CI_result_path, "/all_hypo_bootstrap_res_quants.rds") == FALSE) {
  message("Preparing bootstrapped derived quatitites! Takes 40secs on 8cores")
  
  plan(multiprocess)
  
  all_hypo_bootstrap_res_quants <- (
    future_map_dfr(1:nrow(all_hypo_bootstrap_res), function(x) {
      
      # hypo_name
      hypo_name <- (
        all_hypo_bootstrap_res %.>% 
          slice(., x) %.>% 
          select(., hypo_name) %.>%
          unlist(.)
      )  
      
      
      # pick the parameters   
      param_int <-(
        all_hypo_bootstrap_res %.>% 
          slice(., x) %.>% 
          select(., -c(hypo_name, covar_name)) %.>% 
          unlist(.) %.>% 
          sim_p_vals(.)
      )
      
      if(hypo_name == "gwaning") {
        R0 <- c(param_int, params_for_R0) %.>% calculate_R0_mq(.)$reprodutive_number 
        
        Rp <- c(param_int, params_for_Rp) %.>% calculate_R0_mq(.)$reprodutive_number   
      } else {
        R0 <- c(param_int, params_for_R0) %.>% calculate_R0_mq(.)$reprodutive_number 
        
        Rp <- c(param_int, params_for_Rp) %.>% calculate_R0_mq(.)$reprodutive_number   
      }
      
      
      all_hypo_bootstrap_res %.>% 
        slice(., x) %.>% 
        mutate(., 
               R0 = R0, 
               Rp = Rp, 
               impact = 1-Rp/R0)
      
      
      
    })
  )
  
  save(all_hypo_bootstrap_res_quants, file = "../result_data/CI/all_hypo_bootstrap_res_quants.rds")
    
} else {
  message("loading bootstrapped quantitites from results!")
}

load(paste0(CI_result_path, "/all_hypo_bootstrap_res_quants.rds"))



# get rid of covariate column and summarize 
all_hypo_bootstrap_res_quants_l <- (
  all_hypo_bootstrap_res_quants %.>% 
    mutate(., sigma = 365.25/sigma) %.>% 
    select(., -c(covar_name, `.id`)) %.>%
    gather(., key = "parameter", value = "est_val", -hypo_name)
  )



# CI's find values for the table 
all_hypo_CI_main <- (
  all_hypo_bootstrap_res_quants_l %.>% 
    transmute(., 
              hypo_name = hypo_name, 
              Parameter = parameter, 
              est_val = est_val) %.>% 
    group_by(., hypo_name, Parameter) %.>%
    summarise(., 
              qs = quantile(est_val, c(0.025, 0.975), na.rm = TRUE), 
              prob = c("0.025", "0.975"), 
              .groups = 'drop') %.>%
    ungroup(.) %.>% 
    mutate(., qs = ifelse(qs == 0 | is.infinite(qs) == TRUE, NA, qs)) %.>% 
    spread(., key = prob, value = qs) %.>% 
    filter(., Parameter %in% c("R0", "Rp", "impact", "beta1", 
                               "sigma", "dwan", "epsilon2", "p_intro", "t_intro")) %.>% 
    pivot_wider(., id_cols = Parameter, names_from = hypo_name, values_from = c("0.025", "0.975")) %.>%
    slice(., 6, 7, 4, 2, 3, 1, 8) %.>% 
    select(., 1, 4, 8, 5, 9, 2, 6, 3, 7) %.>% 
    mutate(.,
           `0.025_waning` = ifelse(Parameter == "dwan", 111.4, `0.025_waning`),
           Parameter = case_when(Parameter == "R0"~"$R_0$",
                                    Parameter == "Rp"~"$R_p$",
                                    Parameter == "impact"~"$\\xi$",
                                    Parameter == "beta1"~"$\\beta_1$", 
                                    Parameter == "sigma"~"$\\sigma^{-1}$ (Days)", 
                                    Parameter == "dwan"~ "$\\delta^{-1}$ (Years)",  
                                    Parameter == "epsilon2"~ "$\\epsilon$")
           ) %.>% 
    transmute(., 
              `Parameter/Quantity` = Parameter, 
              `2.5%1` = `0.025_noloss`, 
              `97.5%1` = `0.975_noloss`, 
              `2.5%2` = `0.025_waning`, 
              `97.5%2` = `0.975_waning`, 
              `2.5%3` = `0.025_gwaning`, 
              `97.5%3` = `0.975_gwaning`, 
              `2.5%4` = `0.025_leaky2`, 
              `97.5%4` = `0.975_leaky2`)
    
    
  )
  

CI_kbl <- (
  all_hypo_CI_main %.>%   
  kbl(., 
      align = "c", 
      digits = 4, 
      linesep = "",
      booktabs = T, 
      format = "html", 
      caption = "Model specific 95% Bootstraped Confidence Intervals", 
      escape = FALSE) %.>% 
  kable_styling(.,
                position='left', full_width = F,
                latex_options=c('striped', 'HOLD_position', "scale_down")) %.>%
  add_header_above(., 
                   c(" " = 1, "No Loss" = 2, 
                     "Waning (Exponential)" = 2, "Waning (Erlang, N = 3)" = 2, "Leaky" = 2))
)




# table for the main manuscript
# est_mle_kbl <- make_table(table_data =  table_hypo_compare)




# at the end load the treated covars so that this can they can be referenced while making plots 
# make this happen nicely  
source("../fit/treat_vacc_covar.R", chdir = TRUE)


