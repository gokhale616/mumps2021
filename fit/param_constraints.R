# this script contains parameter constraints for  all the hypotheses

hypo_ind_params <- (
  list(
    lower = c(q_age_1 = 0, q_age_2 = 0, q_age_3 = 0, q_age_4 = 0, q_age_5 = 0, sigma = 365.25/25,
              rho_age_1 = 0, rho_age_2 = 0, rho_age_3 = 0, rho_age_4 = 0, rho_age_5 = 0, rho_age_u = 0, 
              psi_1 = 0, psi_2 = 0, psi_3 = 0, psi_4 = 0, psi_5 = 0, psi_u = 0, 
              beta1 = 0),
    
    upper = c(q_age_1 = 1, q_age_2 = 1, q_age_3 = 1, q_age_4 = 1, q_age_5 = 1, sigma = 365.25/12,
              rho_age_1 = 1, rho_age_2 = 1, rho_age_3 = 1, rho_age_4 = 1, rho_age_5 = 1, rho_age_u = 1, 
              psi_1 = 2, psi_2 = 2, psi_3 = 2, psi_4 = 2, psi_5 = 2, psi_u = 2, 
              beta1 = 1)
  )
)




# parameter constraints for hypothesis on waning vaccine derived immunity
param_range_waning <- (
  list(
    lower = c(hypo_ind_params$lower, dwan = 1),
    
    upper = c(hypo_ind_params$upper, dwan = 500)
    )
  )


# parameter constraints for hypothesis on leaky vaccine with imported infectious indiviuduals
param_range_leaky2 <- (
  list(
    lower = c(hypo_ind_params$lower, t_intro = 2000, epsilon2 = 0),
    
    upper = c(hypo_ind_params$upper, t_intro = 2005, epsilon2 = 1)
  )
  
)



