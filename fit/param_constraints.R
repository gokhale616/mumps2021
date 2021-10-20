# this script contains parameter constraints for  all the hypotheses

# parameter constraints for hypothesis on waning vaccine derived immunity
param_range_waning <- (
  list(
    lower = c(q_age_1 = 0, q_age_2 = 0, q_age_3 = 0, q_age_4 = 0, q_age_5 = 0,
              rho_age_1 = 0, rho_age_2 = 0, rho_age_3 = 0, rho_age_4 = 0, rho_age_5 = 0, rho_age_u = 0, 
              psi_1 = 0, psi_2 = 0, psi_3 = 0, psi_4 = 0, psi_5 = 0, psi_u = 0, 
              dwan = 1),
    
    upper = c(q_age_1 = 1, q_age_2 = 1, q_age_3 = 1, q_age_4 = 1, q_age_5 = 1, 
              rho_age_1 = 1, rho_age_2 = 1, rho_age_3 = 1, rho_age_4 = 1, rho_age_5 = 1, rho_age_u = 1, 
              psi_1 = 2, psi_2 = 2, psi_3 = 2, psi_4 = 2, psi_5 = 2, psi_u = 2, 
              dwan = 500)
    )
  )



