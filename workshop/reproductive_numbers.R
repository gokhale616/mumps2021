source("../00/src.R", chdir = TRUE)



#test_R0_p_vals <- c(param_vals_est, N = 100e6, nu = 1/80, p = 0.86, ad = age_class_duration)

# test_R0_p_vals["epsilon2"] <- c(0.5)

quick_select <- function(p_vect_data, pattern) {
  
  p_vect_data %.>% 
    select(., starts_with(pattern)) %.>% 
    unlist(.) 
    
}


# this functions calculates the pre-requisites for the next generateion operator analysis 
gather_R0_prereqs <- function(p_vals) {
  
  with(as.list(p_vals), {
  
  # definitions for convenience
  Births <- (nu*N)  
  delta <- (1/dwan)  
  
  p_vals_as_data <- (
    p_vals %.>% 
    as.list(.) %.>% 
    as_tibble(.))
  
  # age specific poi given contact
  q_mult <- (
    p_vals_as_data %.>% 
    quick_select(., pattern = "q_age")
    )
  
  
  # convert age duration to aging rates
  omega_mult <- 1/(
    p_vals_as_data %.>% 
    quick_select(., pattern = "ad") 
  ) %.>% setNames(., nm = sprintf("omega_%d", 1:5))
  
  # These variables are useful in calculating the 
  # equilibria near the disease free states 
  # (omega_i + delta)
  omega_mult_delta <- (omega_mult + delta)
  # omega_{i-1}
  omega_mult_1 <- c(1, omega_mult[1:4])
  # \prod_{k = 1}^i \frac{omega_{i-1}}{(omega_i + delta)}
  omega_mult_delta_cum_prod <- cumprod(omega_mult_1/omega_mult_delta)
  # \sum_{l = 1}^k \prod_{k = 1}^i \frac{omega_{i-1}}{(omega_i + delta)}
  omega_mult_delta_cum_prod_sum <- cumsum(omega_mult_delta_cum_prod)
  
  # Three state variables at their disease free equilibrium
  S_de  <- (1-(1-alpha)*p)*Births/omega_mult %.>% setNames(., nm = sprintf("S_de%d", 1:5))
  V_de  <- (1-alpha)*p*Births*omega_mult_delta_cum_prod %.>% setNames(., nm = sprintf("V_de%d", 1:5))
  Sv_de <- (1-alpha)*p*Births*(delta/omega_mult)*omega_mult_delta_cum_prod_sum %.>% 
    setNames(., nm = sprintf("Sv_de%d", 1:5))
  
  # Contact matrix
  C_mat <- (
    p_vals_as_data %.>% 
      quick_select(., pattern = "Cv") %.>% 
      matrix(., nrow = 5)
  )

  list(q_vec = q_mult, omega_vec = omega_mult, 
       S_de = S_de, V_de = V_de, Sv_de = Sv_de, 
       N_vec = S_de + V_de + Sv_de,
       C_mat = C_mat, delta = delta 
       )
  
  })
  
}


# utility functions around the contact matrix - Know that this bit of analysis is going to useful for quickly 
# accessing the next generation matrix for use contact matrix, calculate R0 based on the NextGen method and 
# so on, but very basic...

# Dirac delta function 
dirac_delta <- function(i, j) return(ifelse(i == j, 1, 0)) 


calculate_R0_mq <- function(p_vals) {
  
  with(as.list(p_vals), {
    # browser()
    
    # generate the pre-reqs for further analysis 
    R0_prereq <- gather_R0_prereqs(p_vals = p_vals)
    
    # equilibrium states 
    N  <- R0_prereq$N_vec
    S  <- R0_prereq$S_de
    V  <- R0_prereq$V_de
    Sv <- R0_prereq$Sv_de
    
    # aging rates
    omega <- R0_prereq$omega_vec
    # poi 
    q <- R0_prereq$q_vec
    # waning rate
    delta <- R0_prereq$delta
    
    #contact matrix
    Cmat <- R0_prereq$C_mat
    
    # Define emplty matrices of F and V to calculate the the NextGen Matrix  
    NA_mat <- matrix(NA, nrow = nrow(C), ncol = ncol(C)) 
    zero_mat <- matrix(0, nrow = nrow(C), ncol = ncol(C))
    
    
    # Intialize 
    F_mat_nz <- NA_mat
    
    # populate the F_mat  
    for(i in seq_along(q)) {
        for(j in seq_along(q)) {
          F_mat_nz[i,j] <- q[i]*Cmat[i,j]*(S[i] + Sv[i] + epsilon2*V[i])/N[j]
          }
      }
    
    # Calculate the non-zero block of the F_mat matrix
    # Define the F matrix
    F_mat = rbind(cbind(zero_mat, F_mat_nz), cbind(zero_mat, zero_mat))
    
    # Initialize the block matrices of the V matrix, (explain what hey are later) 
    V1 <- NA_mat 
    V2 <- zero_mat
    V3 <- NA_mat
    V4 <- NA_mat
    
    for(i in seq_along(q)) {
      for(j in seq_along(q)) {
        if(i == 1) {
          V1[i,j] <-  (omega[i] + sigma)*dirac_delta(i,j)
          V3[i,j] <- -sigma*dirac_delta(i,j)
          V4[i,j] <-  (omega[i] + gamma)*dirac_delta(i,j)  
        } else {
          V1[i,j] <-  (omega[i] - omega[i-1] + sigma)*dirac_delta(i,j)
          V3[i,j] <- -sigma*dirac_delta(i,j)
          V4[i,j] <-  (omega[i] - omega[i-1] + gamma)*dirac_delta(i,j)  
        }
      }
    }
    
    V_mat = rbind(cbind(V1, V2), cbind(V3, V4))
    
    # browser()
    # calculate the next generation matrix
    K_mat <- F_mat%*%solve(V_mat)
    
    R0 <- max(eigen(K_mat)$values)
    
    other_output <- eigen(K_mat) 
    
    list(reprodutive_number = R0,
         F_mat = F_mat, V_mat = V_mat, K_mat = K_mat, 
         other_output = other_output)
    
  })
  
  
}


#R0_deets <- calculate_R0_mq(p_vals = test_R0_p_vals)$reprodutive_number


