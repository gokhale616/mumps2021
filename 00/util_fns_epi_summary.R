# source("../00/src.R", chdir = TRUE)

# test_R0_p_vals <- c(param_vals_est, N = 100e6, nu = 1/80, p = 0.86, ad = age_class_duration)

# test_R0_p_vals["epsilon2"] <- c(0.5)

quick_select <- function(p_vect_data, pattern) {
  
  p_vect_data %.>% 
    select(., starts_with(pattern)) %.>% 
    unlist(.) 
    
}


# this functions calculates the pre-requisites for the next generation operator analysis - 
# exponential waning model - NEW VERSION CONATAINS 2 VACCINE DOSES
gather_R0_prereqs <- function(p_vals) {
  
  with(as.list(p_vals), {
    # browser()
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
    
    
    # calculate the equilibria for the state variables 
    ev1 <- (1-alpha)*p1
    ev2 <- (1-alpha)*p2
    
    # set three empty vectors to fill with the endemic equilibria 
    S_de <- rep(NA, 5) %.>% setNames(., nm = sprintf("S_de%d", 1:5))
    V_de <- rep(NA, 5) %.>% setNames(., nm = sprintf("V_de%d", 1:5))
    Sv_de <- rep(NA, 5) %.>% setNames(., nm = sprintf("Sv_de%d", 1:5))
    
    # fill in the values ------
    # [0,5)
    S_de[1]  <- (1-ev1)*Births/omega_mult[1]
    V_de[1]  <- ev1*Births/omega_mult_1[1]
    Sv_de[1] <- (delta/omega_mult[1])*V_de[1]
    # [5,15)
    S_de[2]  <- (1-ev2)*omega_mult[1]/omega_mult[2]*S_de[1] 
    V_de[2]  <- ev2*omega_mult_1[1]/omega_mult_delta[2]*S_de[1]
    Sv_de[2] <- (delta*V_de[2] + omega_mult[1]*Sv_de[1])/omega_mult[2]
    
    # [15, 25), [25, 40), >40
    for(i in 3:5) {
      
      S_de[i]  <- omega_mult[i-1]/omega_mult[i]*S_de[i-1]
      V_de[i]  <- omega_mult[i-1]/omega_mult_delta[i]*V_de[i-1]
      Sv_de[i] <- (delta*V_de[i] + omega_mult[i-1]*Sv_de[i-1])/omega_mult[i]  
      
    }
    
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
    # browser()
    list(reprodutive_number = R0,
         F_mat = F_mat, V_mat = V_mat, K_mat = K_mat, 
         other_output = other_output)
    
  })
  
  
}




# this function calculates the effective R0 for this system
calculate_Reff_mq <- function(p_vals, 
                              S_1, S_2, S_3, S_4, S_5, 
                              V_1, V_2, V_3, V_4, V_5, 
                              N_1, N_2, N_3, N_4, N_5, 
                              t) {
  
  with(as.list(c(p_vals, 
                 S_1, S_2, S_3, S_4, S_5, 
                 V_1, V_2, V_3, V_4, V_5, 
                 N_1, N_2, N_3, N_4, N_5, 
                 t)), {
    
    # browser()
    # generate the pre-reqs for further analysis 
    R0_prereq <- gather_R0_prereqs(p_vals = p_vals)
    
    # equilibrium states 
    N  <- c(N_1, N_2, N_3, N_4, N_5)
    S  <- c(S_1, S_2, S_3, S_4, S_5)
    V  <- c(V_1, V_2, V_3, V_4, V_5)
    
    # aging rates
    omega <- R0_prereq$omega_vec
    # poi 
    q <- R0_prereq$q_vec
    # waning rate
    delta <- R0_prereq$delta
    
    #contact matrix
    Cmat <- R0_prereq$C_mat
    
    # Define emplty matrices of F and V to calculate the the NextGen Matrix  
    NA_mat <- matrix(NA, nrow = nrow(Cmat), ncol = ncol(Cmat)) 
    zero_mat <- matrix(0, nrow = nrow(Cmat), ncol = ncol(Cmat))
    
    
    # Intialize 
    F_mat_nz <- NA_mat
    
    # effective leakiness when dealing with two pathogen hypothesis
    if(t < t_intro) {
      epsilon2_eff = 0
    } else {
      epsilon2_eff = epsilon2
    }
    
    #browser()
    # populate the F_mat  - make sure to add the seasonality
    for(i in seq_along(q)) {
      for(j in seq_along(q)) {
        if(i == 1 & j == 1) {
          F_mat_nz[i,j] <- q[i]*Cmat[i,j]*(1-beta1*sin(2*pi*t))*(S[i] + (epsilon1 + epsilon2_eff)*V[i])/N[j]
        } else {
          F_mat_nz[i,j] <- q[i]*Cmat[i,j]*(S[i] + epsilon2*V[i])/N[j]
        }
      }
    }
    #browser()
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
    # return spectral radius aka R0
    max(eigen(K_mat)$values) %.>% 
      ifelse(is.nan(.) == TRUE, NA, .)
    })
  
}

# this function calculates average transmission rate 
calculate_average_beta <- function(p_vals, 
                                   S_1, S_2, S_3, S_4, S_5, 
                                   V_1, V_2, V_3, V_4, V_5, 
                                   I_1, I_2, I_3, I_4, I_5,
                                   N_1, N_2, N_3, N_4, N_5,
                                   t) {
  
  
  with(as.list(c(p_vals, 
                 S_1, S_2, S_3, S_4, S_5, 
                 V_1, V_2, V_3, V_4, V_5, 
                 I_1, I_2, I_3, I_4, I_5,
                 N_1, N_2, N_3, N_4, N_5,
                 t)), {
                   

     # generate the pre-reqs for further analysis 
     R0_prereq <- gather_R0_prereqs(p_vals = p_vals)
     
     # equilibrium states 
     N <- c(N_1, N_2, N_3, N_4, N_5)
     S <- c(S_1, S_2, S_3, S_4, S_5)
     V <- c(V_1, V_2, V_3, V_4, V_5)
     I <- c(I_1, I_2, I_3, I_4, I_5)
     
     # aging rates
     omega <- R0_prereq$omega_vec
     # poi 
     q <- R0_prereq$q_vec
     # waning rate
     delta <- R0_prereq$delta
     
     #contact matrix
     Cmat <- R0_prereq$C_mat                 
     
     # effective leakiness when dealing with the two pathogen hypothesis
     if(t < t_intro) {
       epsilon2_eff = 0
     } else {
       epsilon2_eff = epsilon2
     }
     
     # Calculate the average transmission rate
     # Initialize an NA matrix for your own good! 
     NA_matrix <- matrix(NA, nrow = nrow(Cmat), ncol = ncol(Cmat)) 
     
     beta <-  NA_matrix
     contact_scale_factor <- NA_matrix
     lambda <- NA_matrix
     
     for(i in seq_along(q)) {
       for(j in seq_along(q)) {
         if(i == 1 & j == 1) {
           beta[i,j] <- q[i]*Cmat[i,j]*(1-beta1*sin(2*pi*t))*(S[i] + (epsilon1 + epsilon2_eff)*V[i])*I[j]/N[j]
           contact_scale_factor[i,j] <- (S[i] + (epsilon1 + epsilon2_eff)*V[i])*I[j]/N[j]
           lambda[i,j] <- q[i]*Cmat[i,j]*(1-beta1*sin(2*pi*t))*I[j]/N[j]
          } else {
           beta[i,j] <- q[i]*Cmat[i,j]*(S[i] + epsilon2*V[i])*I[j]/N[j]
           contact_scale_factor[i,j] <- (S[i] + (epsilon1 + epsilon2_eff)*V[i])*I[j]/N[j]
           lambda[i,j] <- q[i]*Cmat[i,j]*(1-beta1*sin(2*pi*t))*I[j]/N[j]
         }
       }
     }  
     
     # calculate the average transmission rate
     average_beta <- sum(beta)/sum(contact_scale_factor)
     
     
     # calculate the mean age at which individuals from i^{th} age cohort get infected
     Ai <- 1/colSums(lambda)
     
     # calculate the overall  mean age of infection 
     mean_age_at_infection <- sum(Ai*I/sum(I))
     
     
     list(B = average_beta, A = mean_age_at_infection)

})
  
}

# this function produces summary measure of simulated epidemics
summarize_epidemiology <- function(traj_cov_data, p_vals) {
  # life expectancy is assumed to be 80 years
  traj_cov_data %.>% 
    map_dfr(1:nrow(.), function(c, x = .) {
      # browser()
      
      int_x <- (
        x %.>% 
          slice(., c) %.>% 
          mutate_all(., .funs = function(x){ifelse(x < 0, 0, x)})
        )
      
      res <-(
        int_x %.>% 
          transmute(., 
                    `.id` = `.id`,
                    year = year,
                    # scale population sizes to represent system compartments per 100000
                    Ns_1 = 1e5/N_1, Ns_2 = 1e5/N_2, Ns_3 = 1e5/N_3,
                    Ns_4 = 1e5/N_4, Ns_5 = 1e5/N_5,
                    # scale the susceptible compartments
                    Ss_1  = S_1*Ns_1, Ss_2  = S_2*Ns_2, Ss_3  = S_3*Ns_3,
                    Ss_4  = S_4*Ns_4, Ss_5  = S_5*Ns_5,
                    # scale the susceptible compartments
                    Sp_1  = S_1/N_1, Sp_2  = S_2/N_2, Sp_3  = S_3/N_3,
                    Sp_4  = S_4/N_4, Sp_5  = S_5/N_5,
                    # scale the vaccinated compartments
                    Vs_1  = V_1*Ns_1, Vs_2  = V_2*Ns_2, Vs_3  = V_3*Ns_3,
                    Vs_4  = V_4*Ns_4, Vs_5  = V_5*Ns_5,
                    # define and scale the recovered compartments
                    Rs_1  = (N_1 - (S_1 + V_1 + I1_1 + I2_1 + E1_1 + E2_1))*Ns_1,
                    Rs_2  = (N_2 - (S_2 + V_2 + I1_2 + I2_2 + E1_2 + E2_2))*Ns_2, 
                    Rs_3  = (N_3 - (S_3 + V_3 + I1_3 + I2_3 + E1_3 + E2_3))*Ns_3,
                    Rs_4  = (N_4 - (S_4 + V_4 + I1_4 + I2_4 + E1_4 + E2_4))*Ns_4,
                    Rs_5  = (N_5 - (S_5 + V_5 + I1_5 + I2_5 + E1_5 + E2_5))*Ns_5,
                    # proportion of I2 age structured
                    I2p_1 = I2_1/(I1_1+I2_1), I2p_2 = I2_2/(I1_2+I2_2), I2p_3 = I2_3/(I1_3+I2_3),
                    I2p_4 = I2_4/(I1_4+I2_4), I2p_5 = I2_5/(I1_5+I2_5),
                    # proportion of total I2
                    I2p = (I2_1 + I2_2 + I2_3 + I2_4 + I2_5)/(I1_1+I2_1 + I1_2+I2_2 + I1_3+I2_3 + I1_4+I2_4 + I1_5+I2_5),
                    # scale the true new cases compartments
                    Cs_1 = (C_1)*Ns_1, Cs_2 = (C_2)*Ns_2, Cs_3 = (C_3)*Ns_3, 
                    Cs_4 = (C_4)*Ns_4, Cs_5 = (C_5)*Ns_5, 
                    # scaled total cases 
                    Cs = (C_1 + C_2 + C_3 + C_4 + C_5)*1e5/(N_1 + N_2 + N_3 + N_4 + N_5),
                    # define the total infectious compartments
                    Is_1 = (I1_1+I2_1)*Ns_1, Is_2 = (I1_2+I2_2)*Ns_2, Is_3 = (I1_3+I2_3)*Ns_3,
                    Is_4 = (I1_4+I2_4)*Ns_4, Is_5 = (I1_5+I2_5)*Ns_5,
                    # total infectious 
                    Is = (I1_1 + I1_2 + I1_3 + I1_4 + I1_5 + I2_1 + I2_2 + I2_3 + I2_4 + I2_5)*1e5/(N_1 + N_2 + N_3 + N_4 + N_5),
                    # calculate effective R0 for the system
                    Reff = calculate_Reff_mq(p_vals = p_vals, 
                                             N_1 = N_1, N_2 = N_2, N_3 = N_3, 
                                             N_4 = N_4, N_5 = N_5, 
                                             S_1 = S_1, S_2 = S_2, S_3 = S_3, 
                                             S_4 = S_4, S_5 = S_5,
                                             V_1 = V_1, V_2 = V_2, V_3 = V_3, 
                                             V_4 = V_4, V_5 = V_5,
                                             t = year))) 
        
      
      
      res
    
    })
  
  
}


# this is a convenience function for to prepare data for plotting
prep_as_plot_trajectory <- function(traj_data, init_year) {
  
  traj_data %.>% 
    filter(., year > init_year) %.>% 
    select(., -`.id`) %.>% 
    gather(., key = "comp", value = "count", -year) %.>% 
    mutate(., 
           comp_exp = str_split_fixed(comp, pattern = "_", n = 2)[,1], 
           age_number = str_split_fixed(comp, pattern = "_", n = 2)[,2], 
    ) %.>% 
    mutate(., 
           age_cohort = case_when(age_number == 1 ~ age_names[1], 
                                  age_number == 2 ~ age_names[2], 
                                  age_number == 3 ~ age_names[3], 
                                  age_number == 4 ~ age_names[4], 
                                  age_number == 5 ~ age_names[5], 
                                  TRUE            ~ age_number) %.>% as_factor(.)
    ) %.>% 
    select(., year, comp_exp, age_cohort, count)
  
}




