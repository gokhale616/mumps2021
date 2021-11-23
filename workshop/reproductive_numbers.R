source("../00/src.R", chdir = TRUE)

# R0 <- function(par) {
#   with(as.list(par %.>% sim_p_vals(.)), 
#        (beta_t*epsilon*nu*kappa)/(delta_i*(delta_v+beta_v*kappa))
#   )
# }




str_extract()

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
