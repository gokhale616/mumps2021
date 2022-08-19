# depends on ./cov_make_pomp_model.R
# defining the model snippets for pomp -----------------------------------------------------------------------

##### deterministic skeleton - for estimation #####
vseir_skel_est_gamma_n_2 <- Csnippet("
    // define  the 
    double delta = 1/dwan; 

    int i,j,k,l;
    int n = agegroups; 
    
    double lambda1;
    double lambda2;
    
    const double *Cvlocal = &Cv1;
    const double *lvlocal = &lcv1;
    const double *qvlocal = &q_age_1;
    
    const double *Slocal  = &S_1;
    const double *E1local = Slocal+n;
    const double *I1local = E1local+n;
    const double *E2local = I1local+n;
    const double *I2local = E2local+n;
    const double *V1local = I2local+n;
    const double *V2local = V1local+n;
    const double *Clocal  = V2local+n;
    
    double *DSlocal  = &DS_1;
    double *DE1local = DSlocal+n;
    double *DI1local = DE1local+n;
    double *DE2local = DI1local+n;
    double *DI2local = DE2local+n;
    double *DV1local = DI2local+n;
    double *DV2local = DV1local+n;
    double *DClocal  = DV2local+n;
    
    #define CM(J,K) Cvlocal[(J)+n*(K)]
    #define LV(K) lvlocal[(K)]
    #define QAGE(K) qvlocal[(K)]
    
    // Special Age-specific Covariates macro definitions  
    
    //#define N_AGE(J)   (__covars[__covindex[0]  + (J)])
    #define MU_AGE(J)  (__covars[__covindex[5]  + (J)])
    
    #define S(K) Slocal[(K)]
    #define E1(K) E1local[(K)]
    #define I1(K) I1local[(K)]
    #define E2(K) E2local[(K)]
    #define I2(K) I2local[(K)]
    #define V1(K) V1local[(K)]
    #define V2(K) V2local[(K)]
    #define C(K) Clocal[(K)]
    
    #define DS(K) DSlocal[(K)]
    #define DE1(K) DE1local[(K)]
    #define DI1(K) DI1local[(K)]
    #define DE2(K) DE2local[(K)]
    #define DI2(K) DI2local[(K)]
    #define DV1(K) DV1local[(K)]
    #define DV2(K) DV2local[(K)]
    #define DC(K)  DClocal[(K)]
    
    /* ===================================================================== */
    /* conditional introduction of Mumps genotype G in a specified age class */
    /* ===================================================================== */
    // define an integer valued position variable
    int ip_intro = floor(p_intro);
    double iota_v[n];
      
    for (i=0; i<n; i++) {
      if (i == ip_intro) {
        iota_v[i] = iota;
      } else {
          iota_v[i] = 0;
      }
    }
    
    /* ============================================================================== */
    /* augment the contact matrix to incorporate seasonality in the [5,15) sge cohort */
    /* ============================================================================== */      
    // declare and initialize
    double augCM[n][n];
    
    // modify to incorporate seasonality  
    for(k=0; k<n; k++) {
      for(l=0; l<n; l++) {
        if(k == 1 && l == 1) {
          augCM[k][l] = CM(k,l)*(1-beta1*sin(2*M_PI*t));
        } else {
          augCM[k][l] = CM(k,l);
        }
      }
    }
  
    
    /* ============================================================================= */
    /* ========================== defining epidemic process ======================== */
    /* ============================ with type 1 mortality ========================== */
    /* ============================================================================= */
    
    for (i = 0; i < n; i++) {
  
      lambda1 = 0.0;
      lambda2 = 0.0;
      
      if(i == 0) {
        
        // births occur only in the first age class - specially treat the [0,5) age cohort
        
        for (j = 0; j < n; j++) {
          
          lambda1 += q*QAGE(i)*augCM[i][j]*(I1(j))/N_AGE(j);
          if (t < t_intro) {
            lambda2 = 0;
          } else {
              lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j))/N_AGE(j);
          }
        }
            
        /* ============================ */
        /* balance transition equations */
        /* ============================ */
        
        // vaccinated and unvaccinated births - calculated using covariate population  
        
        double uv_births = (1-(1-alpha)*p1)*Births;
        double  v_births = (1-alpha)*p1*Births;
        
        DS(i)  = uv_births + 2*delta*V2(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);
        DE1(i) = lambda1*(S(i) + epsilon1*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E1(i);   
        DI1(i) = sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
        DE2(i) = lambda2*(S(i) + epsilon2*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E2(i);
        DI2(i) = sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
        
        DV1(i) = v_births - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V1(i);
        DV2(i) = 3*delta*V1(i) - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V2(i);
      
      } else if (i == 1) {
        
          // individuals age in and out of age cohort after the first
        
          for (j = 0; j < n; j++) {
            
            lambda1 += q*QAGE(i)*augCM[i][j]*(I1(j))/N_AGE(j);
            if (t < t_intro) {
              lambda2 = 0;
            } else {
                lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j))/N_AGE(j);
            }
          }
            
          /* ============================ */
          /* balance tranistion equations */
          /* ============================ */
          
          double uv_grads = LV(i-1)*(1-(1-alpha)*p2)*S(i-1);
          double  v_grads = LV(i-1)*(1-alpha)*p2*S(i-1);
          
          DS(i)  = uv_grads + 2*delta*V2(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);               
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          
          DV1(i) = v_grads + LV(i-1)*V1(i-1) - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V1(i);
          DV2(i) = LV(i-1)*V2(i-1) + 2*delta*V1(i) - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V2(i);
      
      } else {
        
          // individuals age in and out of age cohort after the first
          for (j = 0; j < n; j++) {
            lambda1 += q*QAGE(i)*augCM[i][j]*(I1(j))/N_AGE(j);
            if (t < t_intro) {
              lambda2 = 0;
            } else {
                lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j))/N_AGE(j);
            }
          }
          
          /* ============================ */
          /* balance transition equations */
          /* ============================ */
          
          DS(i)  = LV(i-1)*S(i-1)  + 2*delta*V2(i)    - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);             
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*(V1(i)+V2(i))) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          DV1(i) = LV(i-1)*V1(i-1) - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V1(i);
          DV2(i) = LV(i-1)*V2(i-1) + 2*delta*V1(i) - (epsilon1*lambda1 + epsilon2*lambda2 + 2*delta + LV(i) - MU_AGE(i))*V2(i);
          
      }
      
      /* ========================= */
      /* Keep account of new cases */
      /* ========================= */
      DC(i) = gamma*(I1(i)+I2(i)); 
    }
")


##### initial value distribution - for estimation #####
vseir_init_est_gamma_n_2 <- Csnippet("
  double *S = &S_1;
  double *E1 = &E1_1;
  double *I1 = &I1_1;
  double *E2 = &E2_1;
  double *I2 = &I2_1;
  double *V1 = &V1_1;
  double *V2 = &V2_1;
  
  double *C = &C_1;
  
  int n = agegroups;
  
  #define N_AGE(J)  (__covars[__covindex[0]  + (J)])
  
  for (int i = 0; i < n; i++) {
    S[i]  = nearbyint(S_0*N_AGE(i));
    E1[i] = 0;
    I1[i] = nearbyint(I1_0*N_AGE(i));
    E2[i] = 0;
    I2[i] = nearbyint(I2_0*N_AGE(i));
    
    V1[i]  = 0;
    V2[i]  = 0;
    
    C[i]  = 0;
  }
")


##### parameter names ######
state_names_est_gamma_n_2 <- c(sprintf("S_%d",1:5), 
                               sprintf("E1_%d",1:5), sprintf("I1_%d",1:5),
                               sprintf("E2_%d",1:5), sprintf("I2_%d",1:5),
                               sprintf("V1_%d",1:5), sprintf("V2_%d",1:5), 
                               sprintf("C_%d", 1:5))



# testing parameter values 
##### build the pomp object #####
# Here the model is defined at an annual scale with a time step of that of a day (dt = 1/365.25)
# When extrapolating the simulation, default scale of simulation has been set to annual (temp_scale = 1)
# This can of course be changed as required by passing arguments

make_gamma_n_2_pomp <- function(..., 
                            start_t = 1800, 
                            extrapolate_simulation = FALSE, dt = 1/365.25,  
                            extra_start_t = 1900, extra_end_t = 2100, temp_scale = 1, 
                            covar, 
                            incidence_data = mumps_case_reports) {
  
  message("Assuming gamma (n = 2) distributed waning rate.")
  
  # changing the column names and subsetting to read into pomp
  mumps_inc_data <- (
    incidence_data %.>%  
      select(., -total) 
  )
  
  colnames(mumps_inc_data) <- c("year", sprintf("D_%d", 1:5), "D_u")
  
  # browser()
  
  if(extrapolate_simulation == TRUE & start_t > extra_start_t) {
    stop("'extra_start_t' needs to be greater than 'start_t'")
  }
  
  if(extrapolate_simulation == TRUE) {
    
    message(paste0("Extrapolating simulation at scale: ", temp_scale*365.25, " days") )
    # browser()
      po <- (
        tibble(year = c(seq(extra_start_t, extra_end_t, by = temp_scale)), 
               D_1 = NA, D_2 = NA, D_3 = NA, D_4 = NA, D_5 = NA, 
               D_u = NA) %.>% 
      arrange(., year) %.>% 
      pomp(., 
           times = "year", t0 = start_t, 
           skeleton = vectorfield(vseir_skel_est_gamma_n_2),
           rinit = vseir_init_est_gamma_n_2,
           rmeasure = vseir_rmeas,
           dmeasure = vseir_dmeas,
           accumvars = c(sprintf("C_%d", 1:5)),
           covar = covariate_table(covar, times = "Year", order = "constant"),
           covarnames = c(sprintf("N_%d", 1:5), sprintf("MU_%d", 1:5), 
                          sprintf("IN_%d", 1:5), sprintf("OUT_%d", 1:5), 
                          sprintf("p%d", 1:2), "Births", "eta_a"),
           statenames = state_names_est_gamma_n_2,
           paramnames = param_names_est, ...)
      ) 
  
  } else {
    
    message("Generating simulations at a scale and range identical to the input data")
    po <- (
      mumps_inc_data %.>%
      pomp(., 
           times = "year", t0 = start_t, 
           skeleton = vectorfield(vseir_skel_est_gamma_n_2),
           rinit = vseir_init_est_gamma_n_2,
           rmeasure = vseir_rmeas,
           dmeasure = vseir_dmeas,
           accumvars = sprintf("C_%d", 1:5),
           covar = covariate_table(covar, times = "Year", order = "constant"),
           covarnames = c(sprintf("N_%d", 1:5), sprintf("MU_%d", 1:5), 
                          sprintf("IN_%d", 1:5), sprintf("OUT_%d", 1:5), 
                          sprintf("p%d", 1:2), "Births", "eta_a"),
           statenames = state_names_est_gamma_n_2,
           paramnames = param_names_est, ...)
      ) 
    
  }
  
}





