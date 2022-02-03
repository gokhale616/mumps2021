# defining the model snippets for pomp -----------------------------------------------------------------------

##### determinstic skeleton - for simulation ##### 
vseir_skel_vacc_effective <- Csnippet("
  
    // define  the rate of waning of vaccine derived immunity
    double delta = 1/dwan; 

    int i,j,k,l;
    int n = agegroups; 
    
    double lambda1;
    double lambda_s_v;
    double lambda_w_v;
    
    const double *Cvlocal = &Cv1;
    const double *lvlocal = &lcv1;
    const double *qvlocal = &q_age_1;
    
    const double *Sslocal = &Ss_1;
    const double *Swlocal = Sslocal+n;
    const double *Stlocal = Swlocal+n;
    const double *Eslocal = Stlocal+n;
    const double *Ewlocal = Eslocal+n;
    const double *Etlocal = Ewlocal+n;
    const double *Islocal = Etlocal+n;
    const double *Iwlocal = Islocal+n;
    const double *Itlocal = Iwlocal+n;
    const double *Vlocal  = Itlocal+n;
    const double *Clocal  = Vlocal+n;
    const double *Vslocal = Clocal+n; 
    const double *lambdas_local = Vslocal+n;
    const double *lambdaw_local = lambdas_local+n;
    const double *Nvlocal = lambdaw_local+n; 
    const double *PRlocal = Nvlocal+n; 
    //const double *Hlocal  = PRlocal+n;
    
    
    double *DSslocal = &DSs_1;
    double *DSwlocal = DSslocal+n;
    double *DStlocal = DSwlocal+n;
    double *DEslocal = DStlocal+n;
    double *DEwlocal = DEslocal+n;
    double *DEtlocal = DEwlocal+n;
    double *DIslocal = DEtlocal+n;
    double *DIwlocal = DIslocal+n;
    double *DItlocal = DIwlocal+n;
    double *DVlocal  = DItlocal+n;
    double *DClocal  = DVlocal+n;
    double *DVslocal = DClocal+n; 
    double *Dlambdas_local = DVslocal+n;
    double *Dlambdaw_local = Dlambdas_local+n;
    double *DNvlocal = Dlambdaw_local+n; 
    double *DPRlocal = DNvlocal+n; 
    //double *DHlocal = DPRlocal+n;
    
    #define CM(J,K) Cvlocal[(J)+n*(K)]
    #define LV(K) lvlocal[(K)]
    #define QAGE(K) qvlocal[(K)]
    
    // Special Age-specific Covariates macro definitions  
    
    //#define N_AGE(J)   (__covars[__covindex[0]  + (J)])
    #define MU_AGE(J)  (__covars[__covindex[5]  + (J)])
    
    #define Ss(K) Sslocal[(K)]
    #define Sw(K) Swlocal[(K)]
    #define St(K) Stlocal[(K)]
    #define Es(K) Eslocal[(K)]
    #define Ew(K) Ewlocal[(K)]
    #define Et(K) Etlocal[(K)]
    #define Is(K) Islocal[(K)]
    #define Iw(K) Iwlocal[(K)]
    #define It(K) Itlocal[(K)]
    #define V(K) Vlocal[(K)]
    #define C(K) Clocal[(K)]
    #define Vs(K) Vslocal[(K)]
    #define Lambdas(K) lambdas_local[(K)]
    #define Lambdaw(K) lambdaw_local[(K)]
    #define Nv(K) Nvlocal[(K)]
    #define PR(K) PRlocal[(K)]
    //#define H(K)  Hlocal[(K)]
    
    
    #define DSs(K) DSslocal[(K)]
    #define DSw(K) DSwlocal[(K)]
    #define DSt(K) DStlocal[(K)]
    #define DEs(K) DEslocal[(K)]
    #define DEw(K) DEwlocal[(K)]
    #define DEt(K) DEtlocal[(K)]
    #define DIs(K) DIslocal[(K)]
    #define DIw(K) DIwlocal[(K)]
    #define DIt(K) DItlocal[(K)]
    #define DV(K)  DVlocal[(K)]
    #define DC(K)  DClocal[(K)]
    #define DVs(K) DVslocal[(K)]
    #define DLambdas(K) Dlambdas_local[(K)]
    #define DLambdaw(K) Dlambdaw_local[(K)]
    #define DNv(K) DNvlocal[(K)]
    #define DPR(K) DPRlocal[(K)]
    //#define DH(K)  DHlocal[(K)]
    

    
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
      
      if(i == 0) {
        
        // births occur only in the first age class - specially treat the [0,5) age cohort
        
        for (j = 0; j < n; j++) {
          lambda1  += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
          lambda_s_v += q*QAGE(i)*augCM[i][j]*(Is(j))/N_AGE(j);
          lambda_w_v += q*QAGE(i)*augCM[i][j]*(Iw(j))/N_AGE(j);
          }
            
        /* ============================ */
        /* balance transition equations */
        /* ============================ */
        
        // vaccinated and unvaccinated births - calculated using covariate population  
        
        double uv_births = (1-(1-alpha)*p1)*Births;
        double  v_births = (1-alpha)*p1*Births;
        
        DSs(i) = uv_births  - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);
        DSw(i) = delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);
        DSt(i) = uv_births + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*St(i);
        
        DEs(i) = lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i);   
        DEw(i) = lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i);
        DEt(i) = lambda1*St(i) - (sigma + LV(i) - MU_AGE(i))*Et(i); 
        
        
        DIs(i) = sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
        DIw(i) = sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
        DIt(i) = sigma*Et(i) - (gamma + LV(i) - MU_AGE(i))*It(i);
        
        DV(i)  = v_births - (delta + LV(i) - MU_AGE(i))*V(i);
        
        DVs(i) = delta*V(i); 
        
        DLambdas(i) = lambda_s_v; 
        DLambdaw(i) = lambda_w_v; 
      
      } else if (i == 1) {
        
          // individuals age in and out of age cohort after the first
        
          for (j = 0; j < n; j++) {
            lambda1  += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
            lambda_s_v += q*QAGE(i)*augCM[i][j]*(Is(j))/N_AGE(j);
            lambda_w_v += q*QAGE(i)*augCM[i][j]*(Iw(j))/N_AGE(j);
          }
            
          /* ============================ */
          /* balance tranistion equations */
          /* ============================ */
          
          //double uv_grads_Ss = LV(i-1)*(1-(1-alpha)*p2)*Ss(i-1);
          //double uv_grads_Sw = LV(i-1)*(1-(1-alpha)*p2)*Sw(i-1);
          double    uv_grads = LV(i-1)*(1-(1-alpha)*p2)*St(i-1);
          double     v_grads = LV(i-1)*(1-alpha)*p2*St(i-1);
          
          DSs(i) =  uv_grads  - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);               
          DSw(i) =  delta*V(i) + LV(i-1)*Sw(i-1) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);               
          DSt(i) =  uv_grads + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*St(i); 
          
          DEs(i) = LV(i-1)*Es(i-1) + lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i); 
          DEw(i) = LV(i-1)*Ew(i-1) + lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i); 
          DEt(i) = LV(i-1)*Et(i-1) + lambda1*St(i) - (sigma + LV(i) - MU_AGE(i))*Et(i); 
          
          DIs(i) = LV(i-1)*Is(i-1) + sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
          DIw(i) = LV(i-1)*Iw(i-1) + sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
          DIt(i) = LV(i-1)*It(i-1) + sigma*Et(i) - (gamma + LV(i) - MU_AGE(i))*It(i);
          
          DV(i)  = v_grads + LV(i-1)*V(i-1) - (delta + LV(i) - MU_AGE(i))*V(i);
          
          DVs(i) = delta*V(i);
          
          DLambdas(i) = lambda_s_v; 
          DLambdaw(i) = lambda_w_v; 
      
      } else {
        
          // individuals age in and out of age cohort after the first
          for (j = 0; j < n; j++) {
            lambda1  += q*QAGE(i)*augCM[i][j]*(It(j))/N_AGE(j);
            lambda_s_v += q*QAGE(i)*augCM[i][j]*(Is(j))/N_AGE(j);
            lambda_w_v += q*QAGE(i)*augCM[i][j]*(Iw(j))/N_AGE(j);
          }
          
          /* ============================ */
          /* balance transition equations */
          /* ============================ */
          
          DSs(i) = LV(i-1)*Ss(i-1) - (lambda1 + LV(i) - MU_AGE(i))*Ss(i);
          DSw(i) = LV(i-1)*Sw(i-1) + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*Sw(i);
          DSt(i) = LV(i-1)*St(i-1) + delta*V(i) - (lambda1 + LV(i) - MU_AGE(i))*St(i);
          
          DEs(i) = LV(i-1)*Es(i-1) + lambda1*Ss(i) - (sigma + LV(i) - MU_AGE(i))*Es(i); 
          DEw(i) = LV(i-1)*Ew(i-1) + lambda1*Sw(i) - (sigma + LV(i) - MU_AGE(i))*Ew(i); 
          DEt(i) = LV(i-1)*Et(i-1) + lambda1*St(i) - (sigma + LV(i) - MU_AGE(i))*Et(i);  
          
          DIs(i) = LV(i-1)*Is(i-1) + sigma*Es(i) - (gamma + LV(i) - MU_AGE(i))*Is(i);
          DIw(i) = LV(i-1)*Iw(i-1) + sigma*Ew(i) - (gamma + LV(i) - MU_AGE(i))*Iw(i);
          DIt(i) = LV(i-1)*It(i-1) + sigma*Et(i) - (gamma + LV(i) - MU_AGE(i))*It(i);
          
          DV(i)  = LV(i-1)*V(i-1)  - (delta + LV(i) - MU_AGE(i))*V(i);
          
          DVs(i) = delta*V(i);
          
          DLambdas(i) = lambda_s_v; 
          DLambdaw(i) = lambda_w_v; 
          
          
          
      }
      
      /* ========================= */
      /* Keep account of new cases */
      /* ========================= */
      DC(i)  = gamma*It(i); 
      DNv(i) = N_AGE(i); 
      DPR(i) = Iw(i)/Is(i); 
      //DH(i)  = Lambdaw(i)/Lambdas(i);
      
    }
")



##### initial value distibution - for simulation #####
vseir_init_vacc_effective <- Csnippet("
  double *Ss = &Ss_1;
  double *Sw = &Sw_1;
  double *St = &St_1;
  double *Es = &Es_1;
  double *Ew = &Ew_1;
  double *Et = &Et_1;
  double *Is = &Is_1;
  double *Iw = &Iw_1;
  double *It = &It_1;
  
  double *V = &V_1;
  double *C = &C_1;
  double *Vs = &Vs_1;
  double *lambdas = &lambdas_1;
  double *lambdaw = &lambdaw_1;
  double *Nv = &Nv_1;
  double *PR = &PR_1;
  //double *H = &H_1;
  
  int n = agegroups;
  
  #define N_AGE(J)  (__covars[__covindex[0]  + (J)])
  
  for (int i = 0; i < n; i++) {
    Ss[i]  = nearbyint(S_0*N_AGE(i));
    Sw[i]  = 0;
    St[i]  = Ss[i] + Sw[i];
    
    Es[i] = 0;
    Ew[i] = 0;
    Et[i] = Es[i] + Ew[i];
    
    Is[i] = nearbyint(I_0*N_AGE(i));
    Iw[i] = 0;
    It[i] = Is[i] + Iw[i];
    
    V[i]  = 0;
    C[i]  = 0;
    Vs[i] = 0;
    
    lambdas[i] = 0; 
    lambdaw[i] = 0; 
    
    Nv[i] = 0; 
    PR[i] = 0;
   //H[i] = 0;
  }
")


##### parameter names ######
state_names_vacc_effective <- (
  c(sprintf("Ss_%d",1:5), sprintf("Sw_%d",1:5), sprintf("St_%d",1:5), 
    sprintf("Es_%d",1:5), sprintf("Ew_%d",1:5), sprintf("Et_%d",1:5), 
    sprintf("Is_%d",1:5), sprintf("Iw_%d",1:5), sprintf("It_%d",1:5),
    sprintf("V_%d",1:5), sprintf("C_%d", 1:5),
    sprintf("Vs_%d",1:5), 
    sprintf("lambdas_%d", 1:5), sprintf("lambdaw_%d", 1:5),
    sprintf("Nv_%d",1:5), 
    sprintf("PR_%d",1:5)#, 
    #sprintf("H_%d",1:5)
    )
  )


param_names_vacc_effective <- (
  c("S_0", "I_0",
    "q", sprintf("q_age_%d", 1:5), 
     "gamma", "sigma", "beta1",
     "rho",
     sprintf("rho_age_%d", 1:5), "rho_age_u", 
     "psi",
     sprintf("psi_%d", 1:5), "psi_u", 
     sprintf("Cv%d",1:25),
     sprintf("lcv%d", 1:5), 
     "alpha", "dwan", "agegroups")
  )



# testing parameter values 
##### setting up parameter values ###### setting up default q such that R0 is 10
param_vals_vacc_effective <- (
  c(S_0 = 1e-1, I_0 = 1e-4,
    q = 1, 
    q_age_1 = 1, q_age_2 = 1, q_age_3 = 1, q_age_4 = 1, q_age_5 = 1, 
    gamma = 365.25/5, sigma = 365.25/17, 
    dwan = Inf, 
    rho = 1,
    rho_age_1 = 1, rho_age_2 = 1, rho_age_3 = 1, rho_age_4 = 1, rho_age_5 = 1, rho_age_u = 1, 
    psi = 1, 
    psi_1 = 1, psi_2 = 1, psi_3 = 1, psi_4 = 1, psi_5 = 1, psi_u = 1,
    agegroups = 5, Cv = C, alpha = 0.054, beta1 = 0.11,
    lcv = c(1/5, 1/10, 1/10, 1/15, 1/40))
  )


##### build the pomp object #####
# Here the model is defined at an annual scale with a time step of that of a day (dt = 1/365.25)
# When extrapolating the simulation, default scale of simulation has been set to annual (temp_scale = 1)
# This can of course be changed as required by passing arguments

make_pomp_vacc_effective <- function(..., 
                                     start_t = 1800, 
                                     extrapolate_simulation = TRUE, 
                                     extra_start_t = 1967, extra_end_t = 2018, temp_scale = 1, 
                                     covar) {
  
  
  message("Assuming exponentially distributed waning rate.")
  
  
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
           skeleton = vectorfield(vseir_skel_vacc_effective),
           rinit = vseir_init_vacc_effective,
           accumvars = c(sprintf("C_%d", 1:5), sprintf("Vs_%d", 1:5),
                         sprintf("lambdas_%d", 1:5), sprintf("lambdaw_%d", 1:5),
                         sprintf("Nv_%d", 1:5), sprintf("PR_%d",1:5)# , 
                         ),  
           covar = covariate_table(covar, times = "Year", order = "constant"),
           covarnames = c(sprintf("N_%d", 1:5), sprintf("MU_%d", 1:5), 
                          sprintf("IN_%d", 1:5), sprintf("OUT_%d", 1:5), 
                          sprintf("p%d", 1:2), "Births", "eta_a"),
           statenames = state_names_vacc_effective,
           paramnames = param_names_vacc_effective, 
           #cdir = getwd(), 
           #cfile = "test_vacc_effectiveness",
           ...)
      ) 
  
  } else {
    
    message("This function has no estimation protocol: Please use 'extrapolate_simulation = TRUE'")
  }
  
  }


