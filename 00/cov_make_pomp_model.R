# contact matrix - set the scale to annual mean contacts
C <- contact_matrix*365.25


# defining the model snippets for pomp -----------------------------------------------------------------------

##### stochastic process: demographic and observation #### 
vseir_step <- Csnippet("
  // define  the 
  double delta = 1/dwan; 
  
  int i,j,k,l;
  int n = agegroups; 
  
  double lambda1;
  double lambda2;
  
  const double *Cvlocal = &Cv1;
  const double *lvlocal = &lcv1;
  const double *qvlocal = &q_age_1;
  
  double *Slocal  = &S_1;
  double *E1local = Slocal+n;
  double *I1local = E1local+n;
  double *E2local = I1local+n;
  double *I2local = E2local+n;
  double *Rlocal  = I2local+n;
  double *Vlocal  = Rlocal+n;
  double *Clocal  = Vlocal+n;
  double *Alocal  = Clocal+n;
  double *Glocal  = Alocal+n;
  
  #define CM(J,K) Cvlocal[(J)+n*(K)]
  #define LV(K) lvlocal[(K)]
  #define QAGE(K) qvlocal[(K)]
  
  // Special Age-specific Covariates macro definitions  
  #define IN_AGE(J)  (__covars[__covindex[11]+(J)])
  #define OUT_AGE(J) (__covars[__covindex[16]+(J)])
  
  
  #define S(K) Slocal[(K)]
  #define E1(K) E1local[(K)]
  #define I1(K) I1local[(K)]
  #define E2(K) E2local[(K)]
  #define I2(K) I2local[(K)]
  #define R(K) Rlocal[(K)]
  #define V(K) Vlocal[(K)]
  #define C(K) Clocal[(K)]
  #define A(K) Alocal[(K)]
  #define G(K) Glocal[(K)]
  
  double r_S1[6];
  
  double r_S[5];
  double r_E1[4];
  double r_I1[4];
  double r_E2[4];
  double r_I2[4];
  double r_R[3];
  double r_V[6];
  
  double dN_S1[6];
  
  double dN_S[5];
  double dN_E1[4];
  double dN_I1[4];
  double dN_E2[4];
  double dN_I2[4];
  double dN_R[3];
  double dN_V[6];
  
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

  /* ======================================================= */
  /* define and initialize matrices to store new transitions */ 
  /* ======================================================= */
  
  double fromS[n][5];
  double fromE1[n][4];
  double fromI1[n][4];
  double fromE2[n][4];
  double fromI2[n][4];
  double fromR[n][3];
  double fromV[n][6];
  
  for(i = 0; i < n; i++) {
    //susceptible has 3 flux rates
    for(j=0; j<5; j++) {
      fromS[i][j] = 0;
    }
    
    // all compartments with 2 flux rates
    for(j=0; j<4; j++) {
      fromE1[i][j] = 0;
      fromI1[i][j] = 0;
      fromE2[i][j] = 0;
      fromI2[i][j] = 0;
    }  
    
    // recovered has 1 flux rate  
    for(j=0; j<3; j++) {
      fromR[i][0] = 0;
    }
  
    // vaccinated has 4 flux rate
    for(j = 0; j<6; j++) {
      fromV[i][j] = 0;
    }
    
  }
  
  double uv_births = rpois(Births*(1-(1-alpha)*p1)*dt);
  double  v_births = rpois(Births*(1-alpha)*p1*dt);
  
  
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
            lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j) + iota_v[i])/N_AGE(j);
        }
      }

      // CAUTION :: The following array is to account for the booster coverage 
      
      r_S1[0] = lambda1; r_S1[1] = lambda2; r_S1[2] = (1-(1-alpha)*p2)*LV(i); r_S1[3] = (1-alpha)*p2*LV(i); 
      r_S1[4] = IN_AGE(i); r_S1[5] = OUT_AGE(i); 
      
      r_E1[0] = sigma;  r_E1[1] = LV(i);  r_E1[2] = IN_AGE(i);  r_E1[3] = OUT_AGE(i);
      r_I1[0] = gamma;  r_I1[1] = LV(i);  r_I1[2] = IN_AGE(i);  r_I1[3] = OUT_AGE(i);
      r_E2[0] = sigma;  r_E2[1] = LV(i);  r_E2[2] = IN_AGE(i);  r_E2[3] = OUT_AGE(i);
      r_I2[0] = gamma;  r_I2[1] = LV(i);  r_I2[2] = IN_AGE(i);  r_I2[3] = OUT_AGE(i);
      
      r_R[0]  = LV(i);  r_R[1] = IN_AGE(i); r_R[2] = OUT_AGE(i); 
      
      r_V[0] = epsilon1*lambda1; r_V[1] = epsilon2*lambda2; r_V[2] = delta; r_V[3] = LV(i);
      r_V[4] = IN_AGE(i); r_V[5] = OUT_AGE(i);
      
      // make the general susceptible arrays for this age class zero - because of the vaccination
      //for(k = 0; k < 6; k++) {fromS[i][k] = 0; dN_S[k] = 0; r_S[k] = 0;}
      
      // draws from the susceptible compartment - to accomodate vaccine booster
      reulermultinom(6, S(i), &r_S1[0], dt, &dN_S1[0]);
      
      // draws from the exposed1 class
      reulermultinom(4, E1(i), &r_E1[0], dt, &dN_E1[0]); 
      for(k = 0; k<4; k++) {fromE1[i][k] = dN_E1[k];}
      
      // draws from the infectious1 class
      reulermultinom(4, I1(i), &r_I1[0], dt, &dN_I1[0]);
      for(k=0; k<4; k++) {fromI1[i][k] = dN_I1[k];}
      
      // draws from the exposed2 class
      reulermultinom(4, E2(i), &r_E2[0], dt, &dN_E2[0]);
      for(k=0; k<4; k++) {fromE2[i][k] = dN_E2[k];}
      
      //draws from the infectious2 class
      reulermultinom(4, I2(i), &r_I2[0], dt, &dN_I2[0]);
      for(k=0; k<4; k++) {fromI2[i][k] = dN_I2[k];}
      
      // draws from the recovered class
      reulermultinom(3, R(i), &r_R[0], dt, &dN_R[0]);
      for(k=0; k<3; k++) {fromR[i][k] = dN_R[k];}
      
      // draws from the vaccinated class
      reulermultinom(6, V(i), &r_V[0], dt, &dN_V[0]);
      for(k=0; k<6; k++) {fromV[i][k] = dN_V[k];}
      
      /* ============================ */
      /* balance transition equations */
      /* ============================ */
      
      S(i)  += uv_births + fromV[i][2] - dN_S1[0] - dN_S1[1] - dN_S1[2] - dN_S1[3] + dN_S1[4] - dN_S1[5]; 
    
      E1(i) += dN_S1[0] + fromV[i][0]  - fromE1[i][0] - fromE1[i][1] + fromE1[i][2] - fromE1[i][3];
      
      I1(i) += fromE1[i][0] - fromI1[i][0] - fromI1[i][1] + fromI1[i][2] - fromI1[i][3];
      
      E2(i) += dN_S1[1] + fromV[i][1]  - fromE2[i][0] - fromE2[i][1] + fromE2[i][2] - fromE2[i][3];
      
      I2(i) += fromE2[i][0] - fromI2[i][0] - fromI2[i][1] + fromI2[i][2] - fromI2[i][3];                                    
      
      R(i)  += fromI1[i][0] + fromI2[i][0] - fromR[i][0] + fromR[i][1] - fromR[i][2];   
      
      V(i)  += v_births - fromV[i][0] - fromV[i][1] - fromV[i][2] - fromV[i][3] + fromV[i][4] - fromV[i][5];  
    
    } else if (i == 1) {
      
        // individuals age in and out of age cohort after the first
        
        for (j = 0; j < n; j++) {
          lambda1 += q*QAGE(i)*augCM[i][j]*(I1(j))/N_AGE(j);
          if (t < t_intro) {
            lambda2 = 0;
          } else {
              lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j) + iota_v[i])/N_AGE(j);
          }
        }
      
        r_S[0] = lambda1; r_S[1] = lambda2; r_S[2] = LV(i); r_S[3] = IN_AGE(i); r_S[4] = OUT_AGE(i);
        
        r_E1[0] = sigma; r_E1[1] = LV(i); r_E1[2] = IN_AGE(i);  r_E1[3] = OUT_AGE(i);
        r_I1[0] = gamma; r_I1[1] = LV(i); r_I1[2] = IN_AGE(i);  r_I1[3] = OUT_AGE(i);
        r_E2[0] = sigma; r_E2[1] = LV(i); r_E2[2] = IN_AGE(i);  r_E2[3] = OUT_AGE(i);
        r_I2[0] = gamma; r_I2[1] = LV(i); r_I2[2] = IN_AGE(i);  r_I2[3] = OUT_AGE(i);
        
        r_R[0] = LV(i); r_R[1] = IN_AGE(i); r_R[2] = OUT_AGE(i); 
        
        r_V[0] = epsilon1*lambda1; r_V[1] = epsilon2*lambda2; r_V[2] = delta; r_V[3] = LV(i);
        r_V[4] = IN_AGE(i); r_V[5] = OUT_AGE(i);
        
        // draws from the susceptible classes
        reulermultinom(5, S(i), &r_S[0], dt, &dN_S[0]); 
        for(k=0; k<5; k++) {fromS[i][k] = dN_S[k];}
        
        // draws from the exposed1 classes
        reulermultinom(4, E1(i), &r_E1[0], dt, &dN_E1[0]); 
        for(k=0; k<4; k++) {fromE1[i][k] = dN_E1[k];}
        
        // draws from the infectious1 classes
        reulermultinom(4, I1(i), &r_I1[0], dt, &dN_I1[0]);
        for(k=0; k<4; k++) {fromI1[i][k] = dN_I1[k];}
        
        //draws from the exposed2 classes
        reulermultinom(4, E2(i), &r_E2[0], dt, &dN_E2[0]);
        for(k=0; k<4; k++) {fromE2[i][k] = dN_E2[k];}
        
        //draws from the infectious2 classes
        reulermultinom(4, I2(i), &r_I2[0], dt, &dN_I2[0]);
        for(k=0; k<4; k++) {fromI2[i][k] = dN_I2[k];}
        
        // draws from the recovered classes
        reulermultinom(3, R(i), &r_R[0], dt, &dN_R[0]);
        for(k=0; k<3; k++) {fromR[i][k] = dN_R[k];}
        
        // draws from the vaccinated class
        reulermultinom(6, V(i), &r_V[0], dt, &dN_V[0]);
        for(k=0; k<6; k++) {fromV[i][k] = dN_V[k];}
        

        /* ============================ */
        /* balance tranistion equations */
        /* ============================ */
        
        S(i)  += dN_S1[2] + fromV[i][2] - fromS[i][0] - fromS[i][1] - fromS[i][2] + fromS[i][3] - fromS[i][4];               
        
        E1(i) += fromE1[i-1][1] + fromS[i][0]  + fromV[i][0]  - fromE1[i][0] - 
                   fromE1[i][1] + fromE1[i][2] - fromE1[i][3]; 
        I1(i) += fromI1[i-1][1] + fromE1[i][0] - fromI1[i][0] - fromI1[i][1] + fromI1[i][2] - fromI1[i][3];               
        
        E2(i) += fromE2[i-1][1] + fromS[i][1]  + fromV[i][1]  - fromE2[i][0] - 
                   fromE2[i][1] + fromE2[i][2] - fromE2[i][3];
        I2(i) += fromI2[i-1][1] + fromE2[i][0] - fromI2[i][0] - fromI2[i][1] + fromI2[i][2] - fromI2[i][3];
        
        R(i)  += fromR[i-1][0]  + fromI1[i][0] + fromI2[i][0] - fromR[i][0] + fromR[i][1] - fromR[i][2];                             
        
        V(i)  += fromV[i-1][3]  + dN_S1[3] - fromV[i][0] - fromV[i][1] - 
                 fromV[i][2] - fromV[i][3] + fromV[i][4] - fromV[i][5];                
    
    } else {
      
      // individuals age in and out of age cohort after the first
        
        for (j = 0; j < n; j++) {
          lambda1 += q*QAGE(i)*augCM[i][j]*(I1(j))/N_AGE(j);
          if (t < t_intro) {
            lambda2 = 0;
          } else {
              lambda2 += q*QAGE(i)*augCM[i][j]*(I2(j) + iota_v[i])/N_AGE(j);
          }
        }
        
        r_S[0]  = lambda1; r_S[1]  = lambda2; r_S[2] = LV(i); r_S[3] = IN_AGE(i); r_S[4] = OUT_AGE(i);
        r_E1[0] = sigma;   r_E1[1] = LV(i); r_E1[2] = IN_AGE(i);  r_E1[3] = OUT_AGE(i);
        r_I1[0] = gamma;   r_I1[1] = LV(i); r_I1[2] = IN_AGE(i);  r_I1[3] = OUT_AGE(i);
        r_E2[0] = sigma;   r_E2[1] = LV(i); r_E2[2] = IN_AGE(i);  r_E2[3] = OUT_AGE(i);
        r_I2[0] = gamma;   r_I2[1] = LV(i); r_I2[2] = IN_AGE(i);  r_I2[3] = OUT_AGE(i);
        
        r_R[0]  = LV(i); r_R[1] = IN_AGE(i); r_R[2] = OUT_AGE(i);   
        
        r_V[0]  = epsilon1*lambda1; r_V[1] = epsilon2*lambda2; r_V[2] = delta; r_V[3] = LV(i);
        r_V[4] = IN_AGE(i); r_V[5] = OUT_AGE(i);
        
        // draws from the susceptible classes
        reulermultinom(5, S(i), &r_S[0], dt, &dN_S[0]); 
        for(k=0; k<5; k++) {fromS[i][k] = dN_S[k];}
        
        // draws from the exposed1 classes
        reulermultinom(4, E1(i), &r_E1[0], dt, &dN_E1[0]); 
        for(k=0; k<4; k++) {fromE1[i][k] = dN_E1[k];}
        
        // draws from the infectious1 classes
        reulermultinom(4, I1(i), &r_I1[0], dt, &dN_I1[0]);
        for(k=0; k<4; k++) {fromI1[i][k] = dN_I1[k];}
        
        //draws from the exposed2 classes
        reulermultinom(4,E2(i), &r_I2[0], dt, &dN_E2[0]);
        for(k=0; k<4; k++) {fromE2[i][k] = dN_E2[k];}
        
        // draws from the infectious2 classes
        reulermultinom(4, I2(i), &r_I2[0], dt, &dN_I2[0]);
        for(k=0; k<4; k++) {fromI2[i][k] = dN_I2[k];}
        
        // draws from the recovered classes
        reulermultinom(3, R(i), &r_R[0], dt, &dN_R[0]);
        for(k=0; k<3; k++) {fromR[i][k] = dN_R[k];}
        
        // draws from the vaccinated class
        reulermultinom(6, V(i), &r_V[0], dt, &dN_V[0]);
        for(k=0; k<6; k++) {fromV[i][k] = dN_V[k];}
      
        
        /* ============================ */
        /* balance tranistion equations */
        /* ============================ */
        
        S(i)  += fromS[i-1][2] + fromV[i][2]  - fromS[i][0]  - fromS[i][1] - 
                   fromS[i][2] + fromS[i][2] - fromS[i][3];               
        
        E1(i) += fromE1[i-1][1] + fromS[i][0]  + fromV[i][0]  - fromE1[i][0] - 
                   fromE1[i][1] + fromE1[i][2] - fromE1[i][3]; 
        I1(i) += fromI1[i-1][1] + fromE1[i][0] - fromI1[i][0] - fromI1[i][1] + fromI1[i][2] - fromI1[i][3];    
        
        E2(i) += fromE2[i-1][1] + fromS[i][1]  + fromV[i][1]  - fromE2[i][0] - 
                   fromE2[i][1] + fromE2[i][2] - fromE2[i][3]; 
        I2(i) += fromI2[i-1][1] + fromE2[i][0] - fromI2[i][0] - fromI2[i][1] + fromI2[i][2] - fromI2[i][3];               
        
        R(i)  += fromR[i-1][0]  + fromI1[i][0] + fromI2[i][0] - fromR[i][0] + fromR[i][1] - fromR[i][2];                             
        
        V(i)  += fromV[i-1][3] - fromV[i][0] - fromV[i][1]  - fromV[i][2]  - 
                   fromV[i][3] + fromV[i][4] - fromV[i][5];               
        
    }
    
    /* ========================= */
    /* Keep account of new cases */
    /* ========================= */
    A(i) += fromI1[i][0];
    G(i) += fromI2[i][0];
    
    C(i) += (fromI1[i][0] + fromI2[i][0]); 
  }

")

##### deterministic skeleton - for estimation #####
vseir_skel_est <- Csnippet("
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
    const double *Vlocal  = I2local+n;;
    const double *Clocal  = Vlocal+n;
    
    double *DSlocal  = &DS_1;
    double *DE1local = DSlocal+n;
    double *DI1local = DE1local+n;
    double *DE2local = DI1local+n;
    double *DI2local = DE2local+n;
    double *DVlocal  = DI2local+n;;
    double *DClocal  = DVlocal+n;
    
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
    #define V(K) Vlocal[(K)]
    #define C(K) Clocal[(K)]
    
    #define DS(K) DSlocal[(K)]
    #define DE1(K) DE1local[(K)]
    #define DI1(K) DI1local[(K)]
    #define DE2(K) DE2local[(K)]
    #define DI2(K) DI2local[(K)]
    #define DV(K) DVlocal[(K)]
    #define DC(K) DClocal[(K)]
    

    double r_S1[4];
    
    double r_S[3];
    double r_E1[2];
    double r_I1[2];
    double r_E2[2];
    double r_I2[2];
    double r_V[4];
    
    double dN_S1[4];
    
    double dN_S[3];
    double dN_E1[2];
    double dN_I1[2];
    double dN_E2[2];
    double dN_I2[2];
    double dN_V[4];
    
    
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
        
        DS(i)  = uv_births + delta*V(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);
        DE1(i) = lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i);   
        DI1(i) = sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
        DE2(i) = lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
        DI2(i) = sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
        DV(i)  = v_births - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
      
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
          
          DS(i)  = uv_grads + delta*V(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);               
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          DV(i)  = v_grads + LV(i-1)*V(i-1) - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
      
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
          
          DS(i)  = LV(i-1)*S(i-1)  + delta*V(i)    - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);             
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          DV(i)  = LV(i-1)*V(i-1)  - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
          
      }
      
      /* ========================= */
      /* Keep account of new cases */
      /* ========================= */
      DC(i) = gamma*(I1(i)+I2(i)); 
    }
")



##### determinstic skeleton - for simulation ##### 
vseir_skel_sim <- Csnippet("
  
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
    const double *Rlocal  = I2local+n;
    const double *Vlocal  = Rlocal+n;
    const double *Clocal  = Vlocal+n;
    const double *Alocal  = Clocal+n;
    const double *Glocal  = Alocal+n;
    
    double *DSlocal  = &DS_1;
    double *DE1local = DSlocal+n;
    double *DI1local = DE1local+n;
    double *DE2local = DI1local+n;
    double *DI2local = DE2local+n;
    double *DRlocal  = DI2local+n;
    double *DVlocal  = DRlocal+n;
    double *DClocal  = DVlocal+n;
    double *DAlocal  = DClocal+n;
    double *DGlocal  = DAlocal+n;
    
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
    #define R(K) Rlocal[(K)]
    #define V(K) Vlocal[(K)]
    #define C(K) Clocal[(K)]
    #define A(K) Alocal[(K)]
    #define G(K) Glocal[(K)]
    
    #define DS(K) DSlocal[(K)]
    #define DE1(K) DE1local[(K)]
    #define DI1(K) DI1local[(K)]
    #define DE2(K) DE2local[(K)]
    #define DI2(K) DI2local[(K)]
    #define DR(K) DRlocal[(K)]
    #define DV(K) DVlocal[(K)]
    #define DC(K) DClocal[(K)]
    #define DA(K) DAlocal[(K)]
    #define DG(K) DGlocal[(K)]
    
    double r_S1[4];
    
    double r_S[3];
    double r_E1[2];
    double r_I1[2];
    double r_E2[2];
    double r_I2[2];
    double r_R[1];
    double r_V[4];
    
    double dN_S1[4];
    
    double dN_S[3];
    double dN_E1[2];
    double dN_I1[2];
    double dN_E2[2];
    double dN_I2[2];
    double dN_R[1];
    double dN_V[4];
    
    
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
        
        DS(i)  = uv_births + delta*V(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);
        DE1(i) = lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i);   
        DI1(i) = sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
        DE2(i) = lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
        DI2(i) = sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
        DR(i)  = gamma*(I1(i) + I2(i)) - (LV(i) - MU_AGE(i))*R(i);
        DV(i)  = v_births - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
      
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
          
          DS(i)  = uv_grads + delta*V(i) - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);               
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          DR(i)  = LV(i-1)*R(i-1)  + gamma*(I1(i) + I2(i)) - (LV(i) - MU_AGE(i))*R(i); 
          DV(i)  = v_grads + LV(i-1)*V(i-1) - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
      
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
          
          DS(i)  = LV(i-1)*S(i-1)  + delta*V(i)    - (lambda1 + lambda2 + LV(i) - MU_AGE(i))*S(i);             
          DE1(i) = LV(i-1)*E1(i-1) + lambda1*(S(i) + epsilon1*V(i)) - (sigma + LV(i) - MU_AGE(i))*E1(i); 
          DI1(i) = LV(i-1)*I1(i-1) + sigma*E1(i) - (gamma + LV(i) - MU_AGE(i))*I1(i);
          DE2(i) = LV(i-1)*E2(i-1) + lambda2*(S(i) + epsilon2*V(i)) - (sigma + LV(i) - MU_AGE(i))*E2(i);
          DI2(i) = LV(i-1)*I2(i-1) + sigma*E2(i) - (gamma + LV(i) - MU_AGE(i))*I2(i) + iota_v[i];
          DR(i)  = LV(i-1)*R(i-1)  + gamma*(I1(i) + I2(i)) - (LV(i) - MU_AGE(i))*R(i);
          DV(i)  = LV(i-1)*V(i-1)  - (epsilon1*lambda1 + epsilon2*lambda2 + delta + LV(i) - MU_AGE(i))*V(i);
          
      }
      
      /* ========================= */
      /* Keep account of new cases */
      /* ========================= */
      DA(i) = gamma*I1(i);
      DG(i) = gamma*I2(i);
      
      DC(i) = gamma*(I1(i)+I2(i)); 
    }
")


##### initial value distibution - for estimation #####
vseir_init_est <- Csnippet("
  double *S = &S_1;
  double *E1 = &E1_1;
  double *I1 = &I1_1;
  double *E2 = &E2_1;
  double *I2 = &I2_1;
  double *V = &V_1;
  double *C = &C_1;
  
  int n = agegroups;
  
  #define N_AGE(J)  (__covars[__covindex[0]  + (J)])
  
  for (int i = 0; i < n; i++) {
    S[i]  = nearbyint(S_0*N_AGE(i));
    E1[i] = 0;
    I1[i] = nearbyint(I1_0*N_AGE(i));
    E2[i] = 0;
    I2[i] = nearbyint(I2_0*N_AGE(i));
    V[i]  = 0;
    C[i]  = 0;
  }
")



##### initial value distibution - for simulation #####
vseir_init_sim <- Csnippet("
  double *S = &S_1;
  double *E1 = &E1_1;
  double *I1 = &I1_1;
  double *E2 = &E2_1;
  double *I2 = &I2_1;
  double *R = &R_1;
  double *V = &V_1;
  double *C = &C_1;
  double *A = &A_1;
  double *G = &G_1;
  
  int n = agegroups;
  
  #define N_AGE(J)  (__covars[__covindex[0]  + (J)])
  
  for (int i = 0; i < n; i++) {
    S[i]  = nearbyint(S_0*N_AGE(i));
    E1[i] = 0;
    I1[i] = nearbyint(I1_0*N_AGE(i));
    E2[i] = 0;
    I2[i] = nearbyint(I2_0*N_AGE(i));;
    R[i]  = nearbyint(R_0*N_AGE(i));
    V[i]  = 0;
    C[i]  = 0;
    A[i]  = 0;
    G[i]  = 0;
  }

                       
                       
")

##### measurement model #####

vseir_rmeas <- Csnippet("
  
  //make sure that no soln values used to simulate from the measurement model are negative
  
  double C_1_tmp, C_2_tmp, C_3_tmp, C_4_tmp, C_5_tmp;
  
  if (C_1 < 0) {
    C_1_tmp = 0;
  } else {
    C_1_tmp = C_1;
  }
  
  if (C_2 < 0) {
    C_2_tmp = 0;
  } else {
    C_2_tmp = C_2;
  }
  
  if (C_3 < 0) {
    C_3_tmp = 0;
  } else {
    C_3_tmp = C_3;
  }
  
  if (C_4 < 0) {
    C_4_tmp = 0;
  } else {
    C_4_tmp = C_4;
  }
  
  if (C_5 < 0) {
    C_5_tmp = 0;
  } else {
    C_5_tmp = C_5;
  }
  
  // define scalar multiples to take care of all types of reporting
  double s_1 = rho*rho_age_1*eta_a;
  double s_2 = rho*rho_age_2*eta_a;
  double s_3 = rho*rho_age_3*eta_a;
  double s_4 = rho*rho_age_4*eta_a;
  double s_5 = rho*rho_age_5*eta_a;
  double s_u = rho*rho_age_u*(1-eta_a);
  
  // define scalar multiples to take care of age specific/independent over-dispersion
  double p_sq_1 = (psi*psi)*(psi_1*psi_1);
  double p_sq_2 = (psi*psi)*(psi_2*psi_2);
  double p_sq_3 = (psi*psi)*(psi_3*psi_3);
  double p_sq_4 = (psi*psi)*(psi_4*psi_4);
  double p_sq_5 = (psi*psi)*(psi_5*psi_5);
  double p_sq_u = (psi*psi)*(psi_u*psi_u);
  
  // define the age-specific mean of the reporting distribution
  double m_1 = s_1*C_1_tmp; 
  double m_2 = s_2*C_2_tmp; 
  double m_3 = s_3*C_3_tmp; 
  double m_4 = s_4*C_4_tmp; 
  double m_5 = s_5*C_5_tmp; 
  double m_u = s_u*(C_1_tmp + C_2_tmp + C_3_tmp + C_4_tmp + C_5_tmp);
  
  // define the age-specific variance of the reporting distribution
  double v_1 = m_1*(1 - s_1 + p_sq_1*m_1);
  double v_2 = m_2*(1 - s_2 + p_sq_2*m_2);
  double v_3 = m_3*(1 - s_3 + p_sq_3*m_3);
  double v_4 = m_4*(1 - s_4 + p_sq_4*m_4);
  double v_5 = m_5*(1 - s_5 + p_sq_5*m_5);
  double v_u = m_u*(1 - s_u + p_sq_u*m_u);
  
  if(ISNA(eta_a)) {
    
    D_1  = NA_REAL;
    D_2  = NA_REAL;
    D_3  = NA_REAL;
    D_4  = NA_REAL;
    D_5  = NA_REAL;
    D_u  = NA_REAL;
  
  } else {
      
      // a bit more book-keeeping 
      double tol = 1.0e-18; 
      
      double D_1_tmp = rnorm(m_1, sqrt(v_1)+tol); 
      double D_2_tmp = rnorm(m_2, sqrt(v_2)+tol); 
      double D_3_tmp = rnorm(m_3, sqrt(v_3)+tol); 
      double D_4_tmp = rnorm(m_4, sqrt(v_4)+tol); 
      double D_5_tmp = rnorm(m_5, sqrt(v_5)+tol); 
      double D_u_tmp = rnorm(m_u, sqrt(v_u)+tol); 
      
      if (D_1_tmp < 0) {
        D_1 = 0; 
      } else {
        D_1 = nearbyint(D_1_tmp);
      }
      
      if (D_2_tmp < 0) {
        D_2 = 0; 
      } else {
        D_2 = nearbyint(D_2_tmp);
      }
      
      if (D_3_tmp < 0) {
        D_3 = 0; 
      } else {
        D_3 = nearbyint(D_3_tmp);
      }
      
      if (D_4_tmp < 0) {
        D_4 = 0; 
      } else {
        D_4 = nearbyint(D_4_tmp);
      }
      
      if (D_5_tmp < 0) {
        D_5 = 0; 
      } else {
        D_5 = nearbyint(D_5_tmp);
      }
      
      if (D_u_tmp < 0) {
        D_u = 0; 
      } else {
        D_u = nearbyint(D_u_tmp);
      }
      
  }
")


vseir_dmeas <- Csnippet("
  
  //make sure that no soln values used to simulate from the measurement model are negative
  
  double C_1_tmp, C_2_tmp, C_3_tmp, C_4_tmp, C_5_tmp;
  
  if (C_1 < 0) {
    C_1_tmp = 0;
  } else {
    C_1_tmp = C_1;
  }
  
  if (C_2 < 0) {
    C_2_tmp = 0;
  } else {
    C_2_tmp = C_2;
  }
  
  if (C_3 < 0) {
    C_3_tmp = 0;
  } else {
    C_3_tmp = C_3;
  }
  
  if (C_4 < 0) {
    C_4_tmp = 0;
  } else {
    C_4_tmp = C_4;
  }
  
  if (C_5 < 0) {
    C_5_tmp = 0;
  } else {
    C_5_tmp = C_5;
  }
  
  // define scalar multiples to take care of all types of reporting
  double s_1 = rho*rho_age_1*eta_a;
  double s_2 = rho*rho_age_2*eta_a;
  double s_3 = rho*rho_age_3*eta_a;
  double s_4 = rho*rho_age_4*eta_a;
  double s_5 = rho*rho_age_5*eta_a;
  double s_u = rho*rho_age_u*(1-eta_a);
  
  // define scalar multiples to take care of age specific/independent over-dispersion
  double p_sq_1 = (psi*psi)*(psi_1*psi_1);
  double p_sq_2 = (psi*psi)*(psi_2*psi_2);
  double p_sq_3 = (psi*psi)*(psi_3*psi_3);
  double p_sq_4 = (psi*psi)*(psi_4*psi_4);
  double p_sq_5 = (psi*psi)*(psi_5*psi_5);
  double p_sq_u = (psi*psi)*(psi_u*psi_u);
  
  // define the age-specific mean of the reporting distribution 
  double m_1 = s_1*C_1_tmp; 
  double m_2 = s_2*C_2_tmp; 
  double m_3 = s_3*C_3_tmp; 
  double m_4 = s_4*C_4_tmp; 
  double m_5 = s_5*C_5_tmp; 
  double m_u = s_u*(C_1_tmp + C_2_tmp + C_3_tmp + C_4_tmp + C_5_tmp);
  
  // define the the age-specific variance of the reporting distribution
  double v_1 = m_1*(1 - s_1 + p_sq_1*m_1);
  double v_2 = m_2*(1 - s_2 + p_sq_2*m_2);
  double v_3 = m_3*(1 - s_3 + p_sq_3*m_3);
  double v_4 = m_4*(1 - s_4 + p_sq_4*m_4);
  double v_5 = m_5*(1 - s_5 + p_sq_5*m_5);
  double v_u = m_u*(1 - s_u + p_sq_u*m_u);
  
  // defining dummy variables to take on likelihood values 
  double lik_D_1, lik_D_2, lik_D_3, lik_D_4, lik_D_5, lik_D_u;
  
  // some more book-keeping
  double tol = 1.0e-18;
  
  // assigning zero log-likelihood to missing data
  
  if(ISNA(D_1)) {
    lik_D_1 = (give_log) ? 0:1;
  } else {
      lik_D_1 = dnorm(D_1, m_1, sqrt(v_1)+tol, give_log);  
  }
  
  if(ISNA(D_2)) {
    lik_D_2 = (give_log) ? 0:1;
  } else {
      lik_D_2 = dnorm(D_2, m_2, sqrt(v_2)+tol, give_log);  
  }
  
  if(ISNA(D_3)) {
    lik_D_3 = (give_log) ? 0:1;
  } else {
      lik_D_3 = dnorm(D_3, m_3, sqrt(v_3)+tol, give_log);  
  }
  
  if(ISNA(D_4)) {
    lik_D_4 = (give_log) ? 0:1;
  } else {
      lik_D_4 = dnorm(D_4, m_4, sqrt(v_4)+tol, give_log);  
  }
  
  if(ISNA(D_5)) {
    lik_D_5 = (give_log) ? 0:1;
  } else {
      lik_D_5 = dnorm(D_5, m_5, sqrt(v_5)+tol, give_log);  
  }
  
  if(ISNA(D_u)) {
    lik_D_u = (give_log) ? 0:1;
  } else {
      lik_D_u = dnorm(D_u, m_u, sqrt(v_u)+tol, give_log);  
  }
  
  // calculating the final value of log likelihood
  lik = lik_D_1 + lik_D_2 + lik_D_3 + lik_D_4 + lik_D_5 + lik_D_u;     
  
")  



##### parameter names ######
state_names_est <- c(sprintf("S_%d",1:5), 
                     sprintf("E1_%d",1:5), sprintf("I1_%d",1:5),
                     sprintf("E2_%d",1:5), sprintf("I2_%d",1:5),
                     sprintf("V_%d",1:5), sprintf("C_%d", 1:5))


state_names_sim <- c(state_names_est,
                     sprintf("R_%d",1:5),
                     sprintf("A_%d", 1:5), sprintf("G_%d", 1:5))


param_names_est <- c("S_0", "I1_0", "I2_0",  
                     "q", sprintf("q_age_%d", 1:5), 
                     "gamma", "sigma", "beta1",
                     "rho",
                     sprintf("rho_age_%d", 1:5), "rho_age_u", 
                     "psi",
                     sprintf("psi_%d", 1:5), "psi_u", 
                     sprintf("Cv%d",1:25), "iota", 
                     sprintf("lcv%d", 1:5), 
                     "alpha", "dwan", "t_intro", "p_intro",
                     "epsilon1", "epsilon2",  
                     "agegroups")


param_names_sim <- c("R_0", 
                     param_names_est)




# testing parameter values 
##### setting up parameter values ###### setting up default q such that R0 is 10
param_vals_est <- c(S_0 = 1e-1, I1_0 = 1e-4, I2_0 = 0,
                    q = 1, 
                    q_age_1 = 1, q_age_2 = 1, q_age_3 = 1, q_age_4 = 1, q_age_5 = 1, 
                    gamma = 365.25/5, sigma = 365.25/17, 
                    dwan = Inf, 
                    epsilon1 = 0, epsilon2 = 0, 
                    rho = 1,
                    rho_age_1 = 1, rho_age_2 = 1, rho_age_3 = 1, rho_age_4 = 1, rho_age_5 = 1, rho_age_u = 1, 
                    psi = 1, 
                    psi_1 = 1, psi_2 = 1, psi_3 = 1, psi_4 = 1, psi_5 = 1, psi_u = 1,
                    agegroups = 5, Cv = C, alpha = 0.054, beta1 = 0.11,
                    lcv = c(1/5, 1/10, 1/10, 1/15, 1/40), iota = 1, 
                    t_intro = 3000, p_intro = 6)

param_vals_sim <- c(R_0 = 0, param_vals_est) 


##### build the pomp object #####
# Here the model is defined at an annual scale with a time step of that of a day (dt = 1/365.25)
# When extrapolating the simulation, default scale of simulation has been set to annual (temp_scale = 1)
# This can of course be changed as required by passing arguments

make_pomp <- function(..., 
                      start_t = 1800, 
                      extrapolate_simulation = FALSE, dt = 1/365.25,  
                      extra_start_t = 1900, extra_end_t = 2100, temp_scale = 1, 
                      covar, 
                      incidence_data = mumps_case_reports) {
  
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
           skeleton = vectorfield(vseir_skel_est),
           rinit = vseir_init_est,
           rmeasure = vseir_rmeas,
           dmeasure = vseir_dmeas,
           accumvars = c(sprintf("C_%d", 1:5)),
           covar = covariate_table(covar, times = "Year", order = "constant"),
           covarnames = c(sprintf("N_%d", 1:5), sprintf("MU_%d", 1:5), 
                          sprintf("IN_%d", 1:5), sprintf("OUT_%d", 1:5), 
                          sprintf("p%d", 1:2), "Births", "eta_a"),
           statenames = state_names_est,
           paramnames = param_names_est, ...)
      ) 
  
  } else {
    
    message("Generating simulations at a scale and range identical to the input data")
    po <- (
      mumps_inc_data %.>%
      pomp(., 
           times = "year", t0 = start_t, 
           skeleton = vectorfield(vseir_skel_est),
           rinit = vseir_init_est,
           rmeasure = vseir_rmeas,
           dmeasure = vseir_dmeas,
           accumvars = sprintf("C_%d", 1:5),
           covar = covariate_table(covar, times = "Year", order = "constant"),
           covarnames = c(sprintf("N_%d", 1:5), sprintf("MU_%d", 1:5), 
                          sprintf("IN_%d", 1:5), sprintf("OUT_%d", 1:5), 
                          sprintf("p%d", 1:2), "Births", "eta_a"),
           statenames = state_names_est,
           paramnames = param_names_est, ...)
      ) 
    
  }
  
  }


