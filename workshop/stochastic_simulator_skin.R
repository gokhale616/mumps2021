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
