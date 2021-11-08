#loading packages####
source("../../00/src.R", chdir = TRUE)
# Defing the pomp object -------------------------------------------------------------------------------------

# process model: stochastic implementation 
vsei_stoc_step <- "
  // redefine time by using the user defined variable
  double usr_t = (t+obs_t); 
  
  // make normal random draws for the Brownian bridge
  double dw = rnorm(0, dt*sigma_wn_scale);
  
  // define a temporal boundary to make sure that bridge state doen not update before t = 0. 
  // this will avoid eccentric evaluations
  double r;
  
  // update the brownian bridge here
  // add the curvature 
  if (usr_t < 0 || (bb_end_t - usr_t) < hack) {
    B += 0; 
     r = 0;
  } else {
      B += ((-B)/(bb_end_t - usr_t) + dw); 
      r  = (pow((usr_t/bb_end_t), k) + B);
  }
  

  
  // add a couple of conditions to make sure that the treated bb stays within the interval [0,b]
  // maybe a logit max function can be used here??!
  
  if((bb_end_t - usr_t) < hack || r > b) {
    p = b;
  } else if (r < 0) {
      p =  0;
  } else {
      p = r;  
  }
  
  
  /*===== Prelimnary definitions ======*/
  double r_births[2];
  double dN_births[2];
  
  double r_fromV[1]; 
  double dN_fromV[1]; 
  
  double r_fromS[2]; 
  double dN_fromS[2]; 
  
  double r_fromE[2]; 
  double dN_fromE[2]; 
  
  double r_fromI[2]; 
  double dN_fromI[2]; 
  
  
  /*====== Auxilliary variables =======*/
  double beta0 = (R0*gamma*(sigma+mu))/sigma;
  
  double beta = beta0*(1-beta1*sin(2*M_PI*usr_t));

  double lambda = (beta/pop)*(I + eta); 
  
  
  Reff = (S/pop)*R0*(1-beta1*sin(2*M_PI*usr_t));
  
  /*======== Rate matrix =============*/ 
  /* Birth flux rates - vaccinated and unvaccinated births */
  r_births[0] = nu*pop*((1-alpha)*p); r_births[1] = nu*pop*(1-(1-alpha)*p); 
  
  /* from V compartment */
  r_fromV[0] = mu;
  
  /* from S compartment */
  r_fromS[0] = lambda; r_fromS[1] = mu; 
  
  /* from E compartment */
  r_fromE[0] = sigma; r_fromE[1] = mu;
  
  /* from I compartment */
  r_fromI[0] = gamma; r_fromI[1] = mu;
  
  
  /*======= Drawing Births and importation from a Poisson distribution =====*/ 
  dN_births[0] = rpois(r_births[0]*dt);
  dN_births[1] = rpois(r_births[1]*dt);
  
  /*====== Drawing other transistions from multinomial processes =====*/
  /*Draws: Vaccinated*/
  reulermultinom(1, V, &r_fromV[0], dt, &dN_fromV[0]);
  
  /*Draws: Susceptible*/
  reulermultinom(2, S, &r_fromS[0], dt, &dN_fromS[0]);
  
  /*Draws: Exposed*/
  reulermultinom(2, E, &r_fromE[0], dt, &dN_fromE[0]);
  
  /*Draws: Infectious*/
  reulermultinom(2, I, &r_fromI[0], dt, &dN_fromI[0]);

  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  V += dN_births[0] - dN_fromV[0];

  /*Compartmental shifts: Susceptible*/
  S += dN_births[1] - dN_fromS[0] - dN_fromS[1];
  
  /*Compartmental shifts: Exposed*/
  E += dN_fromS[0] - dN_fromE[0] - dN_fromE[1];
  
  /*Compartmental shifts: Infectious*/
  I += dN_fromE[0] - dN_fromI[0] - dN_fromI[1];
  
  /*True incidence: dummy variable*/
  C += dN_fromI[0];
  
  //Rprintf(\"usr_t = %lg, S = %lg, V = %lg, mu_V = %lg, p = %lg, v_births = %lg, uv_births = %lg\\n\", 
  //          usr_t, S, V, dN_fromV[0], p, dN_births[0], dN_births[1]);
  
"

# Measurement model: Poisson 
# Drawing from the measurement process
vsei_rmeasure <- "
  
  double m = rho*C; 
  
  double v = m*(1 - rho + m*pow(psi,2)); 
  
  // a bit more book-keeeping 
  double tol = 1.0e-18; 
  
  cases = rnorm(m, sqrt(m)+tol);
"

# Evaluating the density from the measurement process
vsei_dmeasure <- "
  double m = rho*C; 
  double v = m*(1 - rho + m*pow(psi,2)); 
  
  // a bit more book-keeeping 
  double tol = 1.0e-18; 
  
  if(ISNA(cases)) {
    lik = (give_log) ? 0:1;
  } else {
      lik = dnorm(cases, m, sqrt(v)+tol, give_log);  
  }
"

vsei_rinit <- "
  
  // latent states
  V = V_0; 
  S = S_0;
  E = E_0;
  I = I_0;
  
  // observable state
  C = 0;
  
  
  // brownian bridge variable
  B = B_0;
  
  // vaccine coverage
  p = p_0; 
  
  // effective reproductive number
  Reff =  (S/pop)*R0;
" 


# State names:
state_names <- c("V", "S", "E", "I", "C", "B", "p", "Reff")

# Parameter names:
rp_names <- c("sigma", "gamma", "R0", "eta", "nu", "mu", "alpha",
              "beta1", "pop", "rho", "psi", "k", "sigma_wn_scale", "bb_end_t", "b", "hack", 
              "V_0", "S_0", "E_0", "I_0", "B_0", "p_0", "obs_t")

# Variables to set to zero at every integration step or every data step?
zero_names <- c("C") 


pop_val <- 219e6

# Using arbitrary values for the parameters (scale of integrator is expressed in years)
rp_vals <- c(gamma = 365.25/5, sigma = 365.25/13, R0 = 10, eta = 1, 
             alpha = 0.0, 
             nu = 1/80, mu = 1/80, beta1 = 0.15, pop = pop_val,  
             rho = 0.04, psi = 0.8,
             k = 0.4, sigma_wn_scale = 365/30, 
             bb_end_t = 16.99522, b = 0.86, hack = 0.50, 
             V_0 = 0,
             S_0 = round(1/10*pop_val), 
             E_0 = round(0.0004003011*pop_val), 
             I_0 = round(0.0001539356*pop_val), 
             B_0 = 0, p_0 = 0, 
             obs_t = 0)



make_pomp_vsei <- function(data) {
  
  data %.>% 
    pomp(.,
         t0 = 0, 
         times = "year",
         rprocess = euler(Csnippet(vsei_stoc_step), delta.t = (1/365.25)), 
         rinit = Csnippet(vsei_rinit), 
         dmeasure = Csnippet(vsei_dmeasure), 
         rmeasure = Csnippet(vsei_rmeasure), 
         accumvars = zero_names, 
         statenames = state_names, 
         paramnames = rp_names,
         params = rp_vals) 
}











