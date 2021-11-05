#loading packages####
source("../../00/src.R", chdir = TRUE)
# Defing the pomp object -------------------------------------------------------------------------------------

# process model: stochastic implementation 
vseir_stoc_step <- "
  
  // make normal random draws for the Brownian bridge
  double dw = rnorm(0, dt*sigma_wn_scale);
  
  // define a temporal boundary to make sure that bridge state doen not update before t = 0. 
  // this will avoid eccentric evaluations
  double r;

  // update the brownian bridge here
  // add the curvature 
  if (t < 0) {
    B += 0; 
     r = 0;
  } else {
      B += ((-B)/(bb_end_t-t) + dw); 
      r  = (pow((t/bb_end_t), k) + B);
  }
  

  
  // add a couple of conditions to make sure that the treated bb stays within the interval [0,b]
  // maybe a logit max function can be used here??!
  
  if((bb_end_t-t) < hack || r > b) {
    p = b;
  } else if (r < 0) {
      p =  0;
  } else {
      p = r;  
  }
  
  
  // total population size
  N = S + E + I + R + V; 
  
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
  
  double r_fromR[1]; 
  double dN_fromR[1]; 
  
  
  /*====== Auxilliary variables =======*/
  double beta0 = (R0*gamma*(sigma+mu))/sigma;
  
  double beta = beta0*(1-beta1*sin(2*M_PI*(t-phi)/1.));

  double lambda = (beta/N)*(I + eta); 
  
  
  Reff = (S/N)*R0*(1-beta1*sin(2*M_PI*t));
  
  /*======== Rate matrix =============*/ 
  /* Birth flux rates - vaccinated and unvaccinated births */
  r_births[0] = nu*N*((1-alpha)*p); r_births[1] = nu*N*(1-(1-alpha)*p); 
  
  /* from V compartment */
  r_fromV[0] = mu;
  
  /* from S compartment */
  r_fromS[0] = lambda; r_fromS[1] = mu; 
  
  /* from E compartment */
  r_fromE[0] = sigma; r_fromE[1] = mu;
  
  /* from I compartment */
  r_fromI[0] = gamma; r_fromI[1] = mu;
  
  /* from R compartment */
  r_fromR[0] = mu;
  
  
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

  /*Draws: Recovered*/
  reulermultinom(1, R, &r_fromR[0], dt, &dN_fromR[0]);

  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  V += dN_births[0] - dN_fromV[0];

  /*Compartmental shifts: Susceptible*/
  S += dN_births[1] - dN_fromS[0] - dN_fromS[1];
  
  /*Compartmental shifts: Exposed*/
  E += dN_fromS[0] - dN_fromE[0] - dN_fromE[1];
  
  /*Compartmental shifts: Infectious*/
  I += dN_fromE[0] - dN_fromI[0] - dN_fromI[1];
  
  /*Compartmental shifts: Recovered*/
  R += dN_fromI[0] - dN_fromR[0];
  
  /*True incidence: dummy variable*/
  C += dN_fromI[0];
  
"

# Measurement model: Poisson 
# Drawing from the measurement process
vseir_rmeasure <- "
  
  double m = rho*C; 
  
  double v = m*(1 - rho + m*pow(psi,2)); 
  
  // a bit more book-keeeping 
  double tol = 1.0e-18; 
  
  cases = rnorm(m, sqrt(m)+tol);
"

# Evaluating the density from the measurement process
vseir_dmeasure <- "
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

vseir_rinit <- "
  
  //initializing the system at the endemic equilibrium
  double beta0 = (R0*gamma*(sigma+mu))/sigma;
  
  double See = pop0/R0;
  double Eee = pop0*((mu*(mu+gamma))/(beta0*sigma))*(R0-1);
  double Iee = pop0*(mu/sigma)*(R0-1);
  double Ree = pop0 - See - Eee - Iee;
  
  // latent states
  V = 0; 
  S = nearbyint(See);
  E = nearbyint(Eee);
  I = nearbyint(Iee);
  R = nearbyint(Ree);
  
  // observable state
  C = 0;
  
  // population size
  N = S + E + I + R + V;
  
  // brownian bridge variable
  B = 0;
  
  // vaccine coverage
  p = 0; 
  
  // effective reproductive number
  Reff =  (S/N)*R0;
" 


# State names:
state_names <- c("V", "S", "E", "I", "R", "C", "N", "B", "p", "Reff")

# Parameter names:
rp_names <- c("sigma", "gamma", "R0", "eta", "nu", "mu", "alpha",
              "beta1", "phi", "pop0", "rho", "psi", "k", "sigma_wn_scale", "bb_end_t", "b", "hack")

# Variables to set to zero at every integration step or every data step?
zero_names <- c("C") 


# Using arbitrary values for the parameters (scale of integrator is expressed in years)
# defaults are visually a good fit to the data! 
rp_vals <- c(gamma = 365.25/5, sigma = 365.25/13, R0 = 10, eta = 1, 
             alpha = 0.0, 
             nu = 1/80, mu = 1/80, beta1 = 0.13, pop0 = 219e6,  
             rho = 0.06, psi = 0.8, phi = 0.25, 
             k = 0.5, sigma_wn_scale = 365/25, 
             bb_end_t = 17, b = 0.86, hack = 0.50)

# c("beta1", "phi", "R0", "sigma", "k", "rho", "psi", "sigma_wn_scale")
# c(0.13, 0.25, 10, 365.25/13, 0.5, 0.06, 0.8, 365.25/25)

make_pomp_vseir <- function(data) {
  
  data %.>% 
    pomp(.,
         t0 = -60, 
         times = "year",
         rprocess = euler(Csnippet(vseir_stoc_step), delta.t = (1/365.25)), 
         rinit = Csnippet(vseir_rinit), 
         dmeasure = Csnippet(vseir_dmeasure), 
         rmeasure = Csnippet(vseir_rmeasure), 
         accumvars = zero_names, 
         statenames = state_names, 
         paramnames = rp_names,
         params = rp_vals) 
}






