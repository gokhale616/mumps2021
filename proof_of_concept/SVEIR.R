#loading packages####
source("../00/src.R", chdir = TRUE)
set.seed(98674881L)
# Defing the pomp object -------------------------------------------------------------------------------------

# process model: stochastic implementation 
sveir_stoc_step <- "
  /*===== Prelimnary definitions ======*/
  double r[11];
  double dN[11];
  double N = S + E + I + R + V;
  
  /*====== Auxilliary variables =======*/
  double delta = 1/dwan;
  
  double p_val;
  if (t < t_vacc) {
    p_val = 0;
  } else {
      p_val = p;
  }
  
  double beta0 = R0*gamma;
  double beta = beta0*(1-beta1*sin(2*M_PI*t));

  double lambda = (beta/N)*(I + iota); 
  
  /*======== Rate matrix =============*/ 
  /*Birth flux rates for the system*/
  r[0] = nu*N*(1-(1-alpha)*p_val);
  r[1] = nu*N*(1-alpha)*p_val;
  
  /*flux rates for the susceptible compartment*/
  r[2] = lambda; r[3] = mu;
  
  /*flux rates for the exposed compartment*/
  r[4] = sigma; r[5] = mu;
  
  /*flux rates for the infectious compartment*/
  r[6] = gamma; r[7] = mu; 

  /*flux rates for the recovered compartment*/
  r[8] = mu;
  
  /*flux rates for the vaccinated compartment*/
  r[9] = delta; r[10] = mu;
  
  /*======= Drawing vaccinated and unvaccinated births from a Poisson distribution =====*/ 
  dN[0] = rpois(r[0]*dt);
  dN[1] = rpois(r[1]*dt);
  
  /*====== Drawing other transistions from multinomial processes =====*/
  /*Draws: Susceptible*/
  reulermultinom(2, S, &r[2], dt, &dN[2]);
  
  /*Draws: Exposed*/
  reulermultinom(2, E, &r[4], dt, &dN[4]);
  
  /*Draws: Infectious*/
  reulermultinom(2, I, &r[6], dt, &dN[6]);

  /*Draws: Recovered*/
  reulermultinom(1, R, &r[8], dt, &dN[8]);
  
  /*Draws: Vaccinated*/
  reulermultinom(2, V, &r[9], dt, &dN[9]);
  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  S += dN[0] + dN[9] - dN[2] - dN[3];
  
  /*Compartmental shifts: Exposed*/
  E += dN[2] - dN[4] - dN[5];
  
  /*Compartmental shifts: Infectious*/
  I += dN[4] - dN[6] - dN[7];
  
  /*Compartmental shifts: Recovered*/
  R += dN[6] - dN[8];
  
  /*Compartmental shifts: Recovered*/
  V += dN[1] - dN[9] - dN[10];
  
  /*True incidence: dummy variable*/
  C += dN[6];
  
  /*Vaccine uptake: dummy variable*/
  P = p_val;
"

# process model: stochastic implementation 
sveir_skel <- "
  /*===== Prelimnary definitions ======*/
  double N = S + E + I + R + V;
  
  /*====== Auxilliary variables =======*/
  double delta = 1/dwan;
  
  double p_val;
  if (t < t_vacc) {
    p_val = 0;
  } else {
      p_val = p;
  }
  
  double beta0 = R0*gamma;
  double beta = beta0*(1-beta1*sin(2*M_PI*t));

  double lambda = (beta/N)*(I + iota); 
  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  DS = nu*N*(1-(1-alpha)*p_val) + delta*V - (lambda + mu)*S;
  
  /*Compartmental shifts: Exposed*/
  DE = lambda*S - (sigma + mu)*E;
  
  /*Compartmental shifts: Infectious*/
  DI = sigma*E - (gamma + mu)*I;
  
  /*Compartmental shifts: Recovered*/
  DR = gamma*I - mu*R;
  
  /*Compartmental shifts: Recovered*/
  DV = nu*N*(1-alpha)*p_val - (delta - mu)*V;
  
  /*True incidence: dummy variable*/
  DC = gamma*I;
  
  /*Vaccine uptake: dummy variable*/
  DP = p_val;
"


sveir_rmeas <- "
  
  //make sure that no soln values used to simulate from the measurement model are negative
  
  double C_tmp;
  
  if (C < 0) {
    C_tmp = 0;
  } else {
    C_tmp = C;
  }
  
  // define scalar multiples to take care of over-dispersion
  double p_sq = psi*psi;
  
  // define the age-specific mean of the reporting distribution 
  double m = C_tmp*rho; 
  
  // define the the age-specific variance of the reporting distribution
  double v = m*(1 - rho + p_sq*m);
  
  // a bit more book-keeeping 
  double tol = 1.0e-18; 
      
  double Cases_tmp = rnorm(m, sqrt(v)+tol); 
      
  if (Cases_tmp < 0) {
    Cases = 0; 
  } else {
    Cases = nearbyint(Cases_tmp);
  }
    
"


sveir_dmeas <- "
  //make sure that no soln values used to simulate from the measurement model are negative
  
  double C_tmp;
  
  if (C < 0) {
    C_tmp = 0;
  } else {
    C_tmp = C;
  }
  
  
  // define scalar multiples to take care of over-dispersion
  double p_sq = psi*psi;
  
  // define the age-specific mean of the reporting distribution 
  double m = C_tmp*rho; 
  
  // define the the age-specific variance of the reporting distribution
  double v = m*(1 - rho + p_sq*m);
  
  
  // some more book-keeping
  double tol = 1.0e-18;
  
  // assigning zero log-likelihood to missing data
  
  if(ISNA(Cases)) {
    lik = (give_log) ? 0:1;
  } else {
      lik = dnorm(Cases, m, sqrt(v)+tol, give_log);  
  }
  
"


sveir_rinit <- "
  S = nearbyint(S_0);
  E = nearbyint(E_0);
  I = nearbyint(I_0);
  R = nearbyint(R_0);
  V = 0;
  C = 0;
  P = 0;
" 

# State names:
state_names <- c("S", "E", "I", "R", "V", "C", "P")

# Parameter names:
rp_names <- c("sigma", "gamma", "R0", "iota", "nu", "mu", 
              "beta1", "dwan", "alpha", "p", "t_vacc", "rho", "psi")

ip_names <- c("S_0", "E_0", "I_0", "R_0")

# Variables to set to zero at every integration step or every data step?
zero_names <- c("C", "P") 


# Using arbitrary values for the parameters (scale of integrator is expressed in years)
rp_vals <- c(gamma = 365.25/5, sigma = 365.25/25, R0 = 13, iota = 1, 
             nu = 1/80, mu = 1/80, beta1 = 0.11, 
             alpha = 0.54, dwan = 100, p = 0.91, t_vacc = 20, 
             rho = 1e-3, psi = 0.1)

# Using values very close to the endemic equilibrium
ip_vals <- c(S_0 = 1e5, E_0 = 0, I_0 = 1e2, R_0 = 899900)*330

# Constructing the pomp simulator
data.frame(Year = seq(0, 45, by = 1), 
           Cases = NA) %>% 
  pomp(t0 = -1000, 
       times = "Year",
       skeleton = vectorfield(Csnippet(sveir_skel)),
       rprocess = euler(Csnippet(sveir_stoc_step), delta.t = 1/1e4), 
       rinit = Csnippet(sveir_rinit), 
       dmeasure = Csnippet(sveir_dmeas), 
       rmeasure = Csnippet(sveir_rmeas), 
       accumvars = zero_names, 
       statenames = state_names, 
       paramnames = c(rp_names, ip_names),
       params = c(rp_vals, ip_vals)) -> sveir_po

data.frame(Year = seq(0, 45, by = 1), 
           Cases = NA) %>% 
  pomp(t0 = -1000, 
       times = "Year",
       skeleton = vectorfield(Csnippet(sveir_skel)),
       rinit = Csnippet(sveir_rinit), 
       dmeasure = Csnippet(sveir_dmeas), 
       rmeasure = Csnippet(sveir_rmeas), 
       accumvars = zero_names, 
       statenames = state_names, 
       paramnames = c(rp_names, ip_names),
       params = c(rp_vals, ip_vals)) -> sveir_det_po
