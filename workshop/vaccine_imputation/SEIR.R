#loading packages####
source("../../00/src.R", chdir = TRUE)
# Defing the pomp object -------------------------------------------------------------------------------------

# process model: stochastic implementation 
seir_stoc_step <- "
  
  N = S + E + I + R; 
  
  /*===== Prelimnary definitions ======*/
  double r[8];
  double dN[8];
  
  
  /*====== Auxilliary variables =======*/
  double beta0 = (R0*gamma*(sigma+mu))/sigma;
  
  double beta = beta0*(1-beta1*sin(2*M_PI*t));

  double lambda = (beta/N)*(I + eta); 
  
  
  Reff = (S/N)*R0*(1-beta1*sin(2*M_PI*t));
  
  
  /*======== Rate matrix =============*/ 
  /*Birth flux rates for the system*/
  r[0] = nu*N;
  
  /*flux rates for the susceptible compartment*/
  r[1] = lambda; r[2] = mu;
  
  /*flux rates for the exposed compartment*/
  r[3] = sigma; r[4] = mu;
  
  /*flux rates for the infectious compartment*/
  r[5] = gamma; r[6] = mu; 

  /*flux rates for the recovered compartment*/
  r[7] = mu;
  
  /*======= Drawing Births and importation from a Poisson distribution =====*/ 
  dN[0] = rpois(r[0]*dt);
  
  /*====== Drawing other transistions from multinomial processes =====*/
  /*Draws: Susceptible*/
  reulermultinom(2, S, &r[1], dt, &dN[1]);
  
  /*Draws: Exposed*/
  reulermultinom(2, E, &r[3], dt, &dN[3]);
  
  /*Draws: Infectious*/
  reulermultinom(2, I, &r[5], dt, &dN[5]);

  /*Draws: Recovered*/
  reulermultinom(1, R, &r[7], dt, &dN[7]);

  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  S += dN[0] - dN[1] - dN[2];
  
  /*Compartmental shifts: Exposed*/
  E += dN[1] - dN[3] - dN[4];
  
  /*Compartmental shifts: Infectious*/
  I += dN[3] - dN[5] - dN[6];
  
  /*Compartmental shifts: Recovered*/
  R += dN[5] - dN[7];
  
  /*True incidence: dummy variable*/
  C += dN[5];
  
"

# Measurement model: Poisson 
# Drawing from the measurement process
seir_rmeasure <- "
  
  Cases = rpois(rho*C);
"
# Evaluating the density from the measurement process
seir_dmeasure <- "
  lik = dpois(Cases, rho*C, give_log);
"

seir_rinit <- "
  
  //initializing the system at the endemic equilibrium
  
  double beta0 = (R0*gamma*(sigma+mu))/sigma;
  
  double See = pop0/R0;
  double Eee = pop0*((mu*(mu+gamma))/(beta0*sigma))*(R0-1);
  double Iee = pop0*(mu/sigma)*(R0-1);
  double Ree = pop0 - See - Eee - Iee;
  
  S = nearbyint(See);
  E = nearbyint(Eee);
  I = nearbyint(Iee);
  R = nearbyint(Ree);
  
  C = 0;
  
  N = S + E + I + R;
  
  /* this is not right! */
  /* this requires a proper working function */
  Reff = 0;
" 
  
# State names:
state_names <- c("S", "E", "I", "R", "C", "N", "Reff")

# Parameter names:
rp_names <- c("sigma", "gamma", "R0", "eta", "nu", "mu", 
              "beta1", "rho", "pop0", "sigma_w")

# Variables to set to zero at every integration step or every data step?
zero_names <- c("C") 


# Using arbitrary values for the parameters (scale of integrator is expressed in years)
rp_vals <- c(gamma = 365.25/17, sigma = 365.25/5, R0 = 10, eta = 1, 
             nu = 1/80, mu = 1/80, beta1 = 0.11, rho = 0.04, sigma_w = 1e-5, 
             pop0 = 330e6)


# Constructing the pomp simulator
seir_po <- (
  data.frame(year = seq(0, 500, by = 1/52), 
             Cases = NA) %>% 
  pomp(t0 = -60, 
       times = "year",
       rprocess = euler(Csnippet(seir_stoc_step), delta.t = (1/100)), 
       rinit = Csnippet(seir_rinit), 
       dmeasure = Csnippet(seir_dmeasure), 
       rmeasure = Csnippet(seir_rmeasure), 
       accumvars = zero_names, 
       statenames = state_names, 
       paramnames = rp_names,
       params = rp_vals) 
  )

# View the pomp object
#spy(seir_po)       

# Simulating the values 
seir_sim <- seir_po %>% 
  simulate(., nsim = 1, format="d", include.data=FALSE) 


seir_sim %.>% 
  as_tibble(.) %>% 
  filter(., year >1) %.>% 
  select(., year, S, E, I, R, N, Cases, Reff) %.>% 
  gather(., key = "Compartment", value = "Comp_count", -c(year), factor_key = TRUE) %.>% 
  group_by(., year, Compartment) %.>% 
  calculate_quantile(., var = Comp_count) %.>% 
  ggplot(., 
         aes(x = year)) +
  geom_line(aes(y = `0.5`), color = "grey30") +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), fill = "grey30", alpha = 0.6) +
  facet_wrap(.~Compartment, scales = "free") +
  # scale_x_continuous(breaks = seq(2, 15, by = 3)) +
  project_theme +
  cap_axes
  
  













