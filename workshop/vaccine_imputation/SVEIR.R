# The goal of this script is to emulate incorporation of demographic and vaccine covariates in a SVEIR model. 
# This is to make easy working with similar co-variates when age-structure is involved. 
# Objectives include succefully including covariates and finding a way to incoroporate initial conditions 
 
##############################################################################################################

#loading packages####
sapply(c("tidyverse", "pomp", "reshape2", "plyr", 
         "ggthemes","grid", "gridExtra", "magrittr"), 
       require, character.only = TRUE)
theme_set(theme_bw())

set.seed(594709947L)

source("util_fns.R")

# loading a dataframe of covariates --------------------------------------------------------------------------

covar <- read_csv(file = "../processed_data/covar_5_1910_2015.csv")
load("../processed_data/mumps_cr_ma.Rdata")

# setting up covariates the to copy the stpes in the bigger model
# break open mu parameter: note the currently all values of mu are negative so it makes no sense to carry out  
# this step. It is only carried out to replicate steps while dealing with a more complex demogrpahic model.
mumps_1977_2018_ma %>%
  transmute(Year = Year, Cases = total) %>% 
  bind_rows(tibble(Year = 1976, Cases = NA)) %>% 
  arrange(Year) -> total_mumps_reports

covar %>% 
  transmute(p = p1, 
            N = combined, 
            delta_N = diff_combined, 
            Births = nu*N,
            mu = (delta_N - Births)/(N), 
            Year = year, 
            nu = nu) -> covar_primary 
  
covar_primary %>%   
  rowwise() %>% 
  do(break_mu_into_2(.$mu)) %>% 
  as_tibble() -> covar_migr
  
covar_migr %>%   
  bind_cols(covar_primary) %>% 
  select(8:9, 1:4, 6) %>% 
  mutate(efflux = -efflux) -> covar_primary

covar_primary %>% 
  gather(key = "Covariate", value = "rates", -Year) %>% 
  mutate(Covariate = case_when(
    Covariate == "efflux" ~ "mu[efflux]", 
    Covariate == "influx" ~ "mu[influx]",
    Covariate == "p" ~ "Coverage",
    Covariate == "N" ~ "Population", 
    TRUE ~ Covariate)) %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = rates), size = 0.8) +
  facet_wrap(.~Covariate, scales = "free", labeller = label_parsed, nrow = 2) +
  labs(y = "") +
  theme(aspect.ratio = 0.5, 
        strip.background = element_blank()) -> covariates_plts

total_mumps_reports %>% 
  ggplot(aes(x = Year)) +
  geom_line(aes(y = Cases), size = 0.8, colour = "#dd1818") +
  theme(aspect.ratio = 0.5) -> case_reports_plts


ggpubr::ggarrange(covariates_plts, case_reports_plts, nrow = 2) -> combined_data_plt

ggsave(combined_data_plt, file = "combined_data_plt.pdf", 
       width = 11.69, height = 9)

# Defing the pomp object -------------------------------------------------------------------------------------
# process model: stochastic implementation 
sveir_stoc_step <- "
  /*===== Prelimnary definitions ======*/
  /*double births;*/
  double r[15];
  double dN[15];
  
  /*====== Auxilliary variables =======*/
  double beta0 = R0/Dinf;
  double beta = beta0*(1-beta1*sin(2*M_PI*t));
  
  double sigma = 1/Dlat;
  double gamma = 1/Dinf;
  double lambda = (beta/N)*(I + eta); 
  
  /*======= Drawing Births from a Poisson distribution =====*/ 
  double uv_births = rpois(Births*(1-(1-alpha)*p)*dt);
  double  v_births = rpois(Births*(1-alpha)*p*dt);  
  
  /*======== Rate matrix =============*/ 
  /*flux rates for the susceptible compartment*/
  r[0] = lambda; r[1] = influx; r[2] = efflux;
  
  /*flux rates for the exposed compartment*/
  r[3] = sigma; r[4] = influx; r[5] = efflux;
  
  /*flux rates for the infectious compartment*/
  r[6] = gamma; r[7] = influx; r[8] = efflux; 
  
  /*flux rates for the recovered compartmemt*/
  r[9] = influx; r[10] = efflux; 
  
  /*flux rates for the vaccinated compartment*/
  r[11] = epsilon*lambda; r[12]  = delta, r[13] = influx; r[14] = efflux;
  
  /*====== Drawing other transistions from multinomial processes =====*/
  
  /*Draws: Susceptible*/
  reulermultinom(3, S, &r[0], dt, &dN[0]);
  
  /*Draws: Exposed*/
  reulermultinom(3, E, &r[3], dt, &dN[3]);
  
  /*Draws: Infectious*/
  reulermultinom(3, I, &r[6], dt, &dN[6]);
 
  /*Draws: Infectious*/
  reulermultinom(2, R, &r[9], dt, &dN[9]);
 
  /*Draws: Vaccinated*/
  reulermultinom(4, V, &r[11], dt, &dN[11]);
  
  /*====== System of Differential equations =====*/
  /*Compartmental shifts: Susceptible*/
  S += uv_births + dN[12] - (dN[0] - dN[1] + dN[2]);
  
  /*Compartmental shifts: Exposed*/
  E += dN[0] + dN[11] - (dN[3] - dN[4] + dN[5]);
  
  /*Compartmental shifts: Infectious*/
  I += dN[3] - (dN[6] - dN[7] + dN[8]);
  
  /*Compartment shifts: Recovered*/
  R += dN[6] + dN[9] - dN[10];  
  
  /*Compartmental shifts: Vaccinated*/
  V += v_births - (dN[11] + dN[12] - dN[13] + dN[14]);
  
  
  /*True incidence: dummy variable*/
  C += dN[6];
"

# Measurement model: Poisson 
# Drawing from the measurement process
sveir_rmeasure <- "Cases = rpois(rho*C);"
# Evaluating the density from the measurement process
sveir_dmeasure <- "lik = dpois(Cases, rho*C, give_log);"

sveir_rinit <- "
  S = N*S_0;
  E = 0;
  I = N*I_0;
  R = N*R_0;
  V = 0;
  
  C = 0;
" 
  
# State names:
state_names <- c("S", "E", "I" , "R", "V", "C")#, "N", "Births")

# Parameter names:
rp_names <- c("Dlat", "Dinf", "delta", "R0", "eta","beta1", "rho", "epsilon", "alpha")

ip_names <- c("S_0", "I_0", "R_0")

# Variables to set to zero at every integration step or every data step?
zero_names <- c("C") 


# Using arbitrary values for the parameters (scale of integrator is expressed in years)
rp_vals <- c(Dlat = 17/365.25, Dinf = 5/365.25, delta = 0/50, R0 = 2, eta = 0, 
             beta1 = 0.11, rho = 0.02, epsilon = 0, alpha = 0.0)

# Using values very close to the endemic equilibrium
ip_vals <- c(S_0 = 1e-1, I_0 = 1e-3, R_0 = 8e-1)



# Constructing the pomp simulator
mumps_1977_2018_ma %>%
  transmute(Year = Year, Cases = total) %>% 
  bind_rows(tibble(Year = 1850:1976, Cases = NA)) %>% 
  arrange(Year) %>% 
  pomp(t0 = 1600, 
       times = "Year",
       rprocess = euler(Csnippet(sveir_stoc_step), delta.t = (1/365.25*24)), 
       rinit = Csnippet(sveir_rinit), 
       dmeasure = Csnippet(sveir_dmeasure), 
       rmeasure = Csnippet(sveir_rmeasure),
       covar = covariate_table(covar_primary, times = "Year", 
                               order = "constant"),
       accumvars = zero_names, 
       statenames = state_names, 
       paramnames = c(rp_names, ip_names),
       params = c(rp_vals, ip_vals)) -> sveir_po


# View the pomp object
# spy(sveir_po)       

# Simulating the values 
sveir_po %>% 
  simulate(nsim = 1000, format="d", include.data=FALSE) -> sveir_sim

c_levels <- c("S", "E", "I", "R", 
              "V", "C", "Cases", "N", "R[N]", "N[M]")


# impute the population to check the class of sero-poitives

covar_primary %>% 
  select(Year, N) %>% 
  slice(rep(1, times = 63), 1:n()) %>% 
  mutate(Year = seq(1850, 2018)) -> annual_pop_used


sveir_sim %>% 
  select(-.id) %>% 
  as_tibble() %>% 
  group_by(Year) %>% 
  summarise_all(list(~mean(.))) %>% 
  right_join(annual_pop_used) %>% 
  mutate(`R[N]` = N - (S + E + I + V), 
         `N[M]` = (S + E + I + V + R)) %>% 
  filter(Year > 1960) %>% 
  gather(key = Compartment, value = Comp_count, -c(Year)) %>%
  mutate(Compartment = factor(Compartment, levels = c_levels), 
         Comp_colour = ifelse(Compartment == "Cases", TRUE, FALSE)) %>% 
  ggplot(aes(x = Year, y = Comp_count, fill = Comp_colour)) +
  labs(y = "Number") +
  geom_area(position = position_dodge(width = 0.8), alpha = 0.8) +
  scale_fill_manual(values = c("#333333", "#dd1818")) +
  #scale_colour_manual(values = c("#203A43", "red")) +
  theme_bw() +
  facet_wrap(.~Compartment, scales = "free", nrow = 4, 
             labeller = label_parsed) +
  theme(aspect.ratio = 0.5, 
        strip.background = element_blank(), 
        legend.position = "none")







