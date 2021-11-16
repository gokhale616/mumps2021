# load the pomp pre-utils and additional package for pre-processing
source("./process_tycho_data.R")
source("./VSEIR.R")

if(FALSE) {
# quick plot to decide cut-offs
mumps_weekly_case_reports %.>% 
  ggplot(., aes(x = year, y = cases)) +
  geom_line() +
  project_theme +
  scale_x_continuous(breaks = seq(0, 20, by = 4.25)) +
  cap_axes
}

# define a pomp object suing incidence data 
po_vseir <- (
  mumps_weekly_case_reports %.>% 
  make_pomp_vseir(.)
  )

# use the parameter vector to to play around the dynamics

# moving the end of the brownian bridge to suite the time data
test_sim <- po_vseir %.>% 
  simulate(., 
           nsim = 500, 
           include.data = TRUE, format = "d", seed = 986747881L) %.>% 
  mutate(., cases = ifelse(year == 0, NA, cases))
  

test_cases <- test_sim %.>%     
  select(., year, `.id`, cases) %.>% 
  mutate(., 
         ds = ifelse(`.id` == "data", "data", "simulation"))




# treat case data for compartment plots 
mumps_weekly_case_reports_treated <- (
  mumps_weekly_case_reports %.>% 
  mutate(., 
         Compartment = "cases")
  )


# if(FALSE){
# other compartments
sim_data <- (
  test_sim %.>% 
  filter(., 
         `.id` != "data") %.>% 
  select(., year, V, S, E, I, R, N, cases, Reff, p) %.>% 
  gather(.,
         key = "Compartment", value = "Comp_count", -c(year), 
         factor_key = TRUE) %.>% 
  group_by(., year, Compartment) %.>% 
  calculate_quantile(., var = Comp_count)
) 




if(FALSE){

# summarizing over the 
annual_p <- (
  sim_data %.>%  
    filter(., 
           Compartment == "p") %.>%
    select(., year, `0.5`) %.>%
    mutate(., 
           Date = mumps_weekly_case_reports_treated$PeriodMidDate, 
           year_discrete = ifelse(is.na(year(Date)) == TRUE, 1967, year(Date)), 
           p = `0.5`) %.>%
    select(.,
           year, year_discrete, p) %.>% 
    group_by(., year_discrete) %.>% 
    mutate(., 
           coverage = median(p), 
           Compartment = "p") %.>% 
    ungroup(.) 
)

sim_data %.>%   
  mutate_at(., 
            .vars = c("0.025", "0.5", "0.975"), 
            .funs = function(x){
              ifelse(.$Compartment %in% c("V", "S", "E", "I", "R", "N"), x/rp_vals["pop0"], x)
            }
  ) %.>% 
  ggplot(., aes(x = year)) +
  geom_line(aes(y = `0.5`, color = "#A83279")) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`, fill = "#A83279"), alpha = 0.6) +
  geom_line(data = mumps_weekly_case_reports_treated, 
            aes(x = year, y = cases, colour = "grey30")) +
  geom_ribbon(data = mumps_weekly_case_reports_treated,
              aes(x = year,
                  ymin = `cases`, ymax = `cases`, fill = "grey30"), alpha = 0.0) +
  geom_point(data = annual_p, aes(x = year, y = coverage, shape = "p"), 
             colour = "grey30", size = 1) +
  labs(y = "Latent state value", 
       x = "Time (Year)") +
  facet_wrap(.~Compartment, scales = "free") +
  scale_x_continuous(breaks = seq(0, 18, by = 6), 
                     limits = c(0, 18)) +
  scale_shape_manual("", values = 21, label = "Annual Coverage") +
  scale_colour_manual("", 
                      values = c("#A83279", "grey30"), 
                      labels = c("Simulation", "Data")) +
  scale_fill_manual("", 
                    values = c("#A83279", "grey30"), 
                    labels = c("Simulation", "Data")) +
  project_theme +
  cap_axes +
  guides(fill = guide_legend(override.aes = list(fill = c("#A83279", NA))))


# extract initial conditions for the one step simulator 

# init_for_vsei <- (
#   sim_data %.>% 
#   filter(., Compartment %in% c("V", "S", "E", "I")) %.>% 
#   select(., year, Compartment, `0.5`) %.>% 
#   mutate(., `0.5` = round(`0.5`)) %.>% 
#   spread(., key  = Compartment, value = `0.5`) %.>% 
#   filter(., year == 0)
#   )
# 
# save(init_for_vsei, file = "./init_for_vsei.rds")


}


  
  
  


