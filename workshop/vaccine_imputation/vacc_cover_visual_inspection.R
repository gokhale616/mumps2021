# load the pomp pre-utils and additional package for pre-processing
source("./VSEIR.R")
library("lubridate")
library("future.apply")

# load and process data for vaccine imputation
mumps_weekly_case_reports <- (
  read_csv("../../raw_data/tycho_20201118-150400.csv") %.>%
    # select relevant variables for forming the time series
    select(.,
           PeriodStartDate, 
           PeriodEndDate,
           CountValue) %.>% 
    # find week midpoints for plotting convenience
    mutate(.,  
           PeriodMidDate = PeriodStartDate + (PeriodEndDate-PeriodStartDate)/2
           ) %.>% 
    select(., 
           -c(PeriodStartDate, PeriodEndDate)) %.>% 
    # summarize case report volume over all states
    group_by(., 
             PeriodMidDate) %.>% 
    summarise(., cases = sum(CountValue, na.rm = TRUE)) %.>% 
    ungroup(.) %.>% 
    arrange(., 
            PeriodMidDate) %.>% 
    # define a period to run the analysis
    filter(., 
           PeriodMidDate > as.Date("1967-12-01") & PeriodMidDate < as.Date("1985-01-01")) %.>% 
    # define a year variable for pomp
    mutate(., 
           year_val = year(PeriodMidDate), 
           week_val = week(PeriodMidDate), 
           year = (year_val-1968)+week_val/52
           ) %.>%  
    select(., PeriodMidDate, year, cases) %.>% 
    # add a missing row for pomp
    bind_rows(., 
              tibble(PeriodMidDate = NA, 
                     year = 0, 
                     cases = NA)
    )  %.>% 
    arrange(., year) %.>% 
    mutate(., 
           year = seq(0, 17, length.out = nrow(.)))
)



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
# }

  
}


  
  
  


