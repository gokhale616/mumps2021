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

# moving the end of the brownian bridge to suite the time data
test_sim <- po_vseir %.>% 
  simulate(., 
           nsim = 50, 
           include.data = TRUE, format = "d", seed = 986747881L) %.>% 
  mutate(., cases = ifelse(year == 0, NA, cases))
  

test_cases <- test_sim %.>%     
  select(., year, `.id`, cases) %.>% 
  mutate(., 
         ds = ifelse(`.id` == "data", "data", "simulation"))


test_cases %.>%  
  ggplot(., aes(x = year, y = cases, group = `.id`, colour = ds)) +
  geom_line() +
  labs(y = "Cases", 
       x = "Time (Year)", 
       colour = "") +
  scale_colour_manual(values = c("grey30", "#A83279")) +
  project_theme + 
  cap_axes



# other compartments
test_sim %.>% 
  filter(., 
         `.id` != "data") %.>% 
  select(., year, V, S, E, I, R, N, cases, Reff, p) %.>% 
  gather(., key = "Compartment", value = "Comp_count", -c(year), factor_key = TRUE) %.>% 
  group_by(., year, Compartment) %.>% 
  calculate_quantile(., var = Comp_count) %.>% 
  ggplot(., 
         aes(x = year)) +
  geom_line(aes(y = `0.5`), color = "#A83279") +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), fill = "#A83279", alpha = 0.6) +
  labs(y = "Latent state value", 
       x = "Time (Year)") +
  facet_wrap(.~Compartment, scales = "free") +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  project_theme +
  cap_axes




