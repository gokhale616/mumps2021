# load the pomp pre-utils and additional package for pre-processing
source("./VSEIR.R")
library("lubridate")

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
           year = (year_val-1968)+week_val/52.25) %.>% 
    select(., PeriodMidDate, year, cases) %.>% 
    # add a missing row for pomp
    bind_rows(., 
              tibble(PeriodMidDate = NA, 
                     year = 0, 
                     cases = NA)
    ) %.>% 
    arrange(., year)
    
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
  select(., -PeriodMidDate) %.>% 
  make_pomp_vseir(.)
  )

# check if the data was properly loaded
# spy(po_vseir)

# generate a particle filter object for the given pomp object
pfilter_vseir <- (
  po_vseir %.>% 
  pfilter(., 
          Np = 1e3, 
          params = rp_vals, 
          dmeasure = Csnippet(vseir_dmeasure),
          paramnames = rp_names, 
          statenames = state_names)
  )


logLik(pfilter_vseir)

plot(pfilter_vseir)



















