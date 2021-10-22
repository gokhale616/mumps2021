# load the setup
# depends on preprocess_demog_data_5ac.R" 

# look at tycho data
non_cumulative_mumps <- read_csv(file = "../raw_data/tycho_20201118-150400.csv")

mumps_by_states <- non_cumulative_mumps %.>% 
  transmute(.,
            State = Admin1Name, 
            Date = PeriodEndDate, 
            PeriodLength = PeriodEndDate-PeriodStartDate,
            Cases = CountValue) 

# producing the a sequence of the state-year to generate missing values 
n_states <- length(unique(mumps_by_states$State))
unique_states <- unique(mumps_by_states$State)

date_range <- c(min(mumps_by_states$Date), max(mumps_by_states$Date))

# producing the a sequence of the dates to generate missing values 
complete_time_state_weekly <- tibble(Date = rep(seq.Date(from = date_range[1], to = date_range[2], 
                                                       by = "week"), each = n_states), 
                                     State = rep(unique_states, 
                                                 times = (difftime(date_range[2], date_range[1], 
                                                                   units = "weeks")+1)
                                                 )
                                     ) %.>%
  mutate(., 
         State = as.character(State))


# Form a data set with missing values but complete week-state records  
complete_mumps_by_states_weekly <- complete_time_state_weekly %.>% 
  left_join(., 
            mumps_by_states, by = c("Date", "State")) %.>%
  mutate(., 
         Year = format(Date, "%Y"), 
         Month = format(Date, "%m"), 
         Day = format(Date, "%d"),
         State = str_to_title(State))


# function to distinguishing missing values from a zero sum
distinct_na_sum <- function(x) {
  
  x_sum <- sum(x, na.rm = TRUE)
  
  ifelse(x_sum == 0, NA, x_sum)
  
}


# summarize mumps cases by year-state 
complete_mumps_by_states_annual <- complete_mumps_by_states_weekly %.>% 
  group_by(., 
           Year, State) %.>% 
  summarise(., 
            Cases = distinct_na_sum(Cases)) %.>%
  ungroup(.) %.>%
  mutate(.,
         # Cases = ifelse(Cases == 0, NA, Cases),
         Year = as.numeric(Year)
         )

# summarize mumps by week (sum over all states)
complete_mumps_us_weekly <- complete_mumps_by_states_weekly %.>%
  transmute(., 
            Date, Cases) %.>%
  group_by(., Date) %.>%
  summarise(.,
            Cases = distinct_na_sum(Cases)) %.>% 
  ungroup(.) %.>%
  mutate(., 
         DateDiff = Date - lag(Date))


# summarize mumps incidence by year (sum overl all states)
complete_mumps_us_annual <- complete_mumps_us_weekly %.>%
  mutate(., 
         Year = format(Date, "%Y")) %.>% 
  group_by(., 
           Year) %.>% 
  summarise(., 
            Cases = distinct_na_sum(Cases)) %.>%
  ungroup(.) %.>%
  mutate(.,
         Year = as.numeric(Year)
  )

if(FALSE) {
# plots stratified by states 
complete_mumps_by_states_annual %.>%
  ggplot(., aes(x = Year, y = Cases)) +
  geom_line(size =  0.8) +
  facet_wrap(.~State) +
  project_theme +
  theme()
  
complete_mumps_by_states_weekly %.>%
  ggplot(., aes(x = Date, y = Cases)) +
  geom_line(size = 0.6) +
  scale_colour_manual(name = "", values = distinctColorPalette(56)) +
  facet_wrap(.~State) +
  project_theme 

    


}

# fusing case data with population sizes and geographical information
# summarize the demog data by states - get rid of the rogue (KR state for 2005) - 
# Find which state is missing a 2005

demog_by_states <- demog_data_1969_2017 %.>% 
  select(., 
         year, state, population) %.>% 
  group_by(., 
           year, state) %.>% 
  summarise(., 
            population = sum(population, na.rm = TRUE)) %.>% 
  ungroup(.) 

# pick the states' full names and abrevation for seasmless integration with the case data    
state_names_abbr <- statepop %.>%
  transmute(., 
            state = abbr, 
            full = full) 


length(state_names_abbr$state)
length(unique(demog_by_states$state))

unique(demog_by_states$state)[which(unique(demog_by_states$state) %nin% unique(state_names_abbr$state))]

# what state is missing a 2005 that must be KR?
state_KR <- demog_by_states %.>% 
  filter(., state == "KR")

correct_state_names <- state_names_abbr$state

funky_state <- sapply(correct_state_names, function(x) {
  
  year_vec <- demog_by_states %.>% 
    filter(., state != "KR" & state == x) %.>%
    select(., year) %.>%
    unlist(.) %.>%
    unname(.)
  
  if(2005 %in% year_vec == TRUE) {
    res <- "not it"
  } else {
    res <- "it"
  }
  
  res
  
}
)

which(funky_state == "it")

# Even without state KR all te states have 2005! -  must be a typo - 
# Getting rid of KR and adding full state names and abbrevation 

demog_by_states_clean <- demog_by_states %.>% 
  filter(., 
         state != "KR") %.>%
  left_join(., 
            state_names_abbr, 
            by = "state") %.>%
  transmute(., 
            Year = year, 
            Abbr = state,
            State = full, 
            Population = population) 

# adding more information about the longitude and latitude
# some information about the map of us
us_states_long_lat <- us_map(region = "states") %.>%  
  transmute(., 
            State = full, Abbr = abbr, 
            X = x, Y = y,
            Group = group, Order = order, Piece = piece, Hole = hole, FIPS = fips) 

unique(demog_by_states_clean$State)[unique(demog_by_states_clean$State)%nin%unique(us_states_long_lat$State)] 

demog_by_states_annually <- demog_by_states_clean %.>%
  left_join(., 
            us_states_long_lat, 
            by = c("State", "Abbr"))

# which states present in the location data are not present in the incidence data?
state_to_filter <- unique(demog_by_states_annually$State)[unique(demog_by_states_annually$State) %nin% 
                                                            unique(complete_mumps_by_states_annual$State)]

clean_demog_by_states_annually <- demog_by_states_annually %.>%
  filter(., State != state_to_filter)

time_to_filter <- min(demog_by_states_annually$Year)

# merge demography and case data 

mumps_demog_geog_annual <- complete_mumps_by_states_annual %.>%
  filter(., 
         Year > time_to_filter-1) %.>%
  left_join(.,
            clean_demog_by_states_annually, 
            by = c("State", "Year")) %.>%
  mutate(., Incidence = (Cases/Population)*1e5)


# trim the weekly data for fitting and anlysis 
mumps_inc_data_weekly <- complete_mumps_us_weekly %.>% 
  filter(., 
         Date > ('1967-12-23')) %.>% 
  mutate(., 
         # define a time variable embrace time period of the covariates
         Year = ((as.numeric(difftime(Date, min(Date), units = "weeks"))-1)/52 + 1967)
         )

if(FALSE) {  
# Only to be carried out if trajectory matching is going to be performed   
# plot and and save incidence in the folder- prcessed data for further analysis
mumps_inc_data_weekly %.>%   
  filter(., 
         Year > 1967) %.>%   
  ggplot(., 
         aes(x = Date, y = sqrt(Cases)))+
  labs(y = expression(sqrt(Cases)), 
       x = "Time (Weeks)")+  
  geom_line(size = 0.8) +
  project_theme


save(mumps_inc_data_weekly, file = "../processed_data/mumps_inc_data_weekly.Rdata")
}









  

