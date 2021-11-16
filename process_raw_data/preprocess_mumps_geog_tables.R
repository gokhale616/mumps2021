# load data from the mmwr tables to see if the data changes 
path <- "../raw_data/mumps_curated_mmwr_tabs/"

use_path <- function(file) {
  paste0(path, file)
}


data_2000_2004 <- (
  "data_2000_2004.xlsx" %.>% 
    use_path(.) %.>% 
    read_xlsx(., na = "NA") %.>% 
    select(., 2:4)
)

# checking for duplicate entries
# dup_2000_2004 <- data_2000_2004 %.>% group_by(., year, area) %.>% mutate(., n = n()) %.>% filter(., n > 1)

data_2005_2007 <- (
  "data_2005_2007.xlsx" %.>% 
    use_path(.) %.>% 
    read_xlsx(., na = "NA") %.>% 
    select(., 2:4)
)

# dup_2005_2007 <- data_2005_2007 %.>% group_by(., year, area) %.>% mutate(., n = n()) %.>% filter(., n > 1)

data_2008_2009 <- (
  "data_2008_2009.xlsx" %.>% 
    use_path(.) %.>%   
    read_xlsx(., na = "N") %.>% 
      select(., 2:4)
)

# dup_2008_2009 <- data_2008_2009 %.>% group_by(., year, area) %.>% mutate(., n = n()) %.>% filter(., n > 1)

data_2010_2015 <- (
  "data_2010_2015.xlsx" %.>% 
    use_path %.>% 
    read_xlsx(., na = "N") %.>% 
      select(., 2:4) %.>% 
      mutate(., 
             reported_cases = as.numeric(reported_cases))
)

# dup_2010_2015 <- data_2010_2015 %.>% group_by(., year, area) %.>% mutate(., n = n()) %.>% filter(., n > 1)

data_2017_2018 <- (
  "data_2017_2018.xlsx" %.>% 
    use_path(.) %.>% 
    read_xlsx(., na = "N") %.>% 
      select(., 2:4) %.>% 
      mutate(., 
             reported_cases = as.numeric(reported_cases))
)

# dup_2017_2018 <- data_2017_2018 %.>% group_by(., year, area) %.>% mutate(., n = n()) %.>% filter(., n > 1)

# make a single dataframe 
mumps_geography <- (
  data_2000_2004 %.>% 
  rbind(., data_2005_2007) %.>% 
  rbind(., data_2008_2009) %.>% 
  rbind(., data_2010_2015) %.>% 
  rbind(., data_2017_2018)
  ) 
  

clean_mmwr_table_data <- (
  mumps_geography %.>% 
  mutate(., 
         area = ifelse(area == "New York (upstate)", "New York (Upstate)", area)) %.>% 
  spread(., key = area, value = reported_cases) %.>%   
  # pivot_wider(., names_from = area, values_from = reported_cases, values_fn = mean) %.>% 
  mutate(., 
         `New York` = `New York City` + `New York (Upstate)`) %.>% 
  select(., -c(`New York City`, `New York (Upstate)`)) %.>% 
  gather(., key = "area", value = "reported_cases", -c(year)) %.>% 
  filter(., area %in% state.name) %.>% 
  transmute(., 
            Year = year, 
            State = area, 
            Cases = reported_cases) %.>% 
  arrange(., State)
  )
  
# get rid of the tycho data after 2000 and use the mmwr data 
complete_mumps_by_states_annual_clean <- (
  complete_mumps_by_states_annual %.>% 
  filter(., Year < 2000) %.>% 
  bind_rows(., clean_mmwr_table_data) %.>% 
  arrange(., State) 
  ) 
  
# combine the data with the demography data to calculate the rate of incidence
mumps_demog_geog_annual <- (
  complete_mumps_by_states_annual_clean %.>%
  filter(., 
         Year > time_to_filter-1) %.>%
  left_join(.,
            clean_demog_by_states_annually, 
            by = c("State", "Year")) %.>%
  mutate(., Incidence = (Cases/Population)*1e5)
  )


save(mumps_demog_geog_annual, file = "../processed_data/mumps_demog_geog_annual.rds")


message("NOTE :: data objects saved in '../processed_data'. To see intermediate objects run scripts individually.")

rm(list = ls())    
    
    







