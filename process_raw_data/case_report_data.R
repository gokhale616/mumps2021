# load raw data 
sheet_name <- excel_sheets("../raw_data/mumps_case_data.xlsx")
  
mumps_1977_1994 <- read_excel("../raw_data/mumps_case_data.xlsx", sheet = "Sheet1")
mumps_1995 <- read_excel("../raw_data/mumps_case_data.xlsx", sheet = "Sheet2")
mumps_1996_2018 <- read_excel("../raw_data/mumps_case_data.xlsx", sheet = "Sheet3")

mumps_1996_2018 %<>% na_if(1.1) 
if(FALSE) {
# plotting time series #######################################################################################
## DELETE this when cleaning
mumps_1977_1994 %>% 
  gather(key = "Age Classes", value = "Cases", -c(Year, Region)) -> mumps_1977_1994_l

mumps_1977_1994_l %>% 
  select(`Age Classes`) %>% 
  unique() %>% t() %>% as.vector() -> levels_1977_1994

mumps_1977_1994_l %>% 
  filter(`Age Classes` != "total") %>% 
  ggplot(aes(x = Year, y = sqrt(Cases))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = c(seq(1977, 1995, by = 4))) +
  labs(y = expression(sqrt(Cases))) +
  facet_wrap(~factor(`Age Classes`, levels = levels_1977_1994)) +
  theme(aspect.ratio = 0.5)

mumps_1996_2018 %>% 
  gather(key = "Age Classes", value = "Cases", -c(Year, region)) -> mumps_1996_2018_l

mumps_1996_2018_l %>% 
  select(`Age Classes`) %>% 
  unique() %>% t() %>% as.vector() -> levels_1996_2018

mumps_1996_2018_l %>% 
  filter(`Age Classes` != "total") %>% 
  ggplot(aes(x = Year, y = sqrt(Cases))) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = c(seq(1996, 2016, by = 4))) +
  labs(y = expression(sqrt(Cases))) +
  facet_wrap(~factor(`Age Classes`, levels = levels_1996_2018)) +
  theme(aspect.ratio = 0.5) 

mumps_1995 %>% 
  gather(key = "Age Classes", value = "Cases", -c(Year, region)) -> mumps_1995_l

mumps_1995_l %>% 
  select(`Age Classes`) %>% 
  unique() %>% t() %>% as.vector() -> levels_1995

}
# condense the raw data to 6 classes #########################################################################
mumps_1977_1994_proc <- (
  mumps_1977_1994 %.>% 
  mutate(.,
         `[0,5)`   = `[0,1)`  + `[1,5)`, 
         `[5,15)`  = `[5,10)` + `[10,15)`, 
         `[15,25)` = `[15,20)`  + `[20,25)`,
         `[25,40)` = `[20,25)` + `[25,30)` + `[30,40)`,
         `>40`     = `[40,50)` + `[50,60)`  + `>60`, 
         Total     = total) %.>% 
  select(.,
         Year, total, `[0,5)`, `[5,15)`, `[15,25)`, `[25,40)`, `>40`, unknown)
  ) 

mumps_1996_2018_proc <- (
  mumps_1996_2018 %.>% 
  mutate(., 
         `[0,5)`   = `[0,1)`  + `[1,5)`, 
         `>40`     = `[40,65)` + `>65`, 
         Total     = total) %.>% 
  select(., 
         Year, total, `[0,5)`, `[5,15)`, `[15,25)`, `[25,40)`, `>40`, unknown) %.>% 
  filter(., 
         Year != 2015)
  ) 

# combine both time series ###################################################################################
mumps_1977_2018_proc <- (
  mumps_1977_1994_proc %.>%
  bind_rows(., 
            mumps_1996_2018_proc)
  ) 



# define missing data point for further analysis
missing_years <- (
  tibble(Year = c(1976, 1995, 2015), 
         total = rep(NA, 3), `[0,5)` = rep(NA, 3), `[5,15)` = rep(NA, 3),
         `[15,25)` = rep(NA, 3), `[25,40)` = rep(NA, 3), `>40` = rep(NA, 3), 
          unknown = rep(NA, 3)
         )
  )

# add missing data to data 
mumps_case_reports <- (
  mumps_1977_2018_proc %.>%
  bind_rows(., 
            missing_years) %.>%
  arrange(., 
          Year)
  )

colnames(mumps_case_reports)[1] <- "year"

save(mumps_case_reports, file = "../processed_data/mumps_case_reports.rds")














