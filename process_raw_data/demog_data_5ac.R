##############################################################################################################
######################################## era: 1969-2017 ######################################################
##############################################################################################################
#naming the column types
ct_values <- list(col_double(), col_character(), col_double(), col_double(), col_double(), col_double(),
                  col_double(), col_double(), col_double(), col_double())

#naming the coulmns for the tibble
cn_values <- c("year", "state", "state_fips_code", "county_fips_code",
               "registry", "race", "origin", "sex", "age", "population")

#declaring the position for every variable
cs_values <- c(1, 5, 7, 9, 12, 14, 15, 16, 17, 19)

ce_values <- c(4, 6, 8, 11, 13, 14, 15, 16, 18, 26)

demog_data_1969_2017 <- (
  read_fwf(file = "../raw_data/us_1969_2017.txt", 
           fwf_positions(start = cs_values, end = ce_values,col_names = cn_values),
           col_types = ct_values,
           progress = show_progress()
           )
  ) 


# summarize for each age group
age_demog_data_1969_2017 <- (
  demog_data_1969_2017 %.>% 
  select(., 
         year, age, population) %.>% 
  group_by(., 
           year, age) %.>% 
  summarise(., 
            population = sum(population)) %.>% 
  ungroup(.)  
  ) 

# spread the dataframe from above 
wide_demog_data_1969_2017 <- (
  age_demog_data_1969_2017 %.>% 
  spread(., 
         key = age, value = population)
) 


# combine data to form five age classes 
wide_dd_5ac_1969_2017 <- bind_cols(wide_demog_data_1969_2017[,1], 
                                   tibble(rowSums(wide_demog_data_1969_2017[,2:6])), 
                                   tibble(rowSums(wide_demog_data_1969_2017[,7:16])), 
                                   tibble(rowSums(wide_demog_data_1969_2017[,17:26])),
                                   tibble(rowSums(wide_demog_data_1969_2017[,27:41])),
                                   tibble(rowSums(
                                      wide_demog_data_1969_2017[,42:ncol(wide_demog_data_1969_2017)])))


colnames(wide_dd_5ac_1969_2017) <- c("year", age_names) 

##############################################################################################################
######################################## era: 1900-1990 ######################################################
##############################################################################################################
# this function makes one dataframe from all the sheets
make_demog_data <- function(counter, sheet_v, file_v, range_v) {
  
  #collecting names of the columns
  cn_vec <- c("age", "pop1", "male1", "female1", "pop2", "male2", "female2", "pop3", "male3", "female3")
  
  read_excel(file_v, sheet = sheet_v[counter], col_names = cn_vec, range = range_v) %>%
    select(c(age, pop1, pop2, pop3)) %>%
    transmute(year = as.numeric(sheet_v[counter]),
              age = age,
              population = pop1) -> demog_data
  
  return(demog_data)
}

# file path mismatch -- find where these files are located before execution  
# Although, final map function works!

##############################################################################################################
########################## make one data from the rest of the time points in the era #########################
##############################################################################################################
#
make_1demog_data <- function(counter_f, fr_ls) {
  # browser()
  # prepare the file
  fr_ls[[1]][counter_f] %>%
    excel_sheets() -> sheets
  
  # collecting names of the columns
  cn_vec <- c("age", "pop1", "male1", "female1", "pop2", "male2", "female2", "pop3", "male3", "female3")
  
  # create a data over all spreadsheets
  map_df(1:10, make_demog_data, sheet_v = sheets, 
         file_v = fr_ls[[1]][counter_f], range_v = fr_ls[[2]][counter_f]) -> demog_data
  
  return(demog_data)
}

# excel cell decriptions for the data 
r_v <- c("A10:J85", "A10:J85", "A10:J85", "A10:J85", 
         "A10:J95", "A10:J95", "A9:J94", "A9:J94")

# list all the files in the folder 
f_v <- list.files("../raw_data/demog_data_1900_1979", full.names = TRUE)

file_list <- list(f_v, r_v)


# combine all data raw data sheets and adjust format
demog_data_1900_1979 <- (
  map_df(1:length(f_v), .f = make_1demog_data, fr_ls = file_list) %.>% 
    mutate(., age  = case_when(age == "75+"~75, 
                               age == "85+"~85, 
                               TRUE ~ as.numeric(age))
    )  
  )


# Keep this code in reserve might be useful
# demog_data_1900_1979 %.>% 
#   filter_all(., any_vars(is.na(.))) 
# 
# # make nomenclature uniform
# demog_data_1900_1979$age <-ifelse(demog_data_1900_1979$age == "75+", 75, demog_data_1900_1979$age) 
# demog_data_1900_1979$age <-ifelse(demog_data_1900_1979$age == "85+", 85, demog_data_1900_1979$age) 



# spread the data into a wide format 
wide_demog_data_1900_1979 <- demog_data_1900_1979 %.>% 
  spread(., key = age, value = population)


# making 5 age classes
wide_dd_5ac_1900_1979 <- (
  bind_cols(tibble(wide_demog_data_1900_1979[,1]), 
            tibble(rowSums(wide_demog_data_1900_1979[,2:6])), 
            tibble(rowSums(wide_demog_data_1900_1979[,7:16])), 
            tibble(rowSums(wide_demog_data_1900_1979[,17:26])),
            tibble(rowSums(wide_demog_data_1900_1979[,27:41])),
            tibble(rowSums(
            wide_demog_data_1900_1979[,42:ncol(wide_demog_data_1900_1979)], na.rm = TRUE))
            )
  ) 

colnames(wide_dd_5ac_1900_1979) <- c("year", age_names) 

#write_csv(wide_dd_10ac_1900_1979, "wide_dd_10ac_1900_1979.csv")


# demography data 1900-2017 (all age_classes)

wide_dd_1980_2017 <- wide_demog_data_1969_2017 %.>% 
  filter(., 
         year > 1979) 

wide_dd_1900_2017 <- bind_rows(wide_demog_data_1900_1979, wide_dd_1980_2017) 

# demography data 1900-2017 (5 age_classes)
wide_dd_5ac_1980_2017  <- wide_dd_5ac_1969_2017 %.>% 
  filter(., 
         year > 1979)  

wide_dd_5ac_1900_2017 <- bind_rows(wide_dd_5ac_1900_1979, wide_dd_5ac_1980_2017) 

wide_dd_5ac_1900_2017$total <- rowSums(wide_dd_5ac_1900_2017[,-1])

# change the age names for making further calulations convenient
colnames(wide_dd_5ac_1900_2017) <- c("year", sprintf("a_%d", 1:5), "total") 

if(FALSE) {
# DELETE when cleaning the codes
wide_dd_5ac_1900_2017 %.>% 
  select(., -total) %.>% 
  gather(., key = age_cohort, value = population, -year) %.>% 
  mutate(., age_cohort = as_factor(age_cohort)) %.>% 
  ggplot(., aes(x = year, y = population)) +
  geom_area(aes(fill = age_cohort)) +
  scale_y_continuous(expand = c(0.00,0)) +
  scale_x_continuous(expand = c(0.00,0)) + 
  scale_fill_viridis_d()+  
  project_theme

}

wide_dd_5ac_1900_2017 <- wide_dd_5ac_1900_2017 %.>% 
  mutate(., 
         diff_a_1 = c(NA, diff(a_1)), 
         diff_a_2 = c(NA, diff(a_2)), 
         diff_a_3 = c(NA, diff(a_3)), 
         diff_a_4 = c(NA, diff(a_4)), 
         diff_a_5 = c(NA, diff(a_5)), 
         diff_total = c(NA, diff(total))
         )  

#reading the rest of the demography data
birth_data_all <- read.csv("../raw_data/birth_data.csv")

birth_data <- birth_data_all %.>% 
  arrange(., 
          Year) %.>% 
  select(., 
         Year, Crude.Birth.Rate) %.>% 
  transmute(., 
            year = Year, 
            nu = Crude.Birth.Rate/1000) 


# filter relevant demographic data
demog_data_1909_2015 <- wide_dd_5ac_1900_2017 %>% 
  filter(year > 1908 & year < 2016)   

# combine with birth-rates
demog_data_1909_2015 <- left_join(demog_data_1909_2015, birth_data)

################################################### adding migration #########################################
# calculating mu_i using spline functions..

mu_i <- function(data = demog_data_1909_2015, 
                 t = 1910:2015, ac, diff_ac, alpha = 0, ac_1, alpha_1) {
  # browser()
  nu_splfun <- splinefun(x = data$year, 
                         y = data$nu)
  
  ac_splfun <- splinefun(x = data$year, 
                         y = as.vector(t(data[,ac])))
  
  ac_1_splfun <- splinefun(x = data$year, 
                           y = as.vector(t(data[,ac_1])))
  
  diff_ac_splfun <- splinefun(x = data$year, 
                              y = as.vector(t(data[,diff_ac])))
  #browser()
  mu_i_t <- (diff_ac_splfun(t) - alpha_1*ac_1_splfun(t-1))/ac_splfun(t) + alpha
  #browser()  
  return(c(NA, mu_i_t))  
  
}


demog_data_1909_2015 %<>% 
  mutate(mu_1 = (diff_a_1 - nu*total)/a_1 + 1/5,
         mu_2 = mu_i(ac = "a_2", ac_1 = "a_1", 
                     diff_ac = "diff_a_2", alpha = 1/10, alpha_1 = 1/5),
         mu_3 = mu_i(ac = "a_3", ac_1 = "a_2", 
                     diff_ac = "diff_a_3", alpha = 1/10, alpha_1 = 1/10),
         mu_4 = mu_i(ac = "a_4", ac_1 = "a_3", 
                     diff_ac = "diff_a_4", alpha = 1/15, alpha_1 = 1/10),
         mu_5 = mu_i(ac = "a_5", ac_1 = "a_4", 
                     diff_ac = "diff_a_5", alpha = 1/40, alpha_1 = 1/15)) 



demog_data_1909_2015 %>% 
  filter(year != 1909) -> demog_data_1910_2015

#change the names of the columns
col_names <- c(age_names, "combined")

colnames(demog_data_1910_2015)[2:7] <- col_names
colnames(demog_data_1910_2015)[8:13] <- paste0("diff_", col_names)
colnames(demog_data_1910_2015)[15:19] <- paste0("mu_", col_names[1:5])


##############################################################################################################
##################### processing the data to take care of the mortality ######################################
##############################################################################################################
# Select the migration rate "mu" inferred in the previous step
covar_mu_l <- demog_data_1910_2015 %.>%
  select(., year, starts_with("mu_")) %.>%
  gather(., key = "age_class", value = "mu", -year) 

# Break the migration rate into two rates - efflux and influx and seperate 
covar_2_mu_l <- covar_mu_l %.>%   
  rowwise(.) %.>% 
  do(., break_mu_into_2(.$mu)) %.>% 
  as_tibble(.) 

covar_mu_l <- covar_mu_l %.>% 
  bind_cols(., covar_2_mu_l)


covar_influx <- covar_mu_l %.>%
  mutate(., 
         age_class = case_when(
           age_class == "mu_[0,5)"   ~ "IN_1",
           age_class == "mu_[5,15)"  ~ "IN_2",
           age_class == "mu_[15,25)" ~ "IN_3",
           age_class == "mu_[25,40)" ~ "IN_4",
           age_class == "mu_>40"     ~ "IN_5")
         ) %.>% 
  select(., 
         -c(mu, efflux)) %.>% 
  spread(., key = age_class, value = influx)


covar_efflux <- covar_mu_l %.>%
  mutate(., 
         age_class = case_when(
           age_class == "mu_[0,5)"   ~ "OUT_1",
           age_class == "mu_[5,15)"  ~ "OUT_2",
           age_class == "mu_[15,25)" ~ "OUT_3",
           age_class == "mu_[25,40)" ~ "OUT_4",
           age_class == "mu_>40"     ~ "OUT_5")) %.>% 
  select(., -c(mu, influx)) %.>% 
  spread(., key = age_class, value = efflux)

# collect all the in a single tibble, rename with manageable column names 
colnames(demog_data_1910_2015)

demog_data <- demog_data_1910_2015 %.>%  
  select(.,
         year, nu, `[0,5)`,`[5,15)`, `[15,25)`, `[25,40)`,`>40`, combined, 
         starts_with("mu_")) %.>%
  transmute(., 
            year = year, 
            N_1 = `[0,5)`, N_2 = `[5,15)`, N_3 = `[15,25)`, N_4 = `[25,40)`, N_5 = `>40`, 
            MU_1 = `mu_[0,5)`, MU_2 = `mu_[5,15)`, MU_3 = `mu_[15,25)`, 
            MU_4 = `mu_[25,40)`, MU_5 = `mu_>40`,   
            Births = combined*nu) %>% 
  right_join(covar_influx, by = "year") %>% 
  right_join(covar_efflux, by = "year") %>% 
  dplyr::slice(1:n(), rep(n(), times = 3)) %>% 
  select(-year) %>% 
  mutate(Year = seq(1910, 2018), 
         OUT_1 = -OUT_1, OUT_2 = -OUT_2, OUT_3 = -OUT_3,
         OUT_4 = -OUT_4, OUT_5 = -OUT_5) 


##############################################################################################################
############ combine the  vaccine data from the latest 5 class co-variate data ###############################
##############################################################################################################

vacc_data <- read_csv(file = "../raw_data/nat_vac_cover.csv")

demog_vacc_data <- demog_data %>% 
  right_join(vacc_data) %>% 
  mutate(year = seq(1910, 2018)) 


##############################################################################################################
######################### combine the observation process covariate ##########################################
##############################################################################################################

load(file = "../processed_data/mumps_case_reports.rds")


# covariate data #
# inferred probabilities of recording an age class (from the data)
eta_a_data <- mumps_case_reports %.>% 
  mutate(.,
         eta_a = 1-unknown/total) %>%
  select(., 
         year, eta_a) 

# until you read Mathieu's manuscript this hacky imputation of rho_age - hold the first value constant
# rho_age_data %>% slice(1)

hacky_eta_a_imputation <- tibble(year = c(1995, 2015), 
                                 eta_a = c(0.9840213, 0.9991823))

eta_a_data <- (
  eta_a_data %.>% 
  filter(., 
         year %nin% c(1995, 2015)) %.>% 
  bind_rows(., 
            hacky_eta_a_imputation)%.>% 
  arrange(., 
          year) %.>% 
  slice(., 
        rep(1, times = 67), 2:n()) %.>%
  mutate(., 
         year = seq(1910, 2018)) %.>% 
  mutate(., 
         eta_a = ifelse(year < 1977, NA, eta_a))
  )
  

# combine into one covar tibble
mumps_covariates <- demog_vacc_data %.>% 
  full_join(., 
            eta_a_data, by = "year")


save(mumps_covariates, file = "../processed_data/mumps_covariates.rds")



##############################################################################################################
################################### configuring the contact matrix ###########################################
##############################################################################################################
# These steps caliberate the UK polymod matrix to the US populationa and correct it for reciprocity 

polymod_UK_5 <- contact_matrix(survey = polymod, countries = "United Kingdom",
                               age.limits = c(0, 5, 15, 25, 40))

#extracting population from 2005

US_pop_2005 <- demog_data %.>%
  filter(., 
         Year == 2005) %.>% 
  select(., 
         starts_with("N_")) %.>% 
  as.numeric(.) 


#total number of contacts scaled to the population of the US
E_polymod_US <- matrix(NA, nrow = nrow(polymod_UK_5$matrix), ncol = ncol(polymod_UK_5$matrix))

for(i in 1:nrow(polymod_UK_5$matrix)) {
  E_polymod_US[i, ] <- US_pop_2005[i]*polymod_UK_5$matrix[i,]} 

dimnames(E_polymod_US) <- list(age_names, age_names)

#making the matrix symmtric:
Es_polymod_US <- 0.5*(E_polymod_US + t(E_polymod_US))

#Average daily contacts for US:
Cs_polymod_US_5 <- matrix(NA, nrow = nrow(polymod_UK_5$matrix), ncol = ncol(polymod_UK_5$matrix))

for(i in 1:nrow(polymod_UK_5$matrix)) {
  Cs_polymod_US_5[i, ] <- 1/US_pop_2005[i]* Es_polymod_US[i,]} 

dimnames(Cs_polymod_US_5) <- list(age_names, age_names)

contact_matrix <- Cs_polymod_US_5

save(contact_matrix, file = "../processed_data/contact_matrix.rds")

#loading both contact matrices:
# load("cm_18_ac.Rdata")
# load("UScontact_matrix.Rdata")

if(FALSE) {
##############################################################################################################
########################################## Configuring Serology Data #########################################
##############################################################################################################

# laod serology data 
read_csv("../raw_data/england_1987_mumps_serology_edmunds.csv", 
         col_types = cols_only(Age = col_double(), 
                               prop_sero_pos = col_double())) -> sero_data

# correct ages 
sero_data %>% 
  mutate(Age = seq(1, 39, by = 1)) %>% 
  add_row(Age = c(0, seq(40, 85)), 
          prop_sero_pos = c(rep(min(.$prop_sero_pos)), 
                            rep(max(.$prop_sero_pos), length(seq(40, 85))))) %>% 
  arrange(Age) -> sero_data_age_correction

# cumlative sero-positivity curve
if(FALSE) {
sero_data_age_correction %>%
  mutate(extrapolation = ifelse(Age < 1 | Age > 39, TRUE, FALSE)) %>%
  ggplot(aes(x = Age, y = prop_sero_pos)) +
  geom_line(aes(colour = extrapolation, group = 1), size = 0.8) +
  geom_point(aes(colour = extrapolation), size = 2.5) +
  scale_colour_manual(name = "Extrapolation", values = c("grey30", "#FFD92F")) +
  labs(y = "Proportion Sero-positive") +
  project_theme +
  theme(legend.position = c(0.8, 0.2))
}  

# Use the demog data for all age classes for resetting seropoitive proportion into 5 age classes 
wide_dd_1900_2017 %>% 
  filter(year < 1968) %>% 
  summarise_all(mean) %>% 
  select(-year) %>% 
  gather(key = "Age", value = "Population") -> mean_pop_prevac_US

if(FALSE) {
mean_pop_prevac_US %>% 
  ggplot(aes(x = Age, y = Population)) +
  geom_point() + 
  labs(y = expression(log[10]~(Population))) +
  scale_y_continuous(trans = "log10") +
  project_theme
}

# Combine the two data streams
mean_pop_prevac_US %>% 
  mutate(Age = seq(0, 85)) %>% 
  right_join(sero_data_age_correction, by = "Age") %>% 
  mutate(SeroPosPop = (Population*prop_sero_pos)) %>% 
  select(Age, Population, SeroPosPop) -> MeanPopSeroPos

# make five age groups
age_names
#[0,5)
MeanPopSeroPos %>% 
  filter(Age < 5) %>% 
  summarise_all(sum) %>% 
  mutate(Age = age_names[1]) -> SeroPosAgeClass1

#[5,15)
MeanPopSeroPos %>% 
  filter(Age > 4 & Age < 15) %>% 
  summarise_all(sum) %>% 
  mutate(Age = age_names[2]) -> SeroPosAgeClass2  

#[15,25)
MeanPopSeroPos %>% 
  filter(Age > 14 & Age < 25) %>% 
  summarise_all(sum) %>% 
  mutate(Age = age_names[3]) -> SeroPosAgeClass3

#[25,40)
MeanPopSeroPos %>% 
  filter(Age > 24 & Age < 40) %>% 
  summarise_all(sum) %>% 
  mutate(Age = age_names[4]) -> SeroPosAgeClass4

#>=40
MeanPopSeroPos %>% 
  filter(Age > 39) %>% 
  summarise_all(sum) %>% 
  mutate(Age = age_names[5]) -> SeroPosAgeClass5


# bind all age_groups
SeroPosAgeClass1 %>% 
  bind_rows(SeroPosAgeClass2) %>% 
  bind_rows(SeroPosAgeClass3) %>% 
  bind_rows(SeroPosAgeClass4) %>% 
  bind_rows(SeroPosAgeClass5) %>% 
  mutate(PropSeroPos = SeroPosPop/Population) %>% 
  select(Age, PropSeroPos) -> SeroPosAgeClass

if(FALSE) {
SeroPosAgeClass %>% 
  mutate(Age = factor(Age, levels = age_names)) %>%   
  ggplot(aes(x = Age, y = PropSeroPos)) +
  geom_line(aes(group = 1), size = 0.8) +
  geom_point(size = 2.5) +
  labs(y = "Proportion Sero-Positives") +
  project_theme
}  
  
# form a vector to use for estimation 
SeroPosAgeClass %>% 
  select(PropSeroPos) %>% 
  unlist() -> prop_sero_pos_data

names(prop_sero_pos_data) <- age_names 

save(prop_sero_pos_data, file = "../processed_data/prop_sero_pos_data.rds")
}

message("NOTE :: data objects saved in '../processed_data'. To see intermediate objects run scripts individually.")

rm(list = ls())



















