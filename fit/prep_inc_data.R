# this subset will be used to fit the model to the data
in_sample_mumps_case_reports <- (
  mumps_case_reports %.>% 
    filter(., year < 2013)
)
