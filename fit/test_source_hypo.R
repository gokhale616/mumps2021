# this script does a quick source of all hypotheses to make sure that files run 
# and store objects in the right location.
# NOTE:: Make sure that DEoptim convergence parameters are calibrated for a quick run
library(wrapr)

test_source <- function(x){print(x); source(x); rm(list = ls())}

list.files(".", 
           pattern = "waning") %.>% 
  lapply(., test_source)

list.files(".", 
           pattern = "leaky2") %.>% 
  lapply(., test_source)

# no-loss models 
source("./constant.R")
source("./rapid.R")
source("./sigmoid.R")
source("./slow.R")

rm(list = ls())
