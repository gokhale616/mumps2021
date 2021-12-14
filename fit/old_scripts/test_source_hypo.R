# this script does a quick source of all hypotheses to make sure that files run 
# and store objects in the right location.
# NOTE:: Make sure that DEoptim convergence parameters are calibrated for a quick run
library(wrapr)
list.files(".", 
           pattern = "waning") %.>% 
  lapply(., source)

list.files(".", 
           pattern = "leaky2") %.>% 
  lapply(., source)

