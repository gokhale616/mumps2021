# source confidence int pre-requisites
source("./CI_prereqs.R")

# generate bootstrapped distribution of the parameter estimates 
waning_bootstrap_res <- (
  mclapply(X = 1:bs_nsim, 
           FUN = trajmatch_protocol, 
           hypo_name = "waning", 
           covar_name = "sigmoid", 
           log_vars = "dwan", 
           logit_vars = NULL, 
           mc.cores = bs_n_cores, 
           mc.set.seed = TRUE) %.>% 
    do.call(rbind, 
            lapply(1:bs_nsim, 
                   FUN = function(x) {.[[x]]})
    )
)

# save the results 
save(waning_bootstrap_res, file = CI_res_path_fn("waning_bootstrap_res.rds"))
