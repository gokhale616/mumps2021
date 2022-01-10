# source confidence int pre-requisites
source("./CI_prereqs.R")

# generate bootstrapped distribution of the parameter estimates 
noloss_bootstrap_res <- (
  mclapply(X = 1:bs_nsim, 
           FUN = trajmatch_protocol, 
           hypo_name = "noloss", 
           covar_name = "constant",
           log_vars = NULL, 
           logit_vars = NULL, 
           mc.cores = bs_n_cores, mc.set.seed = TRUE) %.>% 
    do.call(rbind, 
            lapply(1:bs_nsim, 
                   FUN = function(x) {.[[x]]})
    )
)

# save the results 
save(noloss_bootstrap_res, file = CI_res_path_fn("noloss_bootstrap_res.rds"))
