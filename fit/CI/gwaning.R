# source confidence int pre-requisites
source("./CI_prereqs.R")

# generate bootstrapped distribution of the parameter estimates 
gwaning_bootstrap_res <- (
  mclapply(1:bs_nsim, 
           FUN = trajmatch_protocol, 
           hypo_name = "gwaning", 
           covar_name = "sigmoid",
           log_vars = "dwan", 
           logit_vars = NULL, 
           mc.cores = bs_n_cores, mc.set.seed = TRUE) %.>% 
    do.call(rbind, 
            lapply(1:bs_nsim, 
                   FUN = function(x) {.[[x]]})
    )
)

# save the results
save(gwaning_bootstrap_res, file = CI_res_path_fn("gwaning_bootstrap_res.rds"))