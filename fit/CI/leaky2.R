# source confidence int pre-requisites
source("./CI_prereqs.R")

# modify the defaults to make introduction in the [15,25) year age class
mod_bs_param_vals_est <- param_vals_est
mod_bs_param_vals_est["p_intro"] <- 3

# generate bootstrapped distribution of the parameter estimates 
leaky2_bootstrap_res <- (
  mclapply(1:bs_nsim, 
           FUN = trajmatch_protocol, 
           hypo_name = "leaky2", 
           covar_name = "constant",
           params = mod_bs_param_vals_est, 
           log_vars = "t_intro", 
           logit_vars = "epsilon2", 
           mc.cores = bs_n_cores, mc.set.seed = TRUE) %.>% 
    do.call(rbind, 
            lapply(1:bs_nsim, 
                   FUN = function(x) {.[[x]]})
    )
)

# save the results 
save(leaky2_bootstrap_res, file = CI_res_path_fn(fn = "leaky2_bootstrap_res.rds"))
