# source confidence int pre-requisites
source("./CI_prereq.R")

# generate bootstrapped distribution of the parameter estimates 
op_noisy_bootstrap_res <- (
  mclapply(X = 1:n_sim, 
           FUN = optim_protocol, 
           ts_data = op_noisy_mle_ts,
           mle = op_noisy_mle,
           mc.cores = detectCores(), 
           mc.set.seed = TRUE) %.>% 
    do.call(rbind, 
            lapply(1:n_sim, 
                   FUN = function(x) {.[[x]]})
    )
)

# save the results 
save(op_noisy_bootstrap_res, file = "./op_noisy_bootstrap_res.rds")
