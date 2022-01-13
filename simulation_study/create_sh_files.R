# This function prepares the .sh scripts for submitting estimation scripts to the cluster
create_sh <- function(pattern = "hypo") {
  paste0("#!/bin/bash\n
#SBATCH --job-name=", pattern, "              # Job name
#SBATCH --partition=rohani_p                 # Partition (queue) name
#SBATCH --nodes=1                            # Number of nodes
#SBATCH --ntasks=24                          # Number of cores	
#SBATCH --mem=16gb                            # Job memory request
#SBATCH --time=240:00:00                      # Time limit hrs:min:sec
#SBATCH --output=", pattern, ".%j.out          # Standard output log
#SBATCH --error=", pattern, ".%j.err          # Standard error log          
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dg34432@uga.edu   # Where to send mail	
          
cd $SLURM_SUBMIT_DIR
          
module load R/4.0.0-foss-2019b

R CMD BATCH ", pattern, ".R
")
}  

write_sh <- function(pattern = "hypo" , path) {
  
  fileConn <- file(paste0(path, pattern, ".sh"))
  writeLines(create_sh(pattern = pattern), fileConn)
  close(fileConn)
  
}

# name all the hypotheses - heterogenous mixing hypothesis

n_hypotheses <- c("leaky_vs_waning")

lapply(n_hypotheses, function(x){write_sh(pattern = x, path = "./")})



rm(list = ls())






  