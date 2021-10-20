# This function prepares the .sh scripts for submitting estimation scripts to the cluster
create_sh <- function(pattern = "hypothesis") {
  paste0("#!/bin/bash\n
#SBATCH --job-name=", pattern, "              # Job name
#SBATCH --partition=rohani_p                 # Partition (queue) name
#SBATCH --nodes=1                            # Number of nodes
#SBATCH --ntasks=24                          # Number of cores	
#SBATCH --mem=8gb                            # Job memory request
#SBATCH --time=240:00:00                      # Time limit hrs:min:sec
#SBATCH --output=", pattern, ".%j.out          # Standard output log
#SBATCH --error=", pattern, ".%j.err          # Standard error log          
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dg34432@uga.edu   # Where to send mail	
          
cd $SLURM_SUBMIT_DIR
          
module load R/3.6.2-foss-2019b
          
R CMD BATCH ", pattern, ".R
")
}  

write_sh <- function(pattern = "hypothesis" , path) {
  
  fileConn <- file(paste0(path, pattern, ".sh"))
  writeLines(create_sh(pattern = pattern), fileConn)
  close(fileConn)
  
}

# name all the hypotheses - heterogenous mixing hypothesis

n_hypotheses <- c("waning_lgs", sprintf("leaky2_%d_lgs", 0:4), 
                  sprintf("waning_leaky2_%d_lgs", 0:4))

lapply(n_hypotheses, function(x){write_sh(pattern = x, path = "./unknwn_popn_all_hypo/")})


# name all the hypotheses - homogenous mixing hypothesis

n_hypotheses_homo <- c("waning_lgs", "waning_leaky2_lgs", "leaky2_lgs")

lapply(n_hypotheses_homo, function(x){write_sh(pattern = x, path = "./unknwn_popn_homo_all_hypo/")})








  