#!/bin/bash

#SBATCH --job-name=waning_slow_sigmoid              # Job name
#SBATCH --partition=rohani_p                 # Partition (queue) name
#SBATCH --nodes=1                            # Number of nodes
#SBATCH --ntasks=24                          # Number of cores	
#SBATCH --mem=16gb                            # Job memory request
#SBATCH --time=500:00:00                      # Time limit hrs:min:sec
#SBATCH --output=waning_slow_sigmoid.%j.out          # Standard output log
#SBATCH --error=waning_slow_sigmoid.%j.err          # Standard error log          
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dg34432@uga.edu   # Where to send mail	
          
cd $SLURM_SUBMIT_DIR
          
module load R/4.0.0-foss-2019b

R CMD BATCH waning_slow_sigmoid.R

