#!/bin/bash

#SBATCH --job-name=o_noisy              # Job name
#SBATCH --partition=BATCH                  # Partition (queue) name
#SBATCH --q cider_qos
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks=32                           # Number of cores	
#SBATCH --mem=20gb                           # Job memory request
#SBATCH --time=240:00:00                      # Time limit hrs:min:sec
#SBATCH --output=waning.%j.out         # Standard output log
#SBATCH --error=waning.%j.err          # Standard error log          
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dg34432@uga.edu   # Where to send mail	
          
cd $SLURM_SUBMIT_DIR
          
module load R/4.0.0-foss-2019b

R CMD BATCH o_noisy.R

