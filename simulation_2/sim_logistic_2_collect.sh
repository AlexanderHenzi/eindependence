#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_2_collect"
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=epyc2

#### Your shell commands below this line ####

module load R
R CMD BATCH --no-save --no-restore simulation_2_collect.R

