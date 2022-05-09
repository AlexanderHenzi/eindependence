#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_1"
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=600M
#SBATCH --partition=epyc2
#SBATCH --array=1-9680

#### Your shell commands below this line ####

module load R
R CMD BATCH --no-save --no-restore simulation_1.R

