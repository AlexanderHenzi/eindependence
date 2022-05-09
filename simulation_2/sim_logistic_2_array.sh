#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_2"
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=700M
#SBATCH --partition=epyc2
#SBATCH --array=1-880

#### Your shell commands below this line ####

module load R
R CMD BATCH --no-save --no-restore simulation_2.R

