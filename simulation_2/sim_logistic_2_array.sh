#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_2"
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=700M
#SBATCH --partition=epyc2
#SBATCH --array=1-880

#### Your shell commands below this line ####

module load R
Rscript simulation_2/simulation_2.R $q $corr

