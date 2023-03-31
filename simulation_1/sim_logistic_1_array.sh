#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_1"
#SBATCH --output="outfiles/array_%j.out"
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=600M
#SBATCH --partition=epyc2
#SBATCH --array=1-9680

#### Your shell commands below this line ####

module load R
Rscript simulation_1/simulation_1.R $q $corr $pen

