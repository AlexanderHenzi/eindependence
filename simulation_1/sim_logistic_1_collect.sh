#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_1_collect"
#SBATCH --output="outfiles/collect_%j.out"
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=epyc2

#### Your shell commands below this line ####

module load R
Rscript simulation_1/simulation_1_collect.R $name

