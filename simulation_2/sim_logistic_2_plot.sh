#!/bin/bash
#SBATCH --mail-user=alexander.henzi@stat.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="sim_logistic_2_plot"
#SBATCH --output="outfiles/plot2_%j.out"
#SBATCH --time=00:40:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=epyc2

#### Your shell commands below this line ####

module load R
Rscript simulation_2_plots.R $path

