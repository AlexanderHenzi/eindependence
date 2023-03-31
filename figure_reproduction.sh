#!/bin/bash

# this script reproduces all the figures in the paper and supplementary material
# note: if the array size in sim_logistic_1_array is much smaller than the maximum number of nodes/cores on the cluster,
# there is a lot to be gained from optimizing this script, as all the array jobs currently run consecutively.
#
#

mkdir -p outfiles

# simulation 1 (first job runs simulations on an array, second job collects the data)
sbatch -W --job-name="sim_logistic_1_1" --parsable --export=q=4,corr="pos",pen="FALSE" simulation_1/sim_logistic_1_array.sh
sbatch --job-name="sim_logistic_1_collect_1" --parsable -W --export=name="simulation_1_4posfalse" simulation_1/sim_logistic_1_collect.sh

# repeat for different parameter values
sbatch --job-name="sim_logistic_1_2" -W --export=q=8,corr="pos",pen="FALSE" --parsable simulation_1/sim_logistic_1_array.sh
sbatch --job-name="sim_logistic_1_collect_2" --parsable -W --export=name="simulation_1_8posfalse" simulation_1/sim_logistic_1_collect.sh

sbatch --job-name="sim_logistic_1_3" -W --parsable --export=q=4,corr="negative_cor",pen="FALSE" simulation_1/sim_logistic_1_array.sh
sbatch --job-name="sim_logistic_1_collect_3" --parsable -W --export=name="simulation_1_4negfalse" simulation_1/sim_logistic_1_collect.sh

# produces Table 1, Figures 1,2,S4
sbatch --job-name="sim_logistic_1_plot_1" --parsable --export=eps=0.05,path="simulation_1_4posfalse.rda" simulation_1/sim_logistic_1_plot.sh
# produces Figure S1
sbatch --job-name="sim_logistic_1_plot_2" --parsable --export=eps=0,path="simulation_1_4posfalse.rda" simulation_1/sim_logistic_1_plot.sh
# produces Figures S5, S7
sbatch --job-name="sim_logistic_1_plot_3" --parsable --export=eps=0.05,path="simulation_1_8posfalse.rda" simulation_1/sim_logistic_1_plot.sh
# produces Figures S6, S8
sbatch --job-name="sim_logistic_1_plot_4" --parsable --export=eps=0.05,path="simulation_1_4negfalse.rda" simulation_1/sim_logistic_1_plot.sh

# simulation 2 (workflow similar to simulation 1) 
sbatch -W --job-name="sim_logistic_2_1" --parsable --export=q=4,corr="pos" simulation_2/sim_logistic_2_array.sh
sbatch --job-name="sim_logistic_2_collect_1" --parsable -W --export=name="simulation_2_4pos" simulation_2/sim_logistic_2_collect.sh

sbatch -W --job-name="sim_logistic_2_2" --parsable --export=q=8,corr="pos" simulation_2/sim_logistic_2_array.sh
sbatch --job-name="sim_logistic_2_collect_2" --parsable -W --export=name="simulation_2_8pos" simulation_2/sim_logistic_2_collect.sh

sbatch -W --job-name="sim_logistic_2_3" --parsable --export=q=4,corr="negative_cor" simulation_2/sim_logistic_2_array.sh
sbatch --job-name="sim_logistic_2_collect_3" --parsable -W --export=name="simulation_2_4neg" simulation_2/sim_logistic_2_collect.sh

# produces Figure S3
sbatch --job-name="sim_logistic_2_plot_1" --parsable ==export=path="simulation_2_4pos" simulation_2/sim_logistic_2_plot.sh
# produces Figure S9
sbatch --job-name="sim_logistic_2_plot_1" --parsable ==export=path="simulation_2_8pos" simulation_2/sim_logistic_2_plot.sh
# produces Figure S10
sbatch --job-name="sim_logistic_2_plot_1" --parsable ==export=path="simulation_2_4neg" simulation_2/sim_logistic_2_plot.sh

# comment this line out for troubleshooting
rm -r outfiles
