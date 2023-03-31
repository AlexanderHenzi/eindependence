#!/bin/bash

# Environment variables

# First job - no dependencies
# option 'parsable' formats job id to number (required!)
jid1=$(sbatch --job-name="sim_logistic_1" --parsable simulation_1/sim_logistic_1_array.sh)

# Second job collecting the results of first

jid2=$(sbatch --job-name="sim_logistic_1_collect" --dependency=afterany:$jid1 --export= simulation_1/sim_logistic_1_collect.sh)

sbatch --job-name="sim_logistic_1_plot" --dependency=afterany:$jid2 --parsable --export=eps=0.05,path="simulation_1.rda" simulation_1/sim_logistic_1_plot.sh

