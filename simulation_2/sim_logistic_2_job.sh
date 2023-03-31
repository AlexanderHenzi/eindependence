#!/bin/bash

# Environment variables

# First job - no dependencies
# option 'parsable' formats job id to number (required!)
jid1=$(sbatch --job-name="sim_logistic_2" --parsable simulation_2/sim_logistic_2_array.sh)

# Second job collecting the results of first

jid2=$(sbatch --job-name="sim_logistic_2_collect" --dependency=afterany:$jid1 simulation_2/sim_logistic_2_collect.sh)

sbatch --job-name="sim_logistic_2_plot" --dependency=afterany:$jid2 --parsable --export=eps=path="simulation_2.rda" simulation_2/sim_logistic_2_plot.sh

