#!/bin/bash

# Environment variables

# First job - no dependencies
# option 'parsable' formats job id to number (required!)
jid1=$(sbatch --job-name="sim_logistic_2" --parsable sim_logistic_2_array.sh)

# Second job collecting the results of first

sbatch --job-name="sim_logistic_2_collect" --dependency=afterany:$jid1 sim_logistic_2_collect.sh


