#!/bin/bash

# Environment variables

# First job - no dependencies
# option 'parsable' formats job id to number (required!)
jid1=$(sbatch --job-name="sim_logistic_1" --parsable sim_logistic_1_array.sh)

# Second job collecting the results of first

sbatch --job-name="sim_logistic_1_collect" --dependency=afterany:$jid1 sim_logistic_1_collect.sh


