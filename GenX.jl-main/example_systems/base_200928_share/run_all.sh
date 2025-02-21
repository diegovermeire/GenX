#!/bin/bash

source /etc/profile

module load julia/1.9.2
module load gurobi/gurobi-1000

echo "My SLURM_ARRAY_TASK_ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

julia --project=. Interannual_Variability/INFORMS_runs/base_200928/scenario_runfile.jl $LLSUB_RANK $LLSUB_SIZE
