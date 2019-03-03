#!/bin/bash
# Job name:
#SBATCH --job-name=CollinsDYE
#
# Account:
#SBATCH --account=fc_surfwill
#
# Partition:
#SBATCH --partition=savio
#
# Request one node:
#SBATCH --nodes=1
#
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=20
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=02:00:00
#
# Log file location:
#SBATCH --output="../SLURMlogs/job_%j.out"





## Command(s) to run (example):
module load r/3.5.1
module load r-packages/default

R CMD BATCH ../rScripts/CollinsDYE.R
