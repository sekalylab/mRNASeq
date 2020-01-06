#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=mrna.homolog_%a

# standard output file
#SBATCH --output=mrna.homolog_%a.out

# number of nodes and processors, memory required
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --partition=smp

# dependencies
#SBATCH --depend=afterok:SLURM_JOB_ID

# time requirements
#SBATCH --time=72:00:00

# create array
#SBATCH --array=1-BATCH

# initialize directories
compress=false
pairEnd=false
while getopts d:g: option
do
    case "$option" in
        d) dirData=$OPTARG;;
        g) genome=$OPTARG;;
    esac
done

# launch executable script
cmd="bash mRNA.homolog_seq.sh -d $dirData/raw$SLURM_ARRAY_TASK_ID"
cmd="$cmd -g $genome"
# echo $cmd
eval $cmd
