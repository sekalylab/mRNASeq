#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=mrna.nonhost_%a

# standard output file
#SBATCH --output=mrna.nonhost_%a.out

# number of nodes and processors, memory required
#SBATCH --tasks=24
#SBATCH --mem=120gb
#SBATCH --partition=smp

# time requirements
#SBATCH --time=72:00:00

# dependencies
#SBATCH --depend=afterok:$SLURM_JOB_ID

# create array
#SBATCH --array=1-BATCH

# initialize directories
while getopts d: option
do
    case "$option" in
        d) dirData=$OPTARG;;
    esac
done

# launch executable script
cmd="bash mRNA.nonHost_seq.sh -d $dirData/raw$SLURM_ARRAY_TASK_ID"
# echo $cmd
eval $cmd
