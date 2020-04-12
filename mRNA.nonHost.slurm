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
#SBATCH --tasks=1
#SBATCH --mem=256gb
#SBATCH --partition=smp

# time requirements
#SBATCH --time=72:00:00

# dependencies
#SBATCH --depend=afterok:SLURM_JOB_ID

# create array
#SBATCH --array=1-BATCH

# initialize directories
pairEnd=false
while getopts d:p option
do
    case "$option" in
        d) dirData=$OPTARG;;
	p) pairEnd=true;;
    esac
done

# launch executable script
cmd="bash mRNA.nonHost_seq.sh -d $dirData/raw$SLURM_ARRAY_TASK_ID"

if $pairEnd
then
    cmd="$cmd -p"
fi

# echo $cmd
eval $cmd
