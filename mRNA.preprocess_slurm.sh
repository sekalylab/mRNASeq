#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=mrna.pre_%a

# standard output file
#SBATCH --output=mrna.pre_%a.out

# number of nodes and processors, memory required
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# time requirements
#SBATCH --time=72:00:00

# dependencies
#SBATCH --depend=afterok:$SLURM_JOB_ID

# create array
#SBATCH --array=1-BATCH

# initialize directories
compress=false
pairEnd=false
isoform=false
homolog=false
while getopts d:g:m:a:cpio option
do
    case "$option" in
        d) dirData=$OPTARG;;
        g) genome=$OPTARG;;
	m) mateLength=$OPTARG;;
	a) adapterFile=$OPTARG;;
	c) compress=true;;
	p) pairEnd=true;;
	i) isoform=true;;
	o) homolog=true;;
    esac
done

# launch executable script
cmd="bash mRNA.preprocess_seq.sh -d $dirData/raw$SLURM_ARRAY_TASK_ID"
cmd="$cmd -g $genome"
cmd="$cmd -a $adapterFile"

if [[ ! -z $mateLength ]]
then
    cmd="$cmd -m $mateLength"
fi

if $compress
then
    cmd="$cmd -c"
fi

if $pairEnd
then
    cmd="$cmd -p"
fi

if $isoform
then
    cmd="$cmd -i"
fi

if $homolog
then
    cmd="$cmd -o"
fi

# echo $cmd
eval $cmd
