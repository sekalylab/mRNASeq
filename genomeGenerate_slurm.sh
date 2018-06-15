#!/bin/bash

# user email address
#SBATCH --mail-user=EMAIL

# mail is sent to you when the job starts and when it terminates or aborts
#SBATCH --mail-type=END,FAIL

# name of job
#SBATCH --job-name=genomeGenerate

# standard output file
#SBATCH --output=genomeGenerate.log

# number of nodes and processors, memory required
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb

# time requirements
#SBATCH --time=12:00:00

# launch executable script
compress=false
while getopts d:g:m:c option
do
    case "$option" in
        d) dirData=$OPTARG;;
        g) genome=$OPTARG;;
	m) mateLength=$OPTARG;;
	c) compress=true;;
    esac
done

cmd="bash genomeGenerate.sh -d $dirData -g $genome"

if [[ ! -z $mateLength ]]
then
    cmd="$cmd -m $mateLength"
fi

if $compress
then
    cmd="$cmd -c"
fi

# echo $cmd
eval $cmd

