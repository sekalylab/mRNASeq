#!/bin/bash
# @author Slim Fourati (sxf279@case.edu)
# @author Aarthi Talla (axt427@case.edu)
# @version 0.5

# read input arguments
email="sxf279@case.edu"
compress=false
pairEnd=false
isoform=false
homolog=false
adapterFile="/mnt/projects/SOM_PATH_RXS745U/bin/Trimmomatic-0.36/adapters"
adapterFile="${adapterFile}/TruSeq3-SE.fa"
genome=GRCh38
acceptedGenome=("GRCh38" "Mmul_8" "MacFas5" "GRCm38")
while getopts :d:e:g:a:cpioh option
do
    case "${option}" in
	h) echo "Command: bash mRNA.preprocess_master.sh -d {fastq/directoryfastq} ..."
	    echo "argument: d=[d]irectory with raw data (required)"
	    echo "          g=reference [g]enome"
	    echo "          a=path to FASTA files with [a]dapters to be trimmed"
	    echo "          c=raw files are [c]ompressed as DSRC files"
	    echo "          p=[p]aired-end sequencing"
	    echo "          i=[i]soform transcript/exon counts"
	    echo "          o=use h[o]mology to annotate reads"  
	    echo "          e=[e]mail address"
	    echo "          h=print [h]elp"
	    exit 1;;
	d) dirFastq=$OPTARG;;
	e) email=$OPTARG;;
	g) genome=$OPTARG
	    if [[ ! "${acceptedGenome[@]}" =~ "$genome" ]]
	    then
		echo "Invalid -g argument: choose between ${acceptedGenome[@]}"
		exit 1
	    fi;;
	c) compress=true;;
	p) pairEnd=true;;
	o) homolog=true;;
	i) isoform=true;;
	a) adapterFile=$OPTARG;;
	\?) echo "Invalid option: -$OPTARG"
	    exit 1;;
	:)
	    echo "Option -$OPTARG requires an argument."
	    exit 1;;
    esac
done

# test that directory is provided
if [ -z ${dirFastq+x} ]
then
    echo "error...option -d required."
    exit 1
fi

# test that directory contains seq files
if $compress
then
    suffix="dsrc"
else
    suffix="fq.gz"
fi
nfiles=$(find $dirFastq -name "*_1.$suffix" | wc -l)
if [ $nfiles -lt 1 ]
then
    echo "error...empty input directory"
    exit 1
fi

# initialize directories
# remove trailing back slash 
dirData=$(echo $dirFastq | sed -r 's|/$||g')
dirData=$(echo $dirData | sed -r 's|/[^/]+$||g')

# find number of files and determine batches
# each batch consisting of 8 samples (16 files PE or 8 SE)
batches=$((($nfiles - 1)/8 + 1))

# make directories for every batch and move batches of 8 samples into
# their corresponding batch directory
for i in `seq 1 $batches`
do
    mkdir -p $dirData/raw$i
    if $pairEnd
    then
	find $dirFastq -name "*_1.$suffix" | \
	    sed -r "s/_1.${suffix}/_2.${suffix}/g" | \
	    head -8 | xargs -i mv "{}" "$dirData/raw$i"
    fi
    find $dirFastq -name "*_1.$suffix" | head -8 | \
        xargs -i mv "{}" "$dirData/raw$i"
done

# lauch genome indexing
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    genomeGenerate_slurm.sh
cmd="sbatch genomeGenerate_slurm.sh -d $dirData/raw1 -g $genome"

if [[ ! -z $mateLength ]]
then
    cmd="$cmd -m $mateLength"
fi

if $compress
then
    cmd="$cmd -c"
fi

# echo $cmd
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')

# modify preprocessing slurm script
sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
    mRNA.preprocess_slurm.sh
sed -ri "s|^#SBATCH --array=1-.+$|#SBATCH --array=1-${batches}|g" \
    mRNA.preprocess_slurm.sh
sed -ri "s|^#SBATCH --depend=afterok:.+$|#SBATCH --depend=afterok:${slurmid}|g" \
    mRNA.preprocess_slurm.sh

# lauch preprocessing slurm script
cmd="sbatch mRNA.preprocess_slurm.sh -d $dirData -g $genome -a $adapterFile"

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
slurmid=$(eval $cmd | sed -r 's|Submitted batch job ([0-9]*)|\1|g')

if $homolog
then
    # modify homology slurm script
    sed -ri "s|^#SBATCH --mail-user=.+$|#SBATCH --mail-user=${email}|g" \
	mRNA.homolog_slurm.sh
    sed -ri "s|^#SBATCH --array=1-.+$|#SBATCH --array=1-${batches}|g" \
	mRNA.homolog_slurm.sh
    sed -ri "s|^#SBATCH --depend=afterok:.+$|#SBATCH --depend=afterok:${slurmid}|g" \
	mRNA.homolog_slurm.sh
    # lauch preprocessing slurm script
    cmd="sbatch mRNA.homolog_slurm.sh -d $dirData -g $genome"
    eval $cmd
fi
