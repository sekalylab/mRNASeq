#!/bin/bash

# load packages
module load gcc/6.3.0
module load openmpi/2.0.1
module load python/3.7.0
module load samtools/1.9
module load R/3.5.3

# read arguments
while getopts d:g: option
do
    case "$option" in
	d) dataDir=$OPTARG;;
	g) genome=$OPTARG;;
    esac
done

# set global variables for the script
bin="/mnt/rstor/SOM_PATH_RXS745U/bin"
genomeDir="/mnt/rstor/SOM_PATH_RXS745U/genome/$genome"
gtfFile="$genomeDir/Annotation/genes.gtf"
bedFile="$genomeDir/Annotation/genes.bed"
homologFile="$genomeDir/Annotation/homolog.gtf"
maxProc=8
stranded=reverse

# parallelization code
function queue {
    queue="$queue $1"
    num=$(($num+1))
}

function regeneratequeue {
    oldReQueue=$queue
    queue=""
    num=0
    for pid in $oldReQueue
    do
        if [ -d /proc/$pid ]
        then
            queue="$queue $pid"
            num=$(($num+1))
        fi
    done
}

function checkqueue {
    oldChQueue=$queue
    for pid in $oldChQueue
    do
        if [ ! -d /proc/$pid ]
        then
            regeneratequeue # at least one pid had finished
            break
        fi
    done
}

# select only primary alignment
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: selecting primary alignment..."
    sampleID=$(find $dataDir -name "*.sorted.bam")
    sampleID=$(echo $sampleID | sed -r 's/.sorted.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	$bin/sambamba-0.7.0/sambamba-0.7.0-linux-static \
	    view \
	    -h \
	    -F "not secondary_alignment" \
	    -t 8 \
	    $sample.sorted.bam | \
	    sed -r 's/NH:i:[0-9]+/NH:i:1/g' | \
	    $bin/sambamba-0.7.0/sambamba-0.7.0-linux-static \
		view \
		-f bam \
		-h \
		-o $sample.primary.bam \
		-S \
		-t 8 \
		/dev/stdin
	# deleting unuse files
	rm $sample.sorted.bam
    done
    echo "done"
fi

# 11. RSeQC to infer if dataset is strand-specific
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: infering strand-specific experiment..."
    file=$(find $dataDir -name "*.primary.bam" | head -n 1)
    # extract from infer_experiment output percentage of reads mapping to
    # reference strand (single-end: "++,--" or pair-end "1++,1--,2+-,2-+")
    percentRefStrand=$(grep -E '1\+\+|\"\+\+' $dataDir/rseqc.infer.out | \
        sed -r 's/.+: //g')
    # extract from infer_experiment output percentage of reads mapping to
    # reverse strand (single-end: "+-,-+" or pair-end "1+-,1-+,2++,2--")
    percentReverseStrand=$(grep -E '1\+\-|\"\+\-' $dataDir/rseqc.infer.out | \
        sed -r 's/.+: //g')
    # if one of the percentage above 80% it will be considered stranded 
    if (( $(echo "$percentRefStrand" | awk '{print ($1 >= 0.8)}') ))
    then
        stranded="yes"
    fi
    if (( $(echo "$percentReverseStrand" | awk '{print ($1 >= 0.8)}') ))
    then
        stranded="reverse"
    fi
    echo "done"
fi

# HTSeq to count
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: estimating gene counts..."
    sampleID=$(find $dataDir -name "*.primary.bam")
    sampleID=$(echo $sampleID | sed -r 's/.primary.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	    $bin/HTSeq-0.11.2/htseq-count \
		--mode=union \
		--stranded=$stranded \
		--idattr=gene_id \
		--format=bam \
		--minaqual=0 \
		--samout=$sample.primary.samout \
		--quiet \
		$sample.primary.bam \
		$gtfFile > $sample.primary_counts_gene
    done
    echo "done"
fi

# 13. Modifying counts output
# Removing the last 5 lines from the count table
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: modifying counts output..."
    sampleID=$(find $dataDir -name "*primary_counts_gene")
    sampleID=$(echo $sampleID | sed -r 's/_counts_gene//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	# tail -n 5 ${sample}_counts_gene > $sample.align.stat
	head -n -5 ${sample}_counts_gene > \
	     ${sample}_genecounts
        # deleting unuse files
	rm ${sample}_counts_gene
    done
    echo "done"
fi

# extract reads that mapped to no feature
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: extracting no features reads..."
    sampleID=$(find $dataDir -name "*.primary.bam")
    sampleID=$(echo $sampleID | sed -r 's/.primary.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	    $bin/sambamba-0.7.0/sambamba-0.7.0-linux-static \
		-H $sample.primary.bam > $sample.nofeature.sam
	    $bin/sambamba-0.7.0/sambamba-0.7.0-linux-static \
		$sample.primary.bam |
		paste - <(cut -sf 2 $sample.primary.samout) |
		grep "__no_feature" |
		cut -sf 1-15 >> $sample.nofeature.sam
	    $bin/sambamba-0.7.0/sambamba-0.7.0-linux-static \
		-h -b $sample.nofeature.sam > $sample.nofeature.bam
	    # delete temporary files
	    rm $sample.nofeature.sam
	    rm $sample.primary.bam
	    rm $sample.primary.samout
    done
    echo "done"
fi

# HTSeq to count
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: estimating homolog gene counts..."
    sampleID=$(find $dataDir -name "*.nofeature.bam")
    sampleID=$(echo $sampleID | sed -r 's/.nofeature.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) ) 
    for sample in ${sampleID[@]}
    do
            $bin/HTSeq-0.11.2/htseq-count \
		--mode=union \
		--stranded=$stranded \
		--idattr=gene_id \
		--format=bam \
		--nonunique=all \
		--samout=$sample.homolog.samout \
		--quiet \
		--minaqual=0 \
		$sample.nofeature.bam \
		$homologFile \
		> $sample.homolog_counts_gene
            # delete unused file
	    rm $sample.nofeature.bam
	    rm $sample.homolog_counts_gene
    done
    echo "done"
fi

# combine primary count and homolog samout
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: creating homolog gene counts..."

    sampleID=$(find $dataDir -name "*.homolog.samout")
    sampleID=$(echo $sampleID | sed -r 's/.homolog.samout//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
            Rscript $bin/homolog2counts.R \
                --primary $sample.primary_genecounts \
                --homolog $sample.homolog.samout
            # delete unused file
            rm $sample.primary_genecounts
            rm $sample.homolog.samout
    done
    echo "done"
fi
