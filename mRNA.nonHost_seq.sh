#!/bin/bash

# load packages
module load gcc/6.3.0
module load samtools/1.9

# read arguments
while getopts d: option
do
    case "$option" in
	d) dataDir=$OPTARG;;
    esac
done

# set global variables for the script
bin="/mnt/projects/SOM_PATH_RXS745U/bin"
maxProc=8

# extract host mismatch reads
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: extract mismatch reads..."
    sampleID=$(find $dataDir -name "*_starAligned.out.bam")
    sampleID=$(echo $sampleID | sed -r 's/_starAligned.out.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	samtools view \
		 -h \
		 ${sampleID}_starAligned.out.bam | \
	    grep -v "nM:i:0" | \
	    samtools view -b > ${sampleID}_mm.bam
	samtools fastq \
		 -1 ${sampleID}_mm_1.fq \
		 -2 ${sampleID}_mm_2.fq \
		 -s ${sampleID}_singleton.fq \
		 ${sampleID}_mm.bam
	# deleting unused files
	rm ${sampleID}_mm.bam
	# rm ${sampleID}_starAligned.out.bam 
	rm ${sampleID}_singleton.fq
    done
    echo "done"
fi

# concatenate mismatch and unaligned reads
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: concatenating mismatch and unaligned reads..."
    sampleID=$(find $dataDir -name "*_starAligned.out.bam")
    sampleID=$(echo $sampleID | sed -r 's/_starAligned.out.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	cp ${sampleID}_mm_1.fq ${sampleID}_mm_um_1.fq
	cat ${sampleID}_starUnmapped.out.mate1 >> ${sampleID}_mm_um_1.fq
	cp ${sampleID}_mm_2.fq ${sampleID}_mm_um_2.fq
	cat ${sampleID}_starUnmapped.out.mate2 >> ${sampleID}_mm_um_2.fq
	# deleting unused files
	rm ${sampleID}_mm_1.fq
	rm ${sampleID}_mm_2.fq
	rm ${sampleID}_starUnmapped.out.mate1
	rm ${sampleID}_starUnmapped.out.mate2
    done    
    echo "done"
fi

# align reads to nt database
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: aligning reads to nt db..."
    sampleID=$(find $dataDir -name "*_starAligned.out.bam")
    sampleID=$(echo $sampleID | sed -r 's/_starAligned.out.bam//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do    
	$bin/ncbi-magicblast-1.4.0/bin/magicblast \
	    -query ${sample}_mm_um_1.fq \
	    -query_mate ${sample}_mm_um_2.fq \
	    -db nt \
	    -infmt fastq \
	    -num_threads 24 \
	    -splice false \
	    -perc_identity 100 \
	    -reftype transcriptome \
	    -no_unaligned \
	    -no_discordant \
	    -max_db_word_count 10 \
	    -limit_lookup true \
	    -out ${sample}_magicblast.sam
	# delete unused files
	rm ${sample}_mm_um_1.fq
	rm ${sample}_mm_um_2.fq
    done
    echo "done"
fi
