#!/bin/bash
# @author Slim Fourati
# @author Aarthi Talla
# @version 0.6

# load modules
module load intel/17
module load openmpi/2.0.1
module load python2/2.7.13
module load samtools/1.8
module load STAR/2.7.0e
# module load xz/5.2.2 # module for pysam

# read arguments
compress=false
pairEnd=false
isoform=false
homolog=false
nonHost=false
while getopts d:g:m:a:cponi option
do
    case "$option" in
	d) dirData=$OPTARG;;
	g) genome=$OPTARG;;
	a) adapterFile=$OPTARG;;
	m) mateLength=$OPTARG;;
	c) compress=true;;
	p) pairEnd=true;;
	o) homolog=true;;
	i) isoform=true;;
	n) nonHost=true;;
    esac
done

# set global variables for the script
bin="/mnt/projects/SOM_PATH_RXS745U/bin"
seqDependencies="/mnt/projects/SOM_PATH_RXS745U/genome/$genome"
genomeFasta="$seqDependencies/Sequence/genome.fa"
gtfFile="$seqDependencies/Annotation/genes.gtf"
bedFile="$seqDependencies/Annotation/genes.bed"
maxProc=8
stranded="no"

# setting project/sample directories
dirDiagnostic=$(echo $dirData | sed -r 's|/[^/]+$|/diagnostic_plot|g')

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

# 1. Converting DSRC to FASTQ
if $compress
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: converting DSRC to FASTQ..."
    sampleID=$(find $dirData -name "*.dsrc")
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
        # format dirData/read{mateNumber}_index_{SampleName}.dsrc
        outputName=$(echo $sample | \
            sed -r 's/.dsrc/.fq.gz/g')
        $bin/dsrc-2.0/dsrc d -s -t8 $sample | gzip > $outputName
	# delete dsrc file
        rm $sample
    done
    echo "done"
fi


# 2. Determine mate length
if [ -z $mateLength ]
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: determining mate length..."
    file=$(find $dirData -name "*.fq.gz" | head -n 1)
    mateLength=$(zcat $file | \
        head -n 4000 | \
        awk 'NR%2==0 {print length($1)}' | \
        sort -rn | \
        head -n 1)
    # echo $mateLength
    echo "done"
fi
genomeDir="$seqDependencies/ggOverhang$(($mateLength -1))"


# 3. FastQC quality control
# Create output directory for FastQC and lauch the app
flag=false
if $flag
then
    suffixe="_fastqc"
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: generating quality control report..."
    for sample in $(find $dirData -name "*.fq.gz")
    do
        sampleID=$(echo $sample | sed -r 's/(.?).fq.gz$/\1/')
        # create directory format <sampleName>_<readNb>_fastqc
        mkdir -p $sampleID$suffixe
        # lauch FASTQC with option -o <outputFile> 
        #                          -t <numberOfParallelProcess>
        #                          -q quiet
        perl $bin/FastQC-0.11.8/fastqc $sample \
	    -o $sampleID$suffixe -t 8 -q
    done
    echo "done"
fi


# 4. Read Stats
# Counting total reads and saving totalReads.txt in "diagnostic_plot" folder
flag=true
if $flag
then
    if [ ! -d $dirDiagnostic ]
    then
	mkdir -p $dirDiagnostic
    fi
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: calculating total reads..."
    sampleID=$(find $dirData -name "*_1.fq.gz")
    sampleID=$(echo $sampleID | sed -r 's/_1.fq.gz//g')
    sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
    for sample in ${sampleID[@]}
    do
	count=$(zcat ${sample}_1.fq.gz | wc -l | awk '{print $1/4}')
	if [ $? != 0 ]
	then
            echo -ne "error\n  error counting total reads $sample\n"
            exit
	fi
	sample=$(echo $sample | sed -r 's/.+\///g')
	printf $sample'\t'$count'\t'"TotalReads"'\n' >> \
	    $dirDiagnostic/ReadStats.txt
    done
    echo "done"
fi


# 5. Trimming reads with "Trimmomatic"
# Raw reads must be trimmed/cleaned of adapter contamination.
# The list of adapters in fasta format must be provided.
# A minimum base quality of 5 is required to be maintained at the leading and
# trailing of the reads.
# A minimum survival length of 36bp is required
# Trimmomatic is launched with options:
#   PE                               pair-ends
#   threads                        : <numberOfParallelProcess
#   phred33|phred64                  quality score format
#   ILLUMINACLIP:<adaptaterSeq.fa> : <maxMM>:<minPhredPalind>:<minPhredSingle>:
#                                    <>:<>               
#   LEADING                        : <minQuality>
#   TRAILING                       : <minQuality>
#   SLIDINGWINDOW                  : <windowLength>:<minQuality>
#   MINLEN                         : <minLength>
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: trimming reads..."
    if $pairEnd
    then
	sampleID=$(find $dirData -name "*_[1-2].fq.gz")
	sampleID=$(echo $sampleID | sed -r 's/_[1-2].fq.gz//g')
	sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
	for sample in ${sampleID[@]}
	do
            # lauch trimmomatic
            java -Xmx1g \
		-jar $bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
		-threads $maxProc \
		-phred33 \
		-trimlog ${sample}_trim.log \
		${sample}_1.fq.gz \
		${sample}_2.fq.gz \
		${sample}_1.trimmed.fq.gz \
		${sample}_1.unpaired.fq.gz \
		${sample}_2.trimmed.fq.gz \
		${sample}_2.unpaired.fq.gz \
		ILLUMINACLIP:$adapterFile:2:30:10:8:TRUE \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36 &> ${sample}_trim.err
            if [ $? != 0 ]
            then
		echo -ne "error\n  unable to trim FastQs in directory $dirData"
		exit
            fi
            # deleting unused files
	    rm ${sample}_1.fq.gz
            rm ${sample}_2.fq.gz
            rm ${sample}_1.unpaired.fq.gz
            rm ${sample}_2.unpaired.fq.gz
            rm ${sample}_trim.log
	    rm ${sample}_trim.err
	done
    else
	sampleID=$(find $dirData -name "*_1.fq.gz")
	sampleID=$(echo $sampleID | sed -r 's/_1.fq.gz//g')
	sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
	for sample in ${sampleID[@]}
	do
            # lauch trimmomatic
            java -Xmx1g \
		-jar $bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
		-threads $maxProc \
		-phred33 \
		-trimlog ${sample}_trim.log \
		${sample}_1.fq.gz \
		${sample}_1.trimmed.fq.gz \
		ILLUMINACLIP:$adapterFile:2:36:10 \
		LEADING:10 \
		TRAILING:10 \
		MAXINFO:50:0.97 \
		MINLEN:17 \
		CROP:27 &> ${sample}_trim.err
            if [ $? != 0 ]
            then
		echo -ne "error\n  unable to trim FastQs in directory $dirData"
		exit
            fi
            # deleting unuse files
	    rm ${sample}_1.fq.gz;
            rm ${sample}_trim.log;
	done	
    fi
    echo "done"
fi


# 6. Loading genome in shared memory of the compute node
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: loading genome in memory..."
    STAR --genomeDir $genomeDir \
         --genomeLoad LoadAndExit &>/dev/null                                    
    if [ $? != 0 ]
    then
        echo -ne "error\n  unable to load genome"
        exit
    fi
    currentDir=$(pwd)
    if [ -e "$currentDir/Aligned.out.sam" ] 
    then
	rm $currentDir/Aligned.out.sam
	rm $currentDir/Log.out
	rm $currentDir/Log.progress.out
	rmdir $currentDir/_STARtmp
    fi
    echo "done"
fi


# 7. Alignment with 'STAR'                                                      
# Alignment of the reads to the reference genome are perfermed on the         
# 'trimmed.fastq' files.                                                      
# The genome is loaded into memory for alignment.                             
# After the end of alignment, the genome is removed from shared memory.       
# Lauch STAR with options:                                                    
#   genomeDir         :  <path/to/dir/where genome has been generated>        
#   genomeLoad        : <mode of shared memory usage for the genome files>    
#   readFilesIn       : <mate_1.trimmed.fastq> <mate_2.trimmed.fastq>         
#   runThreadN        : <numberOfParallelProcess>                             
#   outFileNamePrefix : <prefixOutputFiles>                                   
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: alignment of reads..."
    if $pairEnd
    then
	sampleID=$(find $dirData -name "*_[1-2].trimmed.fq.gz")
	sampleID=$(echo $sampleID | sed -r 's/_[1-2].trimmed.fq.gz//g')
	sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
	for sample in ${sampleID[@]}
	do
            STAR --genomeDir $genomeDir \
		 --genomeLoad LoadAndKeep \
		 --readFilesIn ${sample}_1.trimmed.fq.gz \
		 ${sample}_2.trimmed.fq.gz \
		 --readFilesCommand zcat \
	  	 --runThreadN $maxProc \
 		 --outSAMtype BAM Unsorted \
		 --outFileNamePrefix ${sample}_star \
		 --outReadsUnmapped Fastx &>/dev/null
            if [ $? != 0 ]
            then
		echo -ne "error\n  unable to aligned read in directory $sample"
		exit 1
            fi
            # delete trimmed FASTQ
            rm ${sample}_1.trimmed.fq.gz
	    rm ${sample}_2.trimmed.fq.gz
	    rm ${sample}_starLog.out
	    rm ${sample}_starLog.progress.out
	    if ! $isoform
	    then
		rm ${sample}_starSJ.out.tab
	    fi
	    if ! $nonHost
	    then
		rm ${sample}_starUnmapped.out.mate1
		rm ${sample}_starUnmapped.out.mate2
	    fi
	done
    else
	sampleID=$(find $dirData -name "*_1.trimmed.fq.gz")
	sampleID=$(echo $sampleID | sed -r 's/_1.trimmed.fq.gz//g')
	sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
	for sample in ${sampleID[@]}
	do
            STAR --genomeDir $genomeDir \
		 --genomeLoad LoadAndKeep \
		 --readFilesIn ${sample}_1.trimmed.fq.gz \
		 --readFilesCommand zcat \
	         --runThreadN 8 \
       		 --outSAMtype BAM Unsorted \
		 --outFileNamePrefix ${sample}_star \
	         --outReadsUnmapped Fastx \
       		 --outFilterMismatchNoverLmax 0.05 \
		 --outFilterMatchNmin 16 \
		 --outFilterScoreMinOverLread 0 \
		 --outFilterMatchNminOverLread 0 \
		 --alignIntronMax 1 &>/dev/null                    
            if [ $? != 0 ]
            then
		echo -ne "error\n  unable to aligned read in directory $sample"
		exit
            fi
            # delete trimmed FASTQ
            rm ${sample}_1.trimmed.fq.gz
	    # rm ${sample}_starAligned.out.bam
	    rm ${sample}_starLog.out
	    rm ${sample}_starLog.progress.out	
	done
    fi
    echo "done"
fi


# 8. Remove genome from memory
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: removing genomes from memory..."
    STAR --genomeDir $genomeDir \
         --genomeLoad Remove &>/dev/null
    if [ $? != 0 ]
    then
        echo -ne "error\n  unable to remove genome"
        exit
    fi
    # delete STAR temporary files
    currentDir=$(pwd)
    if [ -e $currentDir/Aligned.out.sam ]
    then
	rm $currentDir/Aligned.out.sam
	rm $currentDir/Log.progress.out
	rm $currentDir/Log.out
	rmdir $currentDir/_STARtmp
    fi
    echo "done"
fi


# 9. Trimmed and Mapped read stats
# Counting reads surviving after trimming and reads aligned and writing to
# 'Readstats.txt' in "diagnostic plot" folder.
flag=true
if $flag
then
     currentDate=$(date +"%Y-%m-%d %X")
     echo -ne "$currentDate: calculating surviving and mapped reads..."
     sampleID=$(find $dirData -name "*_starLog.final.out")
     sampleID=$(echo $sampleID | sed -r 's/_starLog.final.out//g')
     sampleID=( $(echo $sampleID | tr ' ' '\n' | sort | uniq) )
     for sample in ${sampleID[@]}
     do
	 surviving=$(sed -n '6p' ${sample}_starLog.final.out | \
	     sed -r 's/.*\|\t(.*)/\1/')
	 Uniqmapped=$(sed -n '9p' ${sample}_starLog.final.out | \
	     sed -r 's/.*\|\t(.*)/\1/')
	 multiMap=$(sed -n '24p' ${sample}_starLog.final.out | \
             sed -r 's/.*\|\t(.*)/\1/')
	 # delete unused files
	 rm ${sample}_starLog.final.out
	 
	 sample=$(echo $sample | sed -r 's/.+\///g')
	 printf $sample'\t'$surviving'\t'"Surviving"'\n' >> \
             $dirDiagnostic/ReadStats.txt
	 printf $sample'\t'$Uniqmapped'\t'"UniqMapped"'\n' >> \
	     $dirDiagnostic/ReadStats.txt
	 printf $sample'\t'$multiMap'\t'"Multimapped"'\n' >> \
	     $dirDiagnostic/ReadStats.txt
     done
     echo "done"
fi


# 10. Order BAM file by name
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: sorting bam file by name..."
    for sample in $(find $dirData -name "*_starAligned.out.bam")
    do
	sampleID=$(echo $sample | sed -r 's/(.?)_starAligned.out.bam$/\1/')
	samtools sort \
		 -n \
		 -@ $maxProc \
		 -o $sampleID.sorted.bam \
		 $sample &>/dev/null
        # delete unsorted bam
        rm $sample
    done
    echo "done"
fi


# 11. RSeQC to infer if dataset is strand-specific
flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: infering strand-specific experiment..."
    file=$(find $dirData -name "*.sorted.bam" | head -n 1)
    $bin/RSeQC-2.6.4/infer_experiment.py \
	-r $bedFile \
	-i $file &> $dirData/rseqc.infer.out
    # extract from infer_experiment output percentage of reads mapping to
    # reference strand (single-end: "++,--" or pair-end "1++,1--,2+-,2-+")
    percentRefStrand=$(grep -E '1\+\+|\"\+\+' $dirData/rseqc.infer.out | \
        sed -r 's/.+: //g')
    # extract from infer_experiment output percentage of reads mapping to
    # reverse strand (single-end: "+-,-+" or pair-end "1+-,1-+,2++,2--")
    percentReverseStrand=$(grep -E '1\+\-|\"\+\-' $dirData/rseqc.infer.out | \
        sed -r 's/.+: //g')
    # if one of the percentage above 80% it will be considered stranded 
    if (( $(echo "$percentRefStrand" | awk '{print ($1 >= 0.75)}') ))
    then
        stranded="yes"
    fi
    if (( $(echo "$percentReverseStrand" | awk '{print ($1 >= 0.75)}') ))
    then
        stranded="reverse"
    fi
    echo "done"
fi


# 12. Estimating transcript abundance with "HTSeq"
# Given a file with aligned reads and the list of genomic features, the number
# of reads mapped to each feature is counted, which accounts for the transcript
# abundance estimation. For each position 'i' in the read, a set S(i) is the
# set of all features overlapping position 'i'. The set 'S' which is the
# intersection of all non empty sets 'S(i)' is considered and the counts are
# generated.
# Launch HTSeq with options:
#   m        : <modeOfCounting>
#   stranded : <data is from strand-specific assay or no>
#   i        : <GTF attribute to be used as feature to count>
if ! $homolog
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: estimating transcript abundance..."
    for sample in $(find $dirData -name "*.sorted.bam")
    do
        (
            sampleID=$(echo $sample | sed -r 's/(.?).sorted.bam$/\1/')
            $bin/HTSeq-0.9.1/htseq-count \
                --mode=union \
                --stranded=$stranded \
                --idattr=gene_id \
		--format=bam \
                --quiet \
                $sample \
                $gtfFile \
                > ${sampleID}_counts_gene
	    if $isoform
	    then
		$bin/HTSeq-0.9.1/htseq-count \
                    --mode=union \
                    --stranded=$stranded \
                    --idattr=transcript_id \
		    --format=bam \
                    --quiet \
                    $sample \
                    $gtfFile \
                    > ${sampleID}_counts_transcript
		python $bin/HTSeq-0.9.1/htseq-count \
                    --mode=union \
                    --stranded=$stranded \
                    --idattr=exon_id \
		    --format=bam \
                    --quiet \
                    $sample \
                    $gtfFile \
                    > ${sampleID}_counts_exon
	    fi
	    rm $sample
	)&
        pid=$!
        queue $pid
        while [ $num -ge $maxProc ]
        do
            checkqueue
            sleep 0.5
        done
    done
    while [ $num -gt 0 ]
    do
        checkqueue
        sleep 0.5
    done
    echo "done"
fi


# 13. Modifying counts output
# Removing the last 5 lines from the count table
if ! $homolog
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: modifying counts output..."
    for sample in $(find $dirData -name "*_counts_gene")
    do
        (
            sampleID=$(echo $sample | sed -r 's/(.?)_counts_gene$/\1/')
	    tail -n 5 ${sampleID}_counts_gene > ${sampleID}.align.stat
	    head -n -5 ${sampleID}_counts_gene > \
		${sampleID}_genecounts
            # deleting unuse files
	    rm ${sampleID}_counts_gene
	    if $isoform
	    then
		head -n -5 ${sampleID}_counts_transcript > \
		    ${sampleID}_transcriptcounts
		head -n -5 ${sampleID}_counts_exon > \
		    ${sampleID}_exoncounts
		rm ${sampleID}_counts_transcript
		rm ${sampleID}_counts_exon
	    fi
	)
    done
    echo "done"
fi
