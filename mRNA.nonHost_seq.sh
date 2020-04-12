#!/bin/bash

# read arguments
pairEnd=false
while getopts d:p option
do
    case "$option" in
	d) dirData=$OPTARG;;
	p) pairEnd=true;;
    esac
done

flag=true
if $flag
then
    currentDate=$(date +"%Y-%m-%d %X")
    echo -ne "$currentDate: aligning with kraken2..."
    sampleID=$(find $dirData -name "*.out.mate1")
    sampleID=$(echo $sampleID | sed -r 's/.out.mate1//g')
    if $pairEnd
    then
	/mnt/rstor/SOM_PATH_RXS745U/bin/kraken2_2.0.8/kraken2 \
	    --db /mnt/rstor/SOM_PATH_RXS745U/genome/kraken \
	    --paired \
	    ${sampleID}.out.mate1 \
	    ${sampleID}.out.mate2 \
	    --output $sampleID.kraken2.out \
	    --report $sampleID.kraken2.report \
	    --report-zero-counts
	# clean up
	rm ${sampleID}.out.mate1
	rm ${sampleID}.out.mate2
    else
	/mnt/rstor/SOM_PATH_RXS745U/bin/kraken2_2.0.8/kraken2 \
            --db /mnt/rstor/SOM_PATH_RXS745U/genome/kraken \
            ${sampleID}.out.mate1 \
            --output $sampleID.kraken2.out \
            --report $sampleID.kraken2.report \
            --report-zero-counts
        # clean up
        rm ${sampleID}.out.mate1
    fi
    echo "done"
fi
