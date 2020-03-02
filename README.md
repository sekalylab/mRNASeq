## Sequencing pipeline 1 (S1)
Contains scripts for mRNA-Seq preprocessing

<!-- badges: start -->
![GitHub last commit](https://img.shields.io/github/last-commit/sekalylab/mRNAseq)
<!-- badges: end -->

## Dependencies
- dsrc (version 2.0.2) [optional]
- FastQC (version 0.11.8) [optional]
- RSeQC (version 3.0.1)
- Trimmomatic (version 0.39)
- STAR aligner (version 2.7.0e)
- HTSeq (version 0.11.2)

## Usage
#### Add python libraries to your PATH
The different applications and libraries used to run the preprocessing pipelines
have been installed on the CWRU high-performance computing (HPC) cluster, in the
shared directory `/mnt/rstor/SOM_PATH_RXS745U`. To be able to run the python
based applications (HTSeq, RSeQC) you will need to append the path to the shared
directory to the environment variable `$PYTHONPATH`.
```bash
export PYTHONPATH="$PYTHONPATH:/mnt/rstor/SOM_PATH_RXS745U/lib/python3.7/site-packages"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/mnt/rstor/SOM_PATH_RXS745U/lib"
```

#### Rename raw sequencing files
The pipeline expect raw sequencing files to be named using the convention
*SampleID*_*mate*.*suffix* where:  
*SampleID* can be any string  
*mate* is the digit [1,2]  
*suffix* is "fq.gz" for FASTQ files and "dsrc" for DSRC compress files  
ex: S01_1.fq.gz

#### Run mRNA-Seq pipeline
The mRNA-Seq pipeline will map sequencing reads to a reference
genome
```bash
bash mRNA.preprocessing_master.sh -d {raw_directory}

arguments:  
d=[d]irectory with raw data (directory; required)  
g=reference [g]enome  
    accepted values: GRC3h8, Mmul_8, GRCm38, MacFas5 
a=FASTA file with [a]dapters sequences (file)  
    default value: TruSeq3-SE.fa  
m=[m]ate length (integer)  
    if empty will be determine automatically  
c=DSRC [c]ompress files given as input  
    argument -c is not given, FASTQ files (.fq.gz) are expected  
p=[p]air-end sequencing files  
    if arugment -p is not given, single-end sequencing files are  
    expected
e=[e]mail address  
i=[i]soform transcript/exon counts  
o=use h[o]mology to annotate reads  
n=align [n]on-host reads to ncbi nt database 
h=print [h]elp
```
