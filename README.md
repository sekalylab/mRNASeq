## Sequencing pipeline
Contains scripts for mRNA-Seq preprocessing

<!-- badges: start -->
![GitHub last commit](https://img.shields.io/github/last-commit/sekalylab/mRNAseq/aws)
<!-- badges: end -->

## Dependencies
- RSeQC (version 4.0.0)  
- STAR aligner (version 2.7.10a)  
- HTSeq (version 2.0.2)  

## Usage
#### Rename raw sequencing files
The pipeline expect raw sequencing files to be named using the convention
*SampleID*_*mate*.*suffix* where:  
*SampleID* can be any string  
*mate* is the digit [1,2]  
*suffix* is "fq.gz" for FASTQ files
ex: S01_1.fq.gz

#### Run mRNA-Seq pipeline
The mRNA-Seq pipeline will map sequencing reads to a reference
genome
```bash
bash mRNA.preprocessing_master.sh -d {raw_directory}

arguments:  
d=[d]irectory with raw data (directory; required)  
g=reference [g]enome  
    accepted values: GRCh38, Mmul_10, Mnem_1  
p=[p]air-end sequencing files  
    if arugment -p is not given, single-end sequencing files are  
    expected
e=[e]mail address 
i=[i]soform transcript/exon counts  
h=print [h]elp
```
