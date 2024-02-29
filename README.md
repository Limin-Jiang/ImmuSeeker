# HLAseeker
This software aims to detect the HLA gene based on RNA-seq.

## Description
The Human Leukocyte Antigen (HLA) locus is associated with a variety of inflammatory conditions. However, traditional HLA detection relies on genotyping with known limitations, prompting our efforts to refine RNA-sequencing-based HLA detection to better define the role of gene expression and disease. We propose an optimized strategy for accurate HLA gene and allele identification, sequence-specific expression, and Bayesian-based zygosity inference. 

![Framework](https://github.com/Limin-Jiang/HLAdetector/blob/main/HLA_Figure1.jpg)


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example](#Example)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

## Installation

###  Dependencies
Make sure the following programs are in your PATH:

Samtools v1.18

bowtie v1.3.1

R  >=v3.6 (With packages: data.table, stringr, docopt)

###  Install
To get started with this project, follow these steps:

```bash
# Clone the repository
git clone https://github.com/Limin-Jiang/HLAseeker.git

# Navigate to the project directory
cd HLAseeker
export PATH=/your/directory/HLAseeker:$PATH
```



## Usage

```bash
Usage: ./HLAseeker [-SwithB | -SwithF] -i <Parameter1> -v <Parameter2> -er <Parameter3> -c <Parameter4> -o <Parameter5>  -p <Parameter6> [-ex]
Options:
  -SwithB  Input alignment BAM file.  
  -SwithF  Input FASTQ file.  
  -i       Set the input directory and filename. Example: 'data/outfile.bam'.  
  -v       Input parameter or file name.  
           If using -SwithB, specify the genome version with '-v hg37' or '-v hg38'.           
           If using -SwithF, specify FASTQ file(s) with '-i file1,file2' (two files) or '-i file' (one file).           
  -o       Set the output directory and filename. Example: 'data/outfile'.  
  -ex      Specify whether to output gene expression values(default: false).  
  -er      Specify the sequencing error ratio (default: 0.02).  
  -c       Specify min number of supported reads (default: 0).  
  -p       Specify the probability of three genotypes with '-p (1/4,1/2,1/4)' (default: (1/3,1/3,1/3)).  
  --help   Display this help message.
```
## Use cases
This is a example for input a bam file:
```bash
./HLAseeker -SwithB -i data/Example.bam -v hg38 -er 0.02 -c 0 -o results/output -ex
```

This is a example for input one or two fq file:
```bash
./HLAseeker -SwithF -i file1.fq,file2.fq -v hg38 -er 0.02 -c 0 -o results/output -ex
```


## Output file descriptions

You will get three file, including:  

1) output_HLAlist.csv

2) output_genetype.csv

3) output.expression.csv


## Contact

Any Comments or questions, Please contact:
Yan Guo, Ph.D, yanguo1978@gmail.com
Limin Jiang, Ph.D, lxj423@miami.com
