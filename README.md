# HLAdetector
This software aims to detect the HLA gene based on RNA-seq.

## Description
The Human Leukocyte Antigen (HLA) system plays a pivotal role in immune response regulation, serving as a key component in the identification and defense against foreign substances, and its understanding holds profound implications for various aspects of human health and disease. Nevertheless, the considerable polymorphism and homology exhibited by HLA genes pose a formidable challenge in their accurate detection through RNA-seq methodologies. The intricacies arising from the high degree of genetic diversity demand sophisticated analytical approaches to discern and characterize HLA gene profiles with precision and reliability. In the course of our research, we have developed a comprehensive pipeline. This project aims to detection and comparative analysis of HLA genes, leveraging RNA-seq in conjunction with ISO-seq and a Contrastive Learning Neural Network. 

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
git clone https://github.com/Limin-Jiang/HLAdetector.git

# Navigate to the project directory
cd HLAdetector
export PATH=/your/directory/HLAdetector:$PATH
```



## Usage

```bash
Usage: ./myscript [-SwithB | -SwithF] -i <Parameter1> -v <Parameter2> -er <Parameter3> -c <Parameter4> -o <Parameter5>  -p <Parameter6> [-ex]
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
## Example
This is a example for input a bam file:
```bash
./myscript -SwithB -i data/Example.bam -v hg38 -er 0.02 -c 0 -o results/output -ex
```

This is a example for input one or two fq file:
```bash
./myscript -SwithF -i file1.fq,file2.fq -v hg38 -er 0.02 -c 0 -o results/output -ex
```

You will get three file, including:  output_HLAlist.csv, output_genetype.csv, output.expression.csv

