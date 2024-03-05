# HLA_autoimmune
This project encompasses all the code utilized in the paper titled "Enhanced HLA Detection and Immune Functionality Analyses Reveal Insights from Autoimmune Diseases and COVID-19 Cohorts."

## Paper Description
The Human Leukocyte Antigen (HLA) locus is associated with a variety of inflammatory conditions. However, traditional HLA detection relies on genotyping with known limitations, prompting our efforts to refine RNA-sequencing-based HLA detection to better define the role of gene expression and disease. We propose an optimized strategy for accurate HLA gene and allele and protein sequence identification, assessment of HLA gene allele-specific, and protein-sequence-specific expression, HLA diversity, Bayesian-based zygosity inference, and contrastive neural network-based HLA haplotype group comparison. 

![Framework](https://github.com/Limin-Jiang/HLA_autoimmune/blob/main/Figure.JPG)

HLA general description. Of the four HLA levels, our method can detect HLA up to the resolution of HLA protein sequence. 
A dendrogram representing all HLA-A sequences, demonstrating the number and hierarchical nature of HLA sequences. 
A graphical representation of our HLA alignment strategy and detection strategies. 

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example](#Example)
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
git clone https://github.com/Limin-Jiang/HLA_autoimmune.git

# Navigate to the project directory
cd HLA_autoimmune
export PATH=/your/directory/HLA_autoimmune:$PATH
```



## Usage

```bash
Usage: ./HLA_autoimmune [-SwithB | -SwithF] -i <Parameter1> -v <Parameter2> -er <Parameter3> -c <Parameter4> -o <Parameter5>  -p <Parameter6> [-ex]
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
  -n       Max mismatches in seed (default: -n 0). 
  --help   Display this help message.
```
## Example
This is a example for input a bam file:
```bash
./HLA_autoimmune -SwithB -i data/Example.bam -v hg38 -er 0.02 -c 0 -o results/output -ex -n 0
```

This is a example for input one or two fq file:
```bash
./HLA_autoimmune -SwithF -i file.fq -v hg38 -er 0.02 -c 0 -o results/output -ex -n 0
```

or

```bash
./HLA_autoimmune -SwithF -i file1.fq,file2.fq -v hg38 -er 0.02 -c 0 -o results/output -ex -n 0
```


## Output File Descriptions

Upon completion of the process, you will receive three files:

### output_HLAlist.csv: 

This file provides information on the detected HLA names, the unique reads supported, gene names, and associated levels.

### output_genetype.csv: 

This file contains the inferred zygosity for each gene, determined through Bayesian analysis.

### output_expression.csv: 

This file presents details on the detected genes, alleles, and protein expression levels.

### output.fq: 

This file includes the reads utilized for HLA detection purposes.




## contact

Any Comments or questions, Please contact:

Yan Guo, Ph.D, yanguo1978@gmail.com

Limin Jiang, Ph.D, lxj423@miami.com
