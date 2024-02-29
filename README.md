# HLA_autoimmune
This project encompasses all the code utilized in the paper titled "Enhanced HLA Detection and Immune Functionality Analyses Reveal Insights from Autoimmune Diseases and COVID-19 Cohorts."

## Paper Description
The Human Leukocyte Antigen (HLA) locus is associated with a variety of inflammatory conditions. However, traditional HLA detection relies on genotyping with known limitations, prompting our efforts to refine RNA-sequencing-based HLA detection to better define the role of gene expression and disease. We propose an optimized strategy for accurate HLA gene and allele and protein sequence identification, assessment of HLA gene allele-specific, and protein-sequence-specific expression, HLA diversity, Bayesian-based zygosity inference, and contrastive neural network-based HLA haplotype group comparison. Thorough evaluation via long-read Iso-seq and repeated RNA-seq on A549 demonstrated the method's proficiency. We applied the method to two original patient cohorts. The first cohort is comprised of patients with multiple sclerosis (n=105), neuromyelitis optica (n=53), and control subjects (n=68). The second cohort comprises 39 patients with COVID-19 infection, classified as non-severe = 14, severe-survived = 16, severe-deceased = 9. This last cohort was primarily of Native American ancestry, with longitudinal RNA-sequencing data collected during the COVID-19 infection. We found that HLA expression significantly contributed to patient prognosis variability, with enhanced outcomes linked to HLA diversity, improving the recognition of the T cell receptor repertoire computed through V(D)J recombination and the SARS-CoV-2 spike protein. HLA diversity strongly enhances HLA binding affinity. During acute infection, specific HLA genes are upregulated as a response. Our findings indicate that higher HLA expression and diversity advantages extend to acute and chronic scenarios. While higher HLA gene expression and diversity typically signify a robust immune system, specific HLA protein sequences predispose individuals to multiple sclerosis and neuromyelitis optica.

![Framework](https://github.com/Limin-Jiang/HLAdetector/blob/main/HLA_Figure1.jpg)

(A). HLA general description. Of the four HLA levels, our method can detect HLA up to the resolution of HLA protein sequence. (B). A dendrogram representing all HLA-A sequences, demonstrating the number and hierarchical nature of HLA sequences. (C). A graphical representation of our HLA alignment strategy and detection strategies. 

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
  --help   Display this help message.
```
## Use cases
This is a example for input a bam file:
```bash
./HLA_autoimmune -SwithB -i data/Example.bam -v hg38 -er 0.02 -c 0 -o results/output -ex
```

This is a example for input one or two fq file:
```bash
./HLA_autoimmune -SwithF -i file1.fq,file2.fq -v hg38 -er 0.02 -c 0 -o results/output -ex
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




## Contact

Any Comments or questions, Please contact:

Yan Guo, Ph.D, yanguo1978@gmail.com

Limin Jiang, Ph.D, lxj423@miami.com
