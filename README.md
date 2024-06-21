# ImmuSeeker
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
- [Additional Information](#Additional)
- [Contact](#contact)

## Installation ImmuSeeker

### Via docker 

```bash
sudo docker pull lxj423/immuseeker:latest
```


#### Example

cd your/directory/ (that directory includes your data folder/file)

if the Example.bam is in the folder'./data', please run the following command. 
docker run -v ./data:/ImmuSeeker_data -it immuseeker -HLA -SwithB -i Example.bam -v hg38 -er 0.02 -c 0 -o output -ex -n 0 -dvr -pt


### Manual installation


####  Prerequisites
Make sure the following programs are in your PATH:
- Samtools (Copyright (C) 2008-2024 Genome Research Ltd.)
```bash
apt -y install samtools
```

- bowtie

```bash
apt -y install bowtie
```

- R  >=v3.6 ( With packages: data.table, tidyverse, ggraph, igraph, docopt, phyloseq, stringr)
```bash
apt -y install r-base
Rscript -e "install.packages(c('data.table', 'dplyr', 'ggraph', 'igraph', 'docopt','stringr'))"
Rscript -e "install.packages('permute', repos='https://cloud.r-project.org/')"
Rscript -e "install.packages('https://github.com/vegandevs/vegan/archive/refs/tags/v2.6-4.tar.gz', repos = NULL, type = 'source')"
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('phyloseq')"
```


####  Install
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
Usage: ImmuSeeker [-HLA | -KIR] [-SwithB | -SwithF] -i <Parameter1> -v <Parameter2> -c <Parameter3> -n <Parameter4>  -p <Parameter5> -p1 <Parameter6> -o <Parameter7>  -er <Parameter8>  [-ex] [-pt] [-dvr] [-dve]
Options:
  -HLA     Invoke the HLA calling process.
  -KIR     Invoke the KIR calling process.
  -SwithB  Input alignment BAM file.
  -SwithF  Input FASTQ file.
  -i       Set the input directory and filename. Example: 'inputfile.bam'.
  -v       Input parameter or file name.
           If using -SwithB, specify the genome version with '-v hg37' or '-v hg38'.
           If using -SwithF, specify FASTQ file(s) with '-i file1,file2' (two files) or '-i file' (one file).
  -c       Specify min number of supported unique reads (default: 5).
  -n       Max mismatches in seed (default: -n 0).
  -p       Set a noninformative flat prior to allow the data to have a strong influence on the posterior distribution. (default: -n (1/3,1/3,1/3)).
  -p1      Set the probability of observing an allele in genotype (default: 1/2).
  -o       Set the output directory and filename. Example: 'outfileName'.
  -ex      Specify whether to output gene expression values(default: false).
  -pt      Specify whether to output phylogenetic tree for HLAs(default: false).
  -dvr     Specify whether to output diversity anlysis result based on the number of unique reads.(default: false).
  -dve     Specify whether to output diversity anlysis result based on gene expression(default: false).
  -er      Specify the sequencing error ratio (default: 0.02).
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


## Additional information

This project also encompasses an analysis of HLA diversity alongside a Contrastive Neural Network-based comparison of HLA haplotype groups. The code implementations and sample datasets are stored within the directory labeled "Example_Analysis_CL_and_diversity.



## Contact

Any Comments or questions, Please contact:

Yan Guo, Ph.D, yanguo1978@gmail.com

Limin Jiang, Ph.D, lxj423@miami.com
