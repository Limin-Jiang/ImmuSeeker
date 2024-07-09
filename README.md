# ImmuSeeker
This project encompasses all the code utilized in the paper titled "Deep mining of immune-related gene family signatures through lineage reconstruction ."

## Paper Description
Immune gene family signatures, such as Human Leukocyte Antigen (HLA) and Killer-cell Immunoglobulin-like Receptor (KIR), can be extracted from sequencing data. Previous methods offer limited scope, functionality, and preservation. By examining the phylogenetic structures of HLA and KIR, we developed ImmuSeeker for accurate HLA/KIR gene and allele and protein sequence identification, assessment of HLA/KIR gene allele-specific, and protein-sequence-specific expression, diversity, Bayesian-based zygosity inference, and contrastive neural network-based HLA/KIR haplotype group comparison. 

![Framework](https://github.com/Limin-Jiang/HLA_autoimmune/blob/main/Figure.JPG)


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

### Manual installation

####  Prerequisites

- Samtools (Copyright (C) 2008-2024 Genome Research Ltd.)
```bash
apt -y install samtools
```

- bowtie

```bash
apt -y install bowtie
```

- R  >=v3.6 ( With packages: data.table, dplyr, ggraph, igraph, docopt, phyloseq, stringr)
```bash
apt -y install r-base
Rscript -e "install.packages(c('data.table', 'dplyr', 'ggraph', 'igraph', 'docopt','stringr'), repos='https://cran.r-project.org')"
Rscript -e "install.packages('permute', repos='https://cloud.r-project.org/')"
Rscript -e "install.packages('https://github.com/vegandevs/vegan/archive/refs/tags/v2.6-4.tar.gz', repos = NULL, type = 'source')"
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('phyloseq')"
```

####  Install
To get started with this project, follow these steps:

```bash
# Clone the repository
git clone https://github.com/Limin-Jiang/ImmuSeeker.git

# Navigate to the project directory
cd ImmuSeeker
export PATH=/your/directory/ImmuSeeker:$PATH
```



## Usage

```bash
Usage: ImmuSeeker [-HLA | -KIR] [-SwithB | -SwithF] -i <Parameter1> -v <Parameter2> -c <Parameter3> -n <Parameter4>  -p <Parameter5> -p1 <Parameter6> -o <Parameter7>  -er <Parameter8> -ex <Parameter9> -pt <Parameter10> -dv <Parameter11> -nr <Parameter12>
Options:
	-HLA     Invoke the HLA calling process.
	-KIR     Invoke the KIR calling process.
	-SwithB  Input alignment BAM file.
	-SwithF  Input FASTQ file.
	-i       Set the input directory and filename. Example: 'your/input/directory/inputfile.bam'.
	-v       Input parameter or file name. 
	   If using -SwithB, specify the genome version with '-v hg37' or '-v hg38'. 
	   If using -SwithF, specify FASTQ file(s) with '-i file1,file2' (two files) or '-i file' (one file).
	-c       Specify min number of supported unique reads (default: 5).
	-n       Max mismatches in seed (default: 0).
	-p       Set a noninformative flat prior to allow the data to have a strong influence on the posterior distribution. (default: -p '(1/3,1/3,1/3)').
	-p1      Set the probability of observing an allele in genotype (default: 1/2).
	-o       Set the output directory and filename. Example: 'your/output/directory/outfile'.
	-ex      Specify whether to output gene expression values(default: -ex false).
	-pt      Specify whether to output phylogenetic tree for HLAs (default: false).
	-dv      Specify whether to output diversity analysis results based on the number of unique reads ('dvr') or gene expression ('dve'). (default:  false).
	-nr      Specify the type of reads to account for when calculating gene expression. Choose between HLA reads ('HLAn') or total reads ('Totaln'). Calculating using total reads requires more time. (default: 'HLAn').
	-er      Specify the sequencing error ratio (default: 0.02).
	--help   Display this help message.
```
## Example

#### Example based on docker

cd your/directory/ (that directory includes your data folder/file)

if the Example.bam is in the folder'./data', please run the following command. 
```bash
sudo docker run -v ./data:/ImmuSeeker_data -it immuseeker -HLA -SwithB -i 'Example.bam' -o Example_output -v hg38 -c 0 -n 0  -p '(1/3,1/3,1/3)' -p1 0.5 -er 0.02 -ex false -pt false -dv dvr -nr HLAn

```

#### Example in local
This is an example for input a bam file to detect HLA:
```bash
./ImmuSeeker -HLA -SwithB -i 'Example.bam' -o Example_output -v hg38
```

This is an example for input one or two fq file to detect HLA:
```bash
./ImmuSeeker -HLA -SwithF -i Example.fq -er 0.02 -c 0 -o output -ex true -dv dve -pt true
```

or

```bash
./ImmuSeeker -SwithF -i file1.fq,file2.fq -er 0.02 -c 0 -o output -ex true -dv dve -pt true
```

This is an example for detecting KIR:
```bash
./ImmuSeeker -KIR -SwithB -i 'Example.bam' -o Example_output
./ImmuSeeker -KIR -SwithF -i Example.fq -o Example_output
```


## Output File Descriptions

Upon completion of the process, you will receive multiple files:

### output_list.csv: 

This file provides information on the detected HLA names, the unique reads supporting them, gene names, and associated levels.

### output_genetype.csv: 

This file contains the inferred zygosity for each gene, determined through Bayesian analysis.

### output_expression.csv: 

This file presents details on the detected genes, alleles, and protein expression levels.

### output_diversity_ex.csv/output_diversity_re.csv: 

These files present details on the diversity analysis results based on gene expression or the number of supporting reads.

### output_evolution_graphs.pdf: 

This file presents the phylogenetic tree for detected HLAs.

### output.fq: 

This file includes the reads utilized for HLA detection purposes.

## Additional information

This project also encompasses an analysis of HLA diversity alongside a Contrastive Neural Network-based comparison of HLA haplotype groups. The code implementations and sample datasets are stored within the directory labeled "Example_Analysis_CL_and_diversity."


## Contact

Any Comments or questions, Please contact:

Yan Guo, Ph.D, yanguo1978@gmail.com

Limin Jiang, Ph.D, lxj423@miami.com
