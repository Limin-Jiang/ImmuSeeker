# Use an official Python runtime as a parent image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get -y install libssl-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install libfontconfig1-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install libcairo2-dev 
RUN apt-get -y install libharfbuzz-dev libfribidi-dev
RUN apt-get -y install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev


RUN apt -y install samtools
RUN apt -y install bowtie
RUN apt -y install r-base

RUN R -e "install.packages(c('data.table', 'dplyr', 'ggraph', 'igraph', 'docopt','stringr'))"
RUN R -e "install.packages('permute', repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('https://github.com/vegandevs/vegan/archive/refs/tags/v2.6-4.tar.gz', repos = NULL, type = 'source')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('phyloseq')"

# Set the working directory
WORKDIR /app

RUN mkdir /ImmuSeeker_data

# Copy the current directory contents into the container at /app
COPY . /app

# Define default command to run when the container starts
#CMD ["bash"]
COPY ImmuSeeker /usr/local/bin/ImmuSeeker
RUN chmod +x /usr/local/bin/ImmuSeeker

#CMD ["R"]

ENTRYPOINT ["ImmuSeeker"]