#!/bin/bash

# Usage: ./unpack-sra.sh sra-files

# Copy sra file from MetaCentrum
scp tom_zoo@tilia.metacentrum.cz:/storage/du-cesnet/tape_tape/archive/VO_metacentrum/home/tom_zoo/sra/"$*".sra .

# Unpack it into fastq
fasterq-dump --split-files SRR*

# Remove original sra
rm "$*".sra

# Check fastq quality to storage
fastqc "$*"_[1,2].fastq
mv "$*".html /mnt/data/fastqc
mv "$*".zip /mnt/data/fastqc

# Forward fastq for trimming:
mv "$*"_[1,2].fastq ../../bbmap
cd ../../bbmap

# Trim adapters
./trim-adapters.sh "$*"

# Trim right end
./trim-right-end.sh "$*"
mv "$*"_[1,2].fastq /mnt/data/fastq/
rm "$*"_[1,2]-trim-adapters.fastq

# Pack and transport to MitoFinder
gzip "$*"_[1,2]-trim-adapters-right.fastq
mv "$*"_[1,2]-trim-adapters-right.fastq.gz /home/najer/go/src/github.com/sylabs/singularity/MitoFinder_container
cd /home/najer/go/src/github.com/sylabs/singularity/MitoFinder_container

# Run MitoFinder
./MitoFinder_v1.4 -j ./MitoFinder_v1.4 -j "$*" -1 "$*"_1-trim-adapters-right.fastq.gz -2 "$*"_2-trim-adapters-right.fastq.gz -r Amyrsidea-minuta-MH001227.1.gb -o 5 -p 5 -m 10

exit
