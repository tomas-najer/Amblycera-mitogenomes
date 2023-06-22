#!/bin/bash

sequence=$1 # SRR number
subsampling=$2 # for each fq file, e.g. 2000000 or 4000000

echo $sequence $subsampling

cd /mnt/data/fastq-trimmed
cp $sequence*fastq.gz ~/MitoFinder_container
cd ~/MitoFinder_container
gunzip *gz

shortsubsample=$(echo $subsampling | sed 's/000000//g')

seqtk sample -s100 $sequence*1-trim.fastq $subsampling > $shortsubsample\mil-$sequence-1.fq
seqtk sample -s100 $sequence*2-trim.fastq $subsampling > $shortsubsample\mil-$sequence-2.fq
rm *fastq

./MitoFinder_v1.4 \
	-j $shortsubsample\mil-$sequence \
	-1 $shortsubsample\mil-$sequence-1.fq \
	-2 $shortsubsample\mil-$sequence-2.fq \
	-r all-publ-lice.gb \
	--metaspades -o 5 -p 16 -m 200

rm *$sequence*fq

cd $shortsubsample\mil-$sequence/$shortsubsample\mil-$sequence*metaspades
mv contigs.fasta $shortsubsample\mil-$sequence-contigs.fasta
cp $shortsubsample\mil-$sequence-contigs.fasta ~/TCSF_IMRA
cd ../..
mv $shortsubsample\mil-$sequence* /mnt/data/mitofinder-assembly-results
cd ~/TCSF_IMRA

./TCSF-2.7.1.bash  \
	-i $shortsubsample\mil-$sequence-contigs.fasta \
	-refB all-publ-lice.fasta \
	-o $shortsubsample\mil-$sequence \
	-c 16

mkdir /mnt/data/tcsf-results/$shortsubsample\mil
mv $shortsubsample\mil-$sequence* /mnt/data/tcsf-results/$shortsubsample\mil

exit

# usage: ./assembly.sh sequence subsampling
