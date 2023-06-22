#!/bin/bash


sequence=$1 # SRR number
contig=$2 # number of contig
subsampling=$3 # 2mil or 10mil

echo $sequence
echo $subsampling

cd /mnt/data/tcsf-results/TCSF-final-results/$subsampling
cp $subsampling-$sequence-$contig-tcsf.fasta ~/awa-master

cd /mnt/data/fastq-trimmed
cp $sequence*-trim.fastq.gz ~/awa-master

cd ~/awa-master
gunzip *.gz
# mv $sequence-$contig-tcsf.fasta $subsampling-$sequence-$contig-tcsf.fasta

function part1 {
	for k in 0 5 10 15 20 25 30 35 40; do
		python3 awa-trim.py -s $subsampling-$sequence-$contig-tcsf.fasta -w ${k} -l 0 -o $subsampling-$sequence-$contig-${k}-mer.fa > $subsampling-$sequence-$contig-${k}-trim.log
	done
}

function part2 {
	./interleave_fastq.sh $sequence*1-trim.fastq $sequence*2-trim.fastq > $sequence.fq
	rm *fastq
}

function part3 {
	for k in 0 5 10 15 20 25 30 35 40; do
		python3 awa-map.py -i $subsampling-$sequence-$contig-${k}-mer -f $subsampling-$sequence-$contig-${k}-mer.fa -r $sequence.fq -s $subsampling-$sequence-$contig-${k}-mer.sam -p 16 -v &> $subsampling-$sequence-$contig-${k}-mer.log
	done
}

function part4 {
	for k in 0 5 10 15 20 25 30 35 40; do
		echo "$subsampling-$sequence-$contig-${k}-mer.log" >> $subsampling-$sequence-$contig-final.map.log
		tail -n 7 $subsampling-$sequence-$contig-${k}-mer.log >> $subsampling-$sequence-$contig-final.map.log
	done
}

part2 # Prepare reads for awa-map
part1 # Clip putative circular sequences from within the given scaffold
part3 # Get stats for each putative fragment
part4 # Concatenate the stats into one file

cd /mnt/data/awa-results
mkdir $sequence

cd ~/awa-master
mv *log /mnt/data/awa-results/$sequence
mv *mer.fa /mnt/data/awa-results/$sequence
mv *mer.sam /mnt/data/awa-results/$sequence
rm *$sequence* *bt2* *fq

exit
