#!/bin/bash

cd /home/tom/Desktop

missinggenes=$(cat missing-genes.txt | cut -d " " -f 2 | sort | uniq | tr "\n" " ")
echo $missinggenes
for gene in $missinggenes; do
	samples=$(cat missing-genes.txt | grep $gene | cut -d " " -f 1 | tr "\n" " ")
	for sample in $samples; do
		contigs=$(cat final-contigs-for-mitos.fa | grep ">" | tr -d ">" | grep $sample | tr "\n" " ")
		for contig in $contigs; do
			cd /home/tom/Desktop/final-annotation
			mkdir $gene
			folders=$(ls | grep $gene | tr "\n" " ")
			cp */$contig*orf.txt $gene
			cp */$contig*-orf-prot.fasta $gene
			cp */$contig*-orf-nucl.fasta $gene
			cp */*-$gene\$-final.fasta $gene
			cp */*-$gene\$-final-AA.fasta $gene
			cd $gene
			cat *-$gene\$-final-AA.fasta > $gene-prot.fasta
			cat $contig*-orf-prot.fasta > $contig-$gene-database-prot.fasta
			#databaseinput=$(ls | grep orf-prot.fasta$)
			#echo database input is $databaseinput
			makeblastdb -in $contig-$gene-database-prot.fasta -dbtype prot	
			blastp -db $contig-$gene-database-prot.fasta -query $gene-prot.fasta > $contig-$gene.blastp
			cd /home/tom/Desktop
		done
	done
done

exit

cp Anvag*orf-prot.fasta ../atp8

cat */*faa | awk '{if (substr($0,1,1)==">"){print "\n"$0} else printf("%s",$0);p++;}END{print "\n"}' | tail -n+2

export PATH=$PATH:/home/tom/ncbi-blast-2.13.0+/bin
makeblastdb -in Zumac-1308-orf-prot.fasta -dbtype prot
blastp -db Zumac-1308-orf-prot.fasta -query atp8-prot.fasta > Zumac-atp8.blastp

export PATH=$PATH:/home/tom/mauve_snapshot_2015-02-13
