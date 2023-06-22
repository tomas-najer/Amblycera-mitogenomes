#!/bin/bash

cd /home/tom/Documents/Amblycera-mt-genomes/csv-final-contigs/non-circular-fragments

for file in Achop-1011.csv  ApspApsp-934.csv  HsspPrgar-587.csv  Latin-928.csv  Mccos-1494.csv  Mccos-2315.csv  Mebal-1135.csv  Zumac-1308.csv  Zumac-2986.csv; do
	fragment=$(echo $file | cut -d "." -f 1)
	louse=$(echo $fragment | cut -d "-" -f 1)
	for row in $(cat $file | grep $louse | tr "\n" " "); do
		start=$(echo $row | cut -d "," -f 2)
		end=$(echo $row | cut -d "," -f 3)
		gene=$(echo $row | cut -d "," -f 4)
		score=$(echo $row | cut -d "," -f 5)
		orientation=$(echo $row | cut -d "," -f 6)
		is_trn=$(echo $gene | grep "trn" | wc -l)
		is_rrn=$(echo $gene | grep "rrn" | wc -l)
		if [[ "$is_trn" == 1 ]]; then
			genetype="tRNA"
		elif [[ "$is_rrn" == 1 ]]; then
			genetype="rRNA"
		else
			genetype="gene"
		fi
		has_ORF=$(echo $genetype | grep "gene" | wc -l)
		if [[ "$has_ORF" == 1 ]]; then
			source="ORF"
		else
			source="mitos"
		fi
		
		echo "$louse	$source	$genetype	$start	$end	$score	$orientation	"."	Name=$gene" | sed "s/\tF\t/\t-\t/g" | sed "s/\tT\t/\t+\t/g" >> $fragment.gff
	done
done

exit
