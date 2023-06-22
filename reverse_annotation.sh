#!/bin/bash

cd /home/tom/Documents/Amblycera-mt-genomes/csv-final-contigs

sample=$1
header=$(cat $sample.csv | head -n 1)
echo $header > standardized_direction/$sample.csv

cox1_start=$(cat $sample.csv | grep cox1 | cut -d "," -f 2)
echo $cox1_start
cox1_end=$(cat $sample.csv | grep cox1 | cut -d "," -f 3)
echo $cox1_end
cox1_orientation=$(cat $sample.csv | grep cox1 | cut -d "," -f 6)
echo $cox1_orientation
if [[ $cox1_orientation == "T" ]]; then
	cox1_new_orientation="F"
	elif [[ $cox1_orientation == "F" ]]; then
	cox1_new_orientation="T"
	else
	cox1_new_orientation="F"
	fi
echo $cox1_new_orientation

echo $sample,$cox1_start,$cox1_end,"cox1",'0',$cox1_new_orientation >> standardized_direction/$sample.csv

for gene in trnC  trnY  atp8  atp6  cox3  nad5  trnI  nad1  trnS1  nad6  rrnL  trnT  trnD  trnR  trnL2  trnL1  trnS2  trnQ  nad3  trnF  trnP  rrnS  trnH  trnK  trnA  trnM  nad2  nad4_1  nad4l  cob  cox2  trnV;do

	echo $gene
	gene_old_start=$(cat $sample.csv | grep $gene | cut -d "," -f 2)
	echo $gene_old_start
	gene_old_end=$(cat $sample.csv | grep $gene | cut -d "," -f 3)
	echo $gene_old_end
	gene_length=$(($gene_old_end - $gene_old_start))
	echo $gene_length
	previous_gene_end=$(cat standardized_direction/$sample.csv | tail -n 1 | cut -d "," -f 3)
	echo $previous_gene_end
	gene_new_start=$(($previous_gene_end + 1))
	echo $gene_new_start
	gene_new_end=$(($gene_new_start + $gene_length))
	echo $gene_new_end
	gene_score=$(cat $sample.csv | grep $gene | cut -d "," -f 5)
	echo $gene_score
	gene_old_orientation=$(cat $sample.csv | grep $gene | cut -d "," -f 6)
	echo $gene_old_orientation
	if [[ $gene_old_orientation == "T" ]]; then
		gene_new_orientation="F"
		elif [[ $gene_old_orientation == "F" ]]; then
		gene_new_orientation="T"
		else
		gene_new_orientation="F"
		fi
	echo $gene_new_orientation
	echo $sample,$gene_new_start,$gene_new_end,$gene,$gene_score,$gene_new_orientation >> standardized_direction/$sample.csv

done

original_geneorder=$(cat $sample.csv | cut -d "," -f 4 | tr "\n" " ")
echo "old: $original_geneorder"
new_geneorder=$(tac $sample.csv | cut -d "," -f 4 | tr "\n" " ")
echo "new: $new_geneorder"

exit
