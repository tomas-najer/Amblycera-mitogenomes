#!/bin/bash

path=$1

cd $path
rm -r nonaligned aligned
mkdir nonaligned aligned

cd ~/scripts

function get-gene-from-mitos {
	./get-gene-from-mitos.sh $sequence $gene $path
}
	

for gene in atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6 rrnL rrnS trnA trnC trnD trnE trnF trnG trnH trnI trnK trnL1 trnL2 trnM trnN trnP trnQ trnR trnS1 trnS2 trnT trnV trnW trnY trnX; do
	echo $gene
	for sequence in $(grep SRR $path/*.fas | cut -d ':' -f 2 | cut -d '-' -f 2 | sort -u | tr "\n" " "); do
		orientation=$(./get-gene-from-mitos.sh $sequence $gene $path | grep ">" | sed 's/ -/0/g' | cut -d ';' -f 3)
		echo $sequence
		if [[ $orientation == 0 ]]; then
			get-gene-from-mitos | grep ">" | sed 's/ -/ +/g' >> $path/nonaligned/$gene-nonalign.fa
			get-gene-from-mitos | grep -A 1 ">" | tail -n+2 | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev >> $path/nonaligned/$gene-nonalign.fa
		else	
			get-gene-from-mitos | grep -A 1 ">" >> $path/nonaligned/$gene-nonalign.fa
		fi
	done
	mafft --auto $path/nonaligned/$gene-nonalign.fa > $path/aligned/$gene-align.fa
done

exit

# reverse-complement script
tr -d "\n " < input.txt | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev

#usage: ./make-genes-alignments.sh path; e.g. ./make-genes-alignments.sh /mnt/data/mitos-results/circular-contigs/zIXn1rTc
