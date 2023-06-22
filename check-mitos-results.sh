#!/bin/bash

path=$1

cd $path

control=$(grep $sequence *.fas | grep $gene)

function check-sequence-for-gene {
	control=$(grep $sequence *.fas | grep $gene | wc -l)
	if [[ "$control" == 0 ]]; then
		echo "0" | tr "\n" "\t"
	else
		echo "1" | tr "\n" "\t"
	fi
}

echo "		atp6	atp8	cob	cox1	cox2	cox3	nad1	nad2	nad3	nad4	nad4l	nad5	nad6	rrnL	rrnS	trnA	trnC	trnD	trnE	trnF	trnG	trnH	trnI	trnK	trnL1	trnL2	trnM	trnN	trnP	trnQ	trnR	trnS1	trnS2	trnT	trnV	trnW	trnY	trnX"

for sequence in $(grep SRR *.fas | cut -d ':' -f 2 | cut -d '-' -f 2 | uniq | tr "\n" " "); do
	echo $sequence | tr "\n" "\t"
	for gene in atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6 rrnL rrnS trnA trnC trnD trnE trnF trnG trnH trnI trnK trnL1 trnL2 trnM trnN trnP trnQ trnR trnS1 trnS2 trnT trnV trnW trnY trnX; do
		check-sequence-for-gene
	done
	echo
done

exit

# usage: ./check-mitos-results.sh /mnt/data/mitos-results/circular-contigs/zIXn1rTc > exp.txt
