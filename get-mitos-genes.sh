
#######

# I need to prepare separate genes, produced by mitos (from TCSF) and aTRAM, so I can align them in mafft + jalview

#!/bin/bash

#######

# each gene will have 4 versions - 2m-mitos, 2m-atram, 10m-mitos, 10m-atram

# we have 10m in sequences: 

# SRR8177109
# SRR8177108
# SRR8175168
# SRR8145996
# SRR8146021
# SRR8172494
# SRR8146022
# SRR8145949
# SRR8145959

# then, we have only 2m:

# SRR8582564 - 1 nice contig from tcsf
# SRR8334234 - 1 nice contig from tcsf
# SRR8334280 - 1 nice contig from tcsf
# SRR8175167 - 3 nice contigs from tcsf
# SRR8145981 - 3 nice contigs from tcsf
# SRR8175169 - 3 nice contigs from tcsf

# from tcsf -> mitos, we will use results from circularized contigs! (after processing by awa-trim.py)

#######

# first, make directories for all the genes

cd /mnt/data/gene-alignments
mkdir Pro Thr cox1 Cys Met Trp Glu Arg His cox3 Ser nad1 cob Leu nad4 nad4L nad2 Gly nad5 Lys Asn Ala Tyr cox2 nad3 Asp Ile Phe 16S 12S atp6 atp8 Gln Val nad6 Leu2 Ser2

# then extract the genes from mitos results, first, create separate faste file for group of genes belonging to one sequence

grep -A 1 SRR8177109 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8177109-all-genes.fa
grep -A 1 SRR8177108 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8177108-all-genes.fa
grep -A 1 SRR8175168 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8175168-all-genes.fa
grep -A 1 SRR8145996 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8145996-all-genes.fa
grep -A 1 SRR8146021 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8146021-all-genes.fa
grep -A 1 SRR8172494 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8172494-all-genes.fa
grep -A 1 SRR8146022 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8146022-all-genes.fa
grep -A 1 SRR8145949 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8145949-all-genes.fa
grep -A 1 SRR8145959 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8145959-all-genes.fa
grep -A 1 SRR8582564 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8582564-all-genes.fa
grep -A 1 SRR8334234 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8334234-all-genes.fa
grep -A 1 SRR8334280 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8334280-all-genes.fa
grep -A 1 SRR8175167 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8175167-all-genes.fa
grep -A 1 SRR8145981 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8145981-all-genes.fa
grep -A 1 SRR8175169 /mnt/data/mitos-results/circular-contigs/zpuCTLKc/*.fas > /mnt/data/gene-alignments/SRR8175169-all-genes.fa

# then, for each combination of sequence and gene I need to get this:

# grep -A 1 nad4 /mnt/data/gene-alignments/SRR8177109-all-genes.fa > /mnt/data/gene-alignments/nad4/SRR8177109-nad4.fa

# this is an immense number of grep scripts, so I cannot write it manually. Let's write set of scripts which can do it. Firstly, for one sequence:

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8177109-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8177109-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8177109-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-1st-paste.txt
sed 's/$/BSRR8177109-/g' /mnt/data/gene-alignments/grep-scripts-SRR8177109-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8177109-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8177109-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8177109-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8177109-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177109-final.txt

# then repeat the same process for all the others:

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8177108-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8177108-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8177108-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-1st-paste.txt
sed 's/$/BSRR8177108-/g' /mnt/data/gene-alignments/grep-scripts-SRR8177108-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8177108-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8177108-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8177108-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8177108-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8177108-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8175168-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8175168-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175168-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-1st-paste.txt
sed 's/$/BSRR8175168-/g' /mnt/data/gene-alignments/grep-scripts-SRR8175168-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175168-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8175168-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8175168-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8175168-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175168-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8145996-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8145996-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145996-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-1st-paste.txt
sed 's/$/BSRR8145996-/g' /mnt/data/gene-alignments/grep-scripts-SRR8145996-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145996-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8145996-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8145996-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8145996-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145996-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8146021-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8146021-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8146021-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-1st-paste.txt
sed 's/$/BSRR8146021-/g' /mnt/data/gene-alignments/grep-scripts-SRR8146021-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8146021-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8146021-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8146021-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8146021-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146021-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8172494-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8172494-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8172494-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-1st-paste.txt
sed 's/$/BSRR8172494-/g' /mnt/data/gene-alignments/grep-scripts-SRR8172494-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8172494-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8172494-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8172494-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8172494-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8172494-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8146022-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8146022-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8146022-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-1st-paste.txt
sed 's/$/BSRR8146022-/g' /mnt/data/gene-alignments/grep-scripts-SRR8146022-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8146022-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8146022-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8146022-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8146022-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8146022-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8145949-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8145949-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145949-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-1st-paste.txt
sed 's/$/BSRR8145949-/g' /mnt/data/gene-alignments/grep-scripts-SRR8145949-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145949-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8145949-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8145949-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8145949-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145949-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8145959-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8145959-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145959-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-1st-paste.txt
sed 's/$/BSRR8145959-/g' /mnt/data/gene-alignments/grep-scripts-SRR8145959-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145959-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8145959-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8145959-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8145959-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145959-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8582564-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8582564-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8582564-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-1st-paste.txt
sed 's/$/BSRR8582564-/g' /mnt/data/gene-alignments/grep-scripts-SRR8582564-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8582564-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8582564-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8582564-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8582564-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8582564-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8334234-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8334234-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8334234-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-1st-paste.txt
sed 's/$/BSRR8334234-/g' /mnt/data/gene-alignments/grep-scripts-SRR8334234-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8334234-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8334234-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8334234-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8334234-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334234-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8334280-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8334280-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8334280-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-1st-paste.txt
sed 's/$/BSRR8334280-/g' /mnt/data/gene-alignments/grep-scripts-SRR8334280-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8334280-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8334280-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8334280-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8334280-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8334280-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8175167-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8175167-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175167-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-1st-paste.txt
sed 's/$/BSRR8175167-/g' /mnt/data/gene-alignments/grep-scripts-SRR8175167-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175167-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8175167-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8175167-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8175167-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175167-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8145981-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8145981-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145981-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-1st-paste.txt
sed 's/$/BSRR8145981-/g' /mnt/data/gene-alignments/grep-scripts-SRR8145981-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8145981-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8145981-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8145981-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8145981-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8145981-final.txt

sed 's/^/grep -A 1 /g' /mnt/data/gene-alignments/list-to-search.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-front.txt
sed 's/$/ BmntBdataBgene-alignmentsBSRR8175169-all-genes.fa > BmntBdataBgene-alignmentsB/g' /mnt/data/gene-alignments/grep-scripts-SRR8175169-front.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-middle.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175169-middle.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-1st-paste.txt
sed 's/$/BSRR8175169-/g' /mnt/data/gene-alignments/grep-scripts-SRR8175169-1st-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-hind.txt
paste /mnt/data/gene-alignments/grep-scripts-SRR8175169-hind.txt /mnt/data/gene-alignments/list-of-genes.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-2nd-paste.txt
sed 's/$/.fa/g' /mnt/data/gene-alignments/grep-scripts-SRR8175169-2nd-paste.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-end.txt
tr "B" "/" < /mnt/data/gene-alignments/grep-scripts-SRR8175169-end.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-trB.txt
tr -d "\t" < /mnt/data/gene-alignments/grep-scripts-SRR8175169-trB.txt > /mnt/data/gene-alignments/grep-scripts-SRR8175169-final.txt

# after that, we need to concatenate all final lists together. For that, we need to create a cat script:

sed 's/^/cat grep-scripts-/g' /mnt/data/gene-alignments/list-of-samples.txt > /mnt/data/gene-alignments/cat-scripts-front.txt
sed 's/$/-final.txt >> BmntBdataBgene-alignmentsBall-grep-scripts.sh/g' /mnt/data/gene-alignments/cat-scripts-front.txt > /mnt/data/gene-alignments/cat-scripts-hind.txt
tr "B" "/" < /mnt/data/gene-alignments/cat-scripts-hind.txt > /mnt/data/gene-alignments/cat-scripts-B.sh

#######

# that is end of the first part. It is already too much code, so now save it and try to run it.

exit

# save as get-mitos-genes-I-2022-01-28.sh, run ./chmod.sh and run the script
