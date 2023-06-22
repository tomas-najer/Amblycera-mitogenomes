library(ape)
#library(phyloch)
library(strap)
library(phylotate)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(ggtree)
library(caper)

setwd('/home/tom/Documents/Amblycera-mt-genomes/phylogenetic_signal')
louse_tree <- read.nexus("Amblycera_concatenated.fasta.timetree.nex")

plot.phylo(louse_tree,cex=0.5)
nodelabels()
louse_tree <- extract.clade(phy=louse_tree, 123)
?plot.phylo
plot.phylo(louse_tree,cex=0.5)
louse_tree <- extract.clade(phy=louse_tree, 142)
plot.phylo(louse_tree,cex=0.5)
louse_tree <- drop.tip(louse_tree, c("Risp","RmspRatuc","FrspCaban"))
x <- read.csv("labels_in_order.csv")
louse_tree <- rename_taxa(louse_tree, x, label, label2)
plot.phylo(louse_tree,cex=0.5)

arch <- read.csv('arch.csv',header=T)
colnames(arch) <- c("Sequence_ID","Structure")
arch
full_data <- comparative.data(louse_tree,arch, Sequence_ID, na.omit = F)
sig_full <- phylo.d(full_data, binvar = Structure)
summary(sig_full)
print(sig_full)
plot(sig_full)
?comparative.data
?phylo.d
