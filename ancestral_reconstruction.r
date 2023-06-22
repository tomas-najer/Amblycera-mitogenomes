library(ape)
#library(phyloch)
library(strap)
library(phylotate)
library(MCMCtreeR)
library(phytools)
library(geiger)
library(ggtree)
library(caper)
library(ggplotify)
library(treeio)
library("readxl")
library("ape")
library("ggplot2")
library("phytools")
library("geiger")
library(ggplot2)
library(gggenes)
library(ggfittext)
library(ggrepel)
library(dplyr)
library("ggtree")
library(ggplotify)
library(treeio)

setwd('/home/tomas/Documents/Amblycera_mitogenomes/Figure_2')
louse_tree <- read.nexus("Amblycera_concatenated.fasta.timetree.nex")

plot.phylo(louse_tree,cex=0.5)
nodelabels()
louse_tree <- extract.clade(phy=louse_tree, 123)
louse_tree <- extract.clade(phy=louse_tree, 142)
louse_tree <- drop.tip(louse_tree, "Risp")
plot.phylo(louse_tree,cex=0.5)

louse_tree$root.time <- 82.9852 #root age
plot.phylo(louse_tree,cex=0.5,no.margin=T)
x <- read.csv("tree_labels_1.csv")
louse_tree <- rename_taxa(louse_tree, x, label, label2)
plot.phylo(louse_tree,cex=0.5)

geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)

genomes_structure <- read.csv('arch.csv',header=T)
genomes_structure
arch <- setNames(as.factor(genomes_structure[,2]),genomes_structure[,1])
arch
to.matrix(arch,levels(arch))
name.check(louse_tree,arch)
structure_matrix<-to.matrix(arch[louse_tree$tip.label],levels(arch))
structure_matrix
geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)
plot.phylo(louse_tree,cex=0.5)
tiplabels(pie=structure_matrix,piecol=c("red","darkgreen"),cex=0.15)
legend("topleft",levels(arch),pch=21,pt.bg=c("red","darkgreen"),
       pt.cex=2.2)

ace(arch,louse_tree,model='ER',type='discrete')
fitER <- ace(arch,louse_tree,model='ER',type='discrete')
ace(arch,louse_tree,model='ARD',type='discrete')
fitARD <- ace(arch,louse_tree,model='ARD',type='discrete')
ace(arch,louse_tree,model='SYM',type='discrete')
fitSYM <- ace(arch,louse_tree,model='SYM',type='discrete')
fitUSR <- ace(arch, louse_tree, model=matrix(c(0,1,0,0),2),type='discrete')
model_geiger_ER <- fitDiscrete(louse_tree,arch,type='discrete',model="ER")
model_geiger_ARD <- fitDiscrete(louse_tree,arch,type='discrete',model="ARD")
model_geiger_SYM <- fitDiscrete(louse_tree,arch,type='discrete',model="SYM")
model_geiger_USR <- fitDiscrete(louse_tree,arch,type='discrete',model=matrix(c(0,1,0,0)))

ERvSYM <- abs(model_geiger_ER$opt$aicc - model_geiger_SYM$opt$aicc)
ERvARD <- abs(model_geiger_ER$opt$aicc - model_geiger_ARD$opt$aicc)
ERvUSR <- abs(model_geiger_ER$opt$aicc - model_geiger_USR$opt$aicc)
SYMvARD <- abs(model_geiger_SYM$opt$aicc - model_geiger_ARD$opt$aicc)

comp <- data.frame(model = c("ER","SYM","ARD","USR"),AIC = c(model_geiger_ER$opt$aic,model_geiger_SYM$opt$aic,model_geiger_ARD$opt$aic,model_geiger_USR$opt$aic),
                   AICc = c(model_geiger_ER$opt$aicc,model_geiger_SYM$opt$aicc,model_geiger_ARD$opt$aicc,model_geiger_USR$opt$aicc))
comp

1-pchisq(2*abs(fitER$loglik - fitARD$loglik), 1)
1-pchisq(2*abs(fitER$loglik - fitUSR$loglik), 1)
1-pchisq(2*abs(fitARD$loglik - fitUSR$loglik), 1)

round(fitER$lik.anc,3)
round(fitARD$lik.anc,3)
round(fitUSR$lik.anc,3)

geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-15,81), lwd=3, width=2,quat.rm = T)
tiplabels(pie=structure_matrix,piecol=c("red","darkgreen"),cex=0.15)
nodelabels(node=1:louse_tree$Nnode+Ntip(louse_tree),
           pie=fitER$lik.anc,piecol=c("red","darkgreen"),cex=0.3)
legend("topleft",levels(arch),pch=21,pt.bg=c("red","darkgreen"),
       pt.cex=2)

ERtrees<-make.simmap(louse_tree,arch,model='ER',nsim=1000)
summary(ERtrees)
ad<-summary(ERtrees,plot=T)
geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.6, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-35,81), lwd=3, width=2,quat.rm = T)
tiplabels(pie=structure_matrix,piecol=c("red","darkgreen"),cex=0.20)
nodelabels(pie=ad$ace,piecol=c("red","darkgreen"),cex=0.4)
legend("topleft",levels(arch),pch=21,pt.bg=c("red","darkgreen"),
       pt.cex=2)

mtrees<-make.simmap(louse_tree,arch,model=matrix(c(0,1,0,0),2),nsim=1000)
summary(mtrees)
pd<-summary(mtrees,plot=T)
pd
geoscalePhylo(tree=ladderize(louse_tree,right=F), units=c("Period","Epoch"), boxes="Period",
              cex.tip=0.6, cex.age=0.7, cex.ts=0.7, label.offset=1, x.lim=c(-35,81), lwd=3, width=2,quat.rm = T)
tiplabels(pie=structure_matrix,piecol=c("red","darkgreen"),cex=0.20)
nodelabels(pie=pd$ace,piecol=c("red","darkgreen"),cex=0.4)
legend("topleft",levels(arch),pch=21,pt.bg=c("red","darkgreen"),
       pt.cex=2)
