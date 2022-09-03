# Based on FiSSE example file

library(ape)
library(phangorn)
library(diversitree)
library(geiger)
library(phytools)


setwd("~/cirsium/traits/")
getwd()

source("traitDependent_functions.R") # Get from https://github.com/macroevolution/fisse

tree = read.tree("cirsium.phyluce95p.ultrametric.tre")

# Check for ultrametricity and if binary; FiSSE has its own function for ultrametricity but this should be faster and fine with slight precision errors 
if(is.binary(tree)) {
	print("TRUE")
	} else {
	tree <- multi2di(tree)
	}

if(is.ultrametric(tree)) {
	print("TRUE")
	} else {
	tree <- force.ultrametric(tree, method="extend")
	}


morpho<- read.csv("cirsium_binary.csv")

## characters 
# head: Head Position: 0 - erect, 1 - pendent
# habit: 0 - acaulescent, 1 - caulescent
# glutinous: PhyllaryGlutinousRidge, 0 - absent, 1 - present
# duration: 0 - perennial, 1 - biennial
# corolla_throat: CorollaLobeToThroatRatio 1 = lobes less than twice as long as throat; 0=lobes twice as long as throat


#throat

throat <- as.numeric(morpho[,2])
names(throat) <- as.character(morpho[,1])
treedata_object <- treedata(tree, throat)
tree <- treedata_object$phy
throat <- throat[tree$tip.label]
res1 <- FISSE.binary(tree, throat)

pval_1tailed <- min(res1$pval, 1-res1$pval)
pval_1tailed
# throat: 0.1768232

#duration

dur <- as.numeric(morpho[,3])
names(dur) <- as.character(morpho[,1])
treedata_object <- treedata(tree, dur)
tree <- treedata_object$phy
dur <- dur[tree$tip.label]
res2 <- FISSE.binary(tree, dur)

pval_1tailed <- min(res2$pval, 1-res2$pval)
pval_1tailed
# duration: 0.2947053

#glutinous

glu <- as.numeric(morpho[,4])
names(glu) <- as.character(morpho[,1])
treedata_object <- treedata(tree, glu)
tree <- treedata_object$phy
glu <- glu[tree$tip.label]
res3 <- FISSE.binary(tree, glu)

pval_1tailed <- min(res3$pval, 1-res3$pval)
pval_1tailed
# glutinous: 0.2827173

#habit

habit <- as.numeric(morpho[,5])
names(habit) <- as.character(morpho[,1])
treedata_object <- treedata(tree, habit)
tree <- treedata_object$phy
habit <- habit[tree$tip.label]
res4 <- FISSE.binary(tree, habit)

pval_1tailed <- min(res4$pval, 1-res4$pval)
pval_1tailed
# habit: 0.4055944

#head

head <- as.numeric(morpho[,6])
names(head) <- as.character(morpho[,1])
treedata_object <- treedata(tree, head)
tree <- treedata_object$phy
head <- head[tree$tip.label]
res5 <- FISSE.binary(tree, head)

pval_1tailed <- min(res5$pval, 1-res5$pval)
pval_1tailed
# habit: 0.3716284

# pollinator
pollinators <- read.csv("pollinator.csv")
poll <- as.numeric(pollinators[,2])
names(poll) <- as.character(pollinators[,1])
treedata_object <- treedata(tree, poll)
tree <- treedata_object$phy
poll <- poll[tree$tip.label]
res6 <- FISSE.binary(tree, poll)

pval_1tailed <- min(res6$pval, 1-res6$pval)
pval_1tailed
#pollinator: 0.06893107

