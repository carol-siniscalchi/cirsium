setwd("~/cirsium/traits/")
getwd()

library(ggplot2)
library(ggpubr)
library(phytools)
library(caper)
library(geiger)

tree = read.tree("cirsium.phyluce95p.ultrametric.tre")
tree$node.label <- NULL

X <- read.csv("traits.csv", row.names=1)
data <- data.frame(X)

CorollaThroat<-as.factor(setNames(X[,1],rownames(X)))
CorollaColor<-as.factor(setNames(X[,2],rownames(X)))
Duration<-as.factor(setNames(X[,3],rownames(X)))
GlutinousRidge<-as.factor(setNames(X[,4],rownames(X)))
Habit<-as.factor(setNames(X[,5],rownames(X)))
HeadPosition<-as.factor(setNames(X[,6],rownames(X)))
InnerPhyllaryColor<-as.factor(setNames(X[,7],rownames(X)))
OuterPhyllaryPosition<-as.factor(setNames(X[,8],rownames(X)))
StyleBranchColor<-as.factor(setNames(X[,9],rownames(X)))

co <- c("#7fc97f","#beaed4","#fdc086","#ffff99")
co <-c("blueviolet", "gold", "hotpink", "firebrick1")

## Fit discrete traits and get AIC
# Corolla Throat
#ER
er.CoTh <- ace(CorollaThroat, tree, type="d", model="ER")
aic.er.CoTh <- AIC(er.CoTh)

#ARD
ard.CoTh <- ace(CorollaThroat, tree, type="d", model="ARD")
aic.ard.CoTh <- AIC(ard.CoTh)

# CorollaColor
#ER
er.CoCo <- ace(CorollaColor, tree, type="d", model="ER")
aic.er.CoCo <- AIC(er.CoCo)

#ARD
ard.CoCo <- ace(CorollaColor, tree, type="d", model="ARD")
aic.ard.Coco <- AIC(ard.CoCo)

# Duration
#ER
er.Du <- ace(Duration, tree, type="d", model="ER")
aic.er.Du <- AIC(er.Du)

#ARD
ard.Du <- ace(Duration, tree, type="d", model="ARD")
aic.ard.Du <- AIC(ard.Du)

# GlutinousRidge
#ER
er.GuRi <- ace(GlutinousRidge, tree, type="d", model="ER")
aic.er.GuRi <- AIC(er.GuRi)

#ARD
ard.GuRi <- ace(GlutinousRidge, tree, type="d", model="ARD")
aic.ard.GuRi <- AIC(ard.GuRi)

# Habit
#ER
er.Ha <- ace(Habit, tree, type="d", model="ER")
aic.er.Ha <- AIC(er.Ha)

#ARD
ard.Ha <- ace(Habit, tree, type="d", model="ARD")
aic.ard.Ha <- AIC(ard.Ha)

# HeadPosition
#ER
er.HePo <- ace(HeadPosition, tree, type="d", model="ER")
aic.er.HePo <- AIC(er.HePo)

#ARD
ard.HePo <- ace(HeadPosition, tree, type="d", model="ARD")
aic.ard.HePo <- AIC(ard.HePo)

# InnerPhyllaryColor
#ER
er.IPhCo <- ace(InnerPhyllaryColor, tree, type="d", model="ER")
aic.er.IPhCo <- AIC(er.IPhCo)

#ARD
ard.IPhCo <- ace(InnerPhyllaryColor, tree, type="d", model="ARD")
aic.ard.IPhCo <- AIC(ard.IPhCo)

# OuterPhyllaryPosition
#ER
er.OPhPo <- ace(OuterPhyllaryPosition, tree, type="d", model="ER")
aic.er.OPhPo <- AIC(er.OPhPo)

#ARD
ard.OPhPo <- ace(OuterPhyllaryPosition, tree, type="d", model="ARD")
aic.ard.OPhPo <- AIC(ard.OPhPo)

# StyleBranchColor
#ER
er.StBrCo <- ace(StyleBranchColor, tree, type="d", model="ER")
aic.er.StBrCo <- AIC(er.StBrCo)

#ARD
ard.StBrCo <- ace(StyleBranchColor, tree, type="d", model="ARD")
aic.ard.StBrCo <- AIC(ard.StBrCo)

### AIC Table ###

aic.ard.CoCo <- 0
trait <- c("CorollaThroat","CorollaColor", "Duration", "GlutinousRidge", 
           "Habit", "HeadPosition", "InnerPhyllaryColor", "OuterPhyllaryPosition", "StyleBranchColor")
ER.model <- c(aic.er.CoTh, aic.er.CoCo, aic.er.Du,  aic.er.GuRi, aic.er.Ha, aic.er.HePo, 
              aic.er.IPhCo, aic.er.OPhPo, aic.er.StBrCo)
ADR.model <- c(aic.ard.CoTh, aic.ard.CoCo, aic.ard.Du,  aic.ard.GuRi, aic.ard.Ha, aic.ard.HePo, 
               aic.ard.IPhCo, aic.ard.OPhPo, aic.ard.StBrCo)
aic.table <- data.frame(trait, ER.model, ADR.model)
write.csv(aic.table, file = "PollinationTraits.aic.table.csv")


#### AIC Table ####

### trait ER.model ADR.model
### 1         CorollaThroat 32.01494  32.31644
### 2          CorollaColor 94.70025  NA
### 3              Duration 81.27504  82.77892
### 4        GlutinousRidge 73.41766  74.32044
### 5                 Habit 26.56948  26.73565
### 6          HeadPosition 37.20885  37.54251
### 7    InnerPhyllaryColor 91.15670 103.69369
### 8 OuterPhyllaryPosition 95.19537  98.71720
### 9      StyleBranchColor 84.50775  99.37997

###### plots ######

co <-c("blueviolet", "gold", "hotpink", "cyan")

pdf(file="traits.rec1.pdf", width = 30, height = 20)
par(mfrow = c(1,3))
plot(tree, cex = 1.5) 
nodelabels(pie=er.CoTh$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(CorollaThroat[tree$tip.label], levels(CorollaThroat)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold"),c("Lobes less than twice as long as throat","Lobes twice as long as throat")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.CoCo$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(CorollaColor[tree$tip.label], levels(CorollaColor)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold", "hotpink", "cyan"),c("red","reddish-purple", "white, pink, purple", "yellow")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.Du$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(Duration[tree$tip.label], levels(Duration)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold"),c("perennial","biennial")),
                  prompt=FALSE,shape="circle",x=0,y=0)
dev.off()

pdf(file="traits.rec2.pdf", width = 30, height = 20)
par(mfrow = c(1,3))
plot(tree, cex = 1.5) 
nodelabels(pie=er.GuRi$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(GlutinousRidge[tree$tip.label], levels(GlutinousRidge)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold"),c("absent","present")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.Ha$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(Habit[tree$tip.label], levels(Habit)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold"),c("caulescent","acaulescent")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.HePo$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(HeadPosition[tree$tip.label], levels(HeadPosition)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold"),c("erect","pending")),
                  prompt=FALSE,shape="circle",x=0,y=0)
dev.off()


pdf(file="traits.rec3.pdf", width = 30, height = 20)
par(mfrow = c(1,3))
plot(tree, cex = 1.5) 
nodelabels(pie=er.IPhCo$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(InnerPhyllaryColor[tree$tip.label], levels(InnerPhyllaryColor)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold", "hotpink", "cyan"),c("green","purple", "red", "pink")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.OPhPo$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(OuterPhyllaryPosition[tree$tip.label], levels(OuterPhyllaryPosition)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold", "hotpink"),c("ascending","spreading", "reflexed")),
                  prompt=FALSE,shape="circle",x=0,y=0)
plot(tree, cex = 1.5)
nodelabels(pie=er.StBrCo$lik.anc, piecol=co, cex=1)
tiplabels(pie=to.matrix(StyleBranchColor[tree$tip.label], levels(StyleBranchColor)),piecol=co,cex=0.5)
add.simmap.legend(colors=setNames(c("blueviolet", "gold", "hotpink", "cyan"),c("red","reddish-purple", "white, pink or purple", "yellow")),
                  prompt=FALSE,shape="circle",x=0,y=0)
dev.off()


### Phylogenetic signal ####
### D statistics -  only binary traits ####

data <- read.csv("traits.csv")
data <- as.data.frame(data)
data_object <- comparative.data(tree, data, Species, vcv=TRUE)

CoTh.d <- phylo.d(data_object, binvar=CorollaThroat)
print(CoTh.d)

# Estimated D :  -4.408385
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.008
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.954

Du.d <- phylo.d(data_object, binvar=Duration)
print(Du.d)

# Estimated D :  -0.2644426
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.072
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.614

GuRi.d <- phylo.d(data_object, binvar=GlutinousRidge)
print(GuRi.d)

# Estimated D :  -3.129107
# Probability of E(D) resulting from no (random) phylogenetic structure :  0
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.986

Ha.d <- phylo.d(data_object, binvar=Habit)
print(Ha.d)

# Estimated D :  1.739055
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.61
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.282

HePo.d <- phylo.d(data_object, binvar=HeadPosition)
print(HePo.d)

# Estimated D :  -1.200083
# Probability of E(D) resulting from no (random) phylogenetic structure :  0.174
# Probability of E(D) resulting from Brownian phylogenetic structure    :  0.713

