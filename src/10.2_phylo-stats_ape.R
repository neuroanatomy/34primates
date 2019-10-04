#
# Phylogenetic tree plot and independent contrast using APE
# Roberto Toro, Katja Heuer 2018
#
# We are using R Version 3.5

library(ape)
library(phytools)
library(corrplot)

# Set working directory
getwd()
dir="data/derived/stats"
setwd(dir)
getwd()

# Create output directory
outputdir <- "2phylo-ape"
dir.create(file.path(dir, outputdir))

# Load 10k tree, force ultrametric by extending branches
tree <- force.ultrametric(read.nexus("../../external/10kTrees/WilsonAndReeder/1000/consensusTree_10kTrees_Primates_Version3.nex"), method="extend")

# Load lookup table for matching phenotypes and Wilson&Reeder names
lut <- read.csv("../../../src/10kTrees_34PrimateSpecies.tsv", sep='\t')

# Load neuroanatomical phenotypes and append Wilson&Reeder names
pheno <- read.csv("stats.csv", sep="\t")
pheno <- merge(pheno,lut,by="SpecimenID")

# Move measurements to logarithms (except AbsGI)
pheno$Surface.Area <- log10(pheno$Surface.Area)
pheno$Volume <- log10(pheno$Volume)
# pheno$AbsGI <- log10(pheno$AbsGI)
pheno$Folding.Length <- log10(pheno$Folding.Length)
pheno$Folding.number <- log10(pheno$Folding.number)
pheno$Delta[pheno$Delta<0] <- 0

# The function `pic.ortho` accepts multiple observations per individual.
# They are grouped by species name, as the tips of the tree
phenolist <- split(pheno,pheno$WilsonReeder.Name)

# Reorder to match tip label order
phenolist <- phenolist[tree$tip.label];

# Get phenotypes
SA <- lapply(phenolist,"[[","Surface.Area")
V <- lapply(phenolist,"[[","Volume")
G <- lapply(phenolist,"[[","AbsGI")
L <- lapply(phenolist,"[[","Folding.Length")
N <- lapply(phenolist,"[[","Folding.number")
W <- lapply(phenolist,"[[","Lambda")
D <- lapply(phenolist,"[[","Delta")

#---------------
# Plot the tree
#---------------

# Plot tree annotated with branch length
pdf('2phylo-ape/2.1.tree.pdf')
plot(tree, cex=0.5)
axisPhylo()
nodelabels(cex=0.5)
title("Branch length (Mya)")
edgelabels(round(tree$edge.length,1), frame="n", adj=c(0.5,-0.2), cex=0.5)
graphics.off()

# Plot tree as variance-covariance matrix
pdf('2phylo-ape/2.2.varcovmatrix.pdf')
heatmap(vcv.phylo(tree))
graphics.off()

#------------------------------------
# Phylogenetic independent contrasts
#------------------------------------

# Compute contrasts
pic.SA <- pic.ortho(SA,tree,intra=TRUE)
pic.V <- pic.ortho(V,tree,intra=TRUE)
pic.G <- pic.ortho(G,tree,intra=TRUE)
pic.L <- pic.ortho(L,tree,intra=TRUE)
pic.N <- pic.ortho(N,tree,intra=TRUE)
pic.W <- pic.ortho(W,tree,intra=TRUE)
pic.D <- pic.ortho(D,tree,intra=TRUE)

# Plot tree annotated with SA and V contrasts
pdf('2phylo-ape/2.3.tree-SA-V.pdf')
plot(tree, cex=0.5)
axisPhylo()
title("Surface area and volume")
nodelabels(round(pic.SA,3),adj=c(-0.05,-0.25), cex=0.5, frame="n")
nodelabels(round(pic.V,3),adj=c(-0.05,1.25),cex=0.5, frame="n")
graphics.off()

# Plot scatterplots for all PICs
pdf('2phylo-ape/2.4.scatterplots-PICs.pdf')
df=data.frame(pic.SA,pic.V,pic.G,pic.L,pic.N,pic.W,pic.D)
colnames(df) <- c("Surface Area","Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth")
pairs(~pic.SA+pic.V+pic.G+pic.L+pic.N+pic.W+pic.D,
      data=df,
      labels=colnames(df),
      panel=panel.smooth)
graphics.off()

# Plot correlation matrix for all PICs
M <- cbind(pic.SA,pic.V,pic.G,pic.L,pic.N,pic.W,pic.D)
R <- cor(M)
rownames(R) <- colnames(R) <- c("Surface Area","Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth")
pdf('2phylo-ape/2.5.correlations.pdf')
corrplot(R, method="ellipse", order="AOE")
graphics.off()

# Compute a linear regression of SA on V
result <- lm(pic.SA ~ pic.V)
sink("2phylo-ape/2.5.lm-SA-V.txt")
summary(result)
sink()

# PICs have expected mean zero: fit a regression with intercept fixed at 0
result <- lm(pic.SA ~ pic.V -1)
sink("2phylo-ape/2.6.lm-SA-V-mean-nointercept.txt")
summary(result)
sink()

print(pheno$SpecimenID)
