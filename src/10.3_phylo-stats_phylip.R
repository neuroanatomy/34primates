#
# Phylogenetic tree plot and independent contrast using APE
# Roberto Toro, Katja Heuer 2018
#
# We are using R Version 3.5

# To get phylip, download from here http://evolution.genetics.washington.edu/phylip/getme-new1.html 
# and install following the instructions here http://evolution.genetics.washington.edu/phylip/install.html.
# Then open 'contrast.app' to give system permissions, and create a symlink to 'contrast' [ln -s /Applications/phylip-3.695/exe/contrast.app/Contents/MacOS/contrast /Applications/phylip-3.695/exe/contrast]

library(ape)
library(phytools)

# Set working directory
dir="data/derived/stats"
setwd(dir)

outputdir <- "3phylo-phylip"
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

#---------------------------------
# Covariation using varCompPhylip
# (which relies on `contrast` from
#  Phylip)
#---------------------------------

if (Sys.info()[["sysname"]]) {
  y <- list()
  for(i in 1:length(SA)) {
    y[[names(SA[i])]]<-matrix(unlist(c(SA[i],V[i],G[i],L[i],N[i],W[i],D[i])),nrow=length(SA[i][[1]]))
  }
  sink("3phylo-phylip/3.1.varCompPhylip.txt")
  varCompPhylip(y,tree, exec="/Applications/phylip-3.695/exe/contrast") # change to the appropriate path
  sink()
}
