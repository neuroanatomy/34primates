#
# Trying Rphylopars
# Roberto Toro, Katja Heuer 2018

# We are using R Versions 3.5

library(Rphylopars)
require(phytools)

# Set working directory
dir="data/derived/stats"
setwd(dir)

# Create output directory
outputdir <- "4phylo-rphylopars"
dir.create(file.path(dir, outputdir))

source('../../../src/catn.R')

# Load 10k tree, force ultrametric by extending branches
tree <- force.ultrametric(read.nexus("../../external/10kTrees/WilsonAndReeder/1000/consensusTree_10kTrees_Primates_Version3.nex"), method="extend")

# Load lookup table for matching phenotypes and Wilson&Reeder names
lut <- read.csv("../../../src/10kTrees_34PrimateSpecies.tsv", sep='\t')

# Load neuroanatomical phenotypes and append Wilson&Reeder names
pheno <- read.csv("stats.csv", sep="\t")
pheno <- merge(pheno,lut,by="SpecimenID")

#-------------------------------
options(width=as.integer(200))

# Move measurements to logarithms (except AbsGI)
pheno$Surface.Area <- log10(pheno$Surface.Area)
pheno$Volume <- log10(pheno$Volume)
# pheno$AbsGI <- log10(pheno$AbsGI)
pheno$Folding.Length <- log10(pheno$Folding.Length)
pheno$Folding.number <- log10(pheno$Folding.number)
pheno$Delta[pheno$Delta<0] <- 0

# Prepare an appropriate data frame: 1st column "specimen", following columns traits
pheno2 <- pheno[c("WilsonReeder.Name","Surface.Area","Volume","AbsGI","Folding.Length","Folding.number","Lambda","Delta")]
colnames(pheno2)[1] <- "species"

print(pheno2)

#-------------------------
# Fit evolutionary models
#-------------------------

# Fit the Brownian motion model
#------------------------------------
# BM=completely unconstrained character evolution
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("fit a Brownian Motion model")
p_BM <- phylopars(trait_data=pheno2, tree=tree)
sink("4phylo-rphylopars/4.1.BM.txt")
p_BM
sink()

# ntips Number of species we have phenotypes for
ntips <- length(tree$tip.label)
# Number of nodes inside the tree (ntips-1) + including the tips
nnodes <- ntips + tree$Nnode

# Put into csv the info for all the estimated phenotypes at each internal node (=all ancestors)
catn("save data in csv")
csv <- cbind(p_BM$anc_recon[(ntips+1):nnodes,])
# lower and upper confidence interval
lci <- csv - sqrt(p_BM$anc_var[(ntips+1):nnodes,])*1.96
uci <- csv + sqrt(p_BM$anc_var[(ntips+1):nnodes,])*1.96
colnames(lci) <- paste(colnames(lci), "-95%CI", sep = "")
colnames(uci) <- paste(colnames(uci), "+95%CI", sep = "")
csv <- cbind(csv,lci)
csv <- cbind(csv,uci)
write.csv(csv, file="4phylo-rphylopars/4.2.ancestral.csv")

myplotanc <- function(tree, tips, tipnames, states, mytitle) {
    names(tips) <- tipnames
    obj <- contMap(tree,tips,method="user",anc.states=states, plot="F")
    obj <- setMap(obj, colors=c("blue","white","red"))
    plot.contMap(obj)
    axisPhylo()
    title(mytitle)
}

catn("plot ancestral values")
tipnames <- rownames(p_BM$anc_recon)[1:ntips]

# Evolution of surface area along the tree
pdf("4phylo-rphylopars/4.3.BM-anc-surface-area.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Surface.Area"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Surface.Area"], "Ancestral surface area reconstruction")
graphics.off()

# Ancestral volume
pdf("4phylo-rphylopars/4.4.BM-anc-volume.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Volume"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Volume"], "Ancestral volume reconstruction")
graphics.off()

# Ancestral absGI
pdf("4phylo-rphylopars/4.5.BM-anc-absgi.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"AbsGI"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"AbsGI"], "Ancestral absGI reconstruction")
graphics.off()

# Ancestral folding length
pdf("4phylo-rphylopars/4.6.BM-anc-folding-length.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Folding.Length"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Folding.Length"], "Ancestral folding length reconstruction")
graphics.off()

# Ancestral folding number
pdf("4phylo-rphylopars/4.7.BM-anc-folding-number.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Folding.number"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Folding.number"], "Ancestral folding number reconstruction")
graphics.off()

# Ancestral fold wavelength
pdf("4phylo-rphylopars/4.8.BM-anc-fold-wavelength.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Lambda"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Lambda"], "Ancestral fold wavelength reconstruction")
graphics.off()

# Ancestral fold depth
pdf("4phylo-rphylopars/4.9.BM-anc-fold-depth.pdf")
myplotanc(tree, p_BM$anc_recon[1:ntips,"Delta"], tipnames, p_BM$anc_recon[(ntips+1):nnodes,"Delta"], "Ancestral fold depth reconstruction")
graphics.off()


# Test Pagel's lambda 1 (no transformation of branch lengths) versus 0 (complete star phylogeny, with all tip branches equal in length and all internal branches of length 0): 
# phylogenetic and phenotypic trait variance-covariance matrices for both cases
# ?? scales the tree between a constant-rates model (??=1) to one where every species is statistically independent of every other species in the tree (??=0).
catn("Fit and compare Pagel's lambda model for lambda=1 (BM) and lambda=0 (star)")
p_lambda <- phylopars(trait_data=pheno2, tree=tree)
sink("4phylo-rphylopars/4.10.test-pagel-lamdba=1-vs-lambda=0.txt")
p_lambda
p_star <- phylopars(trait_data=pheno2, tree=tree, model="star")
p_star
catn('Log-likelihood BM', logLik(p_lambda), sep='\t')
catn('Log-likelihood star',logLik(p_star), sep='\t')
chi_square <- as.double(2*(logLik(p_lambda) - logLik(p_star))) # 2*(logLik_alt - logLik_null)
catn('chi_square', chi_square, sep='\t')
degrees_freedom <- p_lambda$npars - p_star$npars # df = difference in number of model parameters
catn('df', degrees_freedom, sep='\t')
p <- pchisq(q = chi_square,df = degrees_freedom,lower.tail = FALSE) # p-value
catn ('p-value', p, sep='\t')
sink()


# Fit the Ornstein-Uhlenbeck model
#------------------------------------
# Univariate
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("Fit Ornstein-Uhlenbeck model (OU, univariate: a single alpha value shared by all traits)")
p_OU <- phylopars(trait_data=pheno2, tree=tree, model = "OU")
sink("4phylo-rphylopars/4.11.OU-univariate.txt")
p_OU
sink()

# Multivariate
# Save phylogenetic and phenotypic trait variance-covariance matrices, variance explained by phylogeny
catn("Fit multivariate OU (a full matrix of alpha values)")
p_mvOU <- phylopars(trait_data =pheno2,tree = tree,model = "mvOU", full_alpha = FALSE)
sink("4phylo-rphylopars/4.12.OU-multivariate-diagonal.txt")
p_mvOU # Multivariate OU model (with alpha as full matrix)
sink()

catn("Fit multivariate OU (a diagonal matrix of alpha values)")
p_mvOU_diag <- phylopars(trait_data = pheno2,tree = tree, model = "mvOU")
sink("4phylo-rphylopars/4.13.OU-multivariate-full.txt")
p_mvOU_diag
sink()


# Fit the Early-Burst model
#------------------------------------
catn("Fit early-burst model (EB)")
p_EB <- phylopars(trait_data = pheno2,tree = tree,model = "EB")
sink("4phylo-rphylopars/4.14.EB.txt")
p_EB # Estimated trait covariance and EB rate parameter
sink()


# Model selection
#------------------------------------
catn("Model selection using AIC")
sink("4phylo-rphylopars/4.15.model-selection.txt")
catn("Brownian motion", AIC(p_BM), sep='\t')
catn("Ornstein-Uhlenbeck, single alpha", AIC(p_OU), sep='\t')
catn("Ornstein-Uhlenbeck, diagonal alpha matrix", AIC(p_mvOU), sep='\t')
catn("Ornstein-Uhlenbeck, full alpha matrix", AIC(p_mvOU_diag), sep='\t')
catn("Early burst", AIC(p_EB), sep='\t')
sink()


# Plot phenograms
#------------------------------------
catn("Plot phenograms (BM model)")
plotpheno<-function(ph,tr,trait,mytitle) {
  t<-node.depth.edgelength(tree)
  h<-p_BM$anc_recon[,trait]
  plot(c(min(t),max(t)*1.1),c(min(h),max(h)),type="n")
  length(tr$edge)
  for(i in 1:nrow(tr$edge)) {
    ed<-tr$edge[i,]
    segments(t[ed[1]],h[ed[1]],t[ed[2]],h[ed[2]])
  }
  for(i in 1:length(tree$tip.label)) {
    text(t[i],h[i],tree$tip.label[i], adj=c(0,0), cex=0.1)
  }
  title(mytitle)
}

pdf("4phylo-rphylopars/4.16.BM-anc-surface-phenogram.pdf")
plotpheno(pheno2,tree,"Surface.Area","Evolution of surface area")
graphics.off()

pdf("4phylo-rphylopars/4.17.BM-anc-volume-phenogram.pdf")
plotpheno(pheno2,tree,"Volume","Evolution of volume")
graphics.off()

pdf("4phylo-rphylopars/4.18.BM-anc-absgi-phenogram.pdf")
plotpheno(pheno2,tree,"AbsGI","Evolution of AbsGI")
graphics.off()

pdf("4phylo-rphylopars/4.19.BM-anc-folding-length-phenogram.pdf")
plotpheno(pheno2,tree,"Folding.Length","Evolution of folding length")
graphics.off()

pdf("4phylo-rphylopars/4.20.BM-anc-folding-number-phenogram.pdf")
plotpheno(pheno2,tree,"Folding.number","Evolution of folding number")
graphics.off()

pdf("4phylo-rphylopars/4.21.BM-anc-wavelength.pdf")
plotpheno(pheno2,tree,"Lambda","Evolution of folding wavelength")
graphics.off()

pdf("4phylo-rphylopars/4.22.BM-anc-fold-depth-phenogram.pdf")
plotpheno(pheno2,tree,"Delta","Evolution of fold depth")
graphics.off()

