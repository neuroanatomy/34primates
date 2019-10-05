#-------------------------------
# Fitting phylogenetic models using Rphylopars. Different models of phenotypic evolution are fitted to the data 100 times, and their AIC values are stored (smaller values indicate a better fit)
# Roberto Toro, Katja Heuer 2018
#-------------------------------

# We are using R Version 3.5

library(Rphylopars)
require(phytools)
library(future)
library(listenv)

# Set working directory
dir="data/derived/stats"
setwd(dir)

# Create output directory
outputdir <- "5phylo-rphylopars100"
dir.create(outputdir)

source('../../../src/catn.R')

# Load 10k tree block
print("load tree block")
treeblock <<- read.nexus("../../external/10kTrees/WilsonAndReeder/100/TreeBlock_10kTrees_Primates_Version3.nex")

# Load lookup table for matching phenotypes and Wilson&Reeder names
print("load species names")
lut <- read.csv("../../../src/10kTrees_34PrimateSpecies.tsv", sep='\t')

# Load neuroanatomical phenotypes and append Wilson&Reeder names
print("load phenotypes")
pheno <- read.csv("stats.csv", sep="\t")
pheno <- merge(pheno,lut,by="SpecimenID")

#-------------------------------
options(width=as.integer(200))

# Move measurements to logarithms (except AbsGI)
print("log convert when appropriate")
pheno$Surface.Area <- log10(pheno$Surface.Area)
pheno$Volume <- log10(pheno$Volume)
# pheno$AbsGI <- log10(pheno$AbsGI)
pheno$Folding.Length <- log10(pheno$Folding.Length)
pheno$Folding.number <- log10(pheno$Folding.number)
pheno$Delta[pheno$Delta<0] <- 0

# Prepare an appropriate data frame: 1st column "specimen", following columns traits
print("prepare dataframe")
pheno2 <<- pheno[c("WilsonReeder.Name","Surface.Area","Volume","AbsGI","Folding.Length","Folding.number","Lambda","Delta")]
colnames(pheno2)[1] <- "species"

#print(pheno2)


#------------------------------------
# Fit evolutionary models, 100 times
#------------------------------------
print("fit evolutionary model on 100 trees")

# Fit sequentially
fitModels <- listenv()
for(iter in 1:100) {
    print(iter)
    tree <- force.ultrametric( treeblock[[iter]], method="extend")

    p_star <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="star"), error=function(e) 0)
    p_BM <- tryCatch(phylopars(trait_data=pheno2, tree=tree), error=function(e) 0)
    p_OU <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="OU"), error=function(e) 0)
    p_OU_diag <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="mvOU", full_alpha=FALSE), error=function(e) 0)
    p_OU_full <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="mvOU", full_alpha=TRUE), error=function(e) 0)
    p_EB <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="EB"), error=function(e) 0)

    fitModels[[iter]] <- c(iter,
        tryCatch(AIC(p_star), error=function(e) 0),
        tryCatch(AIC(p_BM), error=function(e) 0),
        tryCatch(AIC(p_OU), error=function(e) 0),
        tryCatch(AIC(p_OU_diag), error=function(e) 0),
        tryCatch(AIC(p_OU_full), error=function(e) 0),
        tryCatch(AIC(p_EB), error=function(e) 0),
        tryCatch(p_OU$model$alpha, error=function(e) 0),
        tryCatch(p_EB$model$rate, error=function(e) 0)
    )
}
fitModels <- as.list(fitModels)

# Fit in parallel
# plan(multicore, workers=8)
# fitModels <- listenv()
# for(iter in 1:100) {
#     fitModels[[iter]] %<-% {
#         print("iteration")
#         tree <- force.ultrametric( treeblock[[iter]], method="extend")

#         p_star <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="star"), error=function(e) 0)
#         p_BM <- tryCatch(phylopars(trait_data=pheno2, tree=tree), error=function(e) 0)
#         p_OU <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="OU"), error=function(e) 0)
#         p_OU_diag <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="mvOU", full_alpha=FALSE), error=function(e) 0)
#         p_OU_full <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="mvOU", full_alpha=TRUE), error=function(e) 0)
#         p_EB <- tryCatch(phylopars(trait_data=pheno2, tree=tree, model="EB"), error=function(e) 0)

#         return(c(iter,
#             tryCatch(AIC(p_star), error=function(e) 0),
#             tryCatch(AIC(p_BM), error=function(e) 0),
#             tryCatch(AIC(p_OU), error=function(e) 0),
#             tryCatch(AIC(p_OU_diag), error=function(e) 0),
#             tryCatch(AIC(p_OU_full), error=function(e) 0),
#             tryCatch(AIC(p_EB), error=function(e) 0),
#             tryCatch(p_OU$model$alpha, error=function(e) 0),
#             tryCatch(p_EB$model$rate, error=function(e) 0)
#         ))
#     }
# }
# fitModels <- as.list(fitModels)

sink("5phylo-rphylopars100/100+.csv")
catn("Tree", "Star", "BM", "OU", "OU_diag", "OU_full", "EB", "OU_alpha", "EB_rate", sep=" ")
for(iter in 1:100) {
    catn(fitModels[[iter]])
}
sink()

#------------------------------------
# Get median values
#------------------------------------
print("get median values")
sink("5phylo-rphylopars100/100++.csv") 
results <- read.csv("5phylo-rphylopars100/100+.csv", sep=' ')
catn("Model", "AIC")
catn("Star", median(results$Star))
catn("BM", median(results$BM))
catn("OU", median(results$OU))
catn("OU_diag", median(results$OU_diag))
catn("OU_full", median(results$OU_full))
catn("EB", median(results$EB))
catn("OU_alpha", median(results$OU_alpha))
catn("EB_rate", median(results$EB_rate))
sink()
