#
# Statistics and visualisation of the raw data
# Roberto Toro, Katja Heuer 2018
#

# Install required packages
install.packages("corrplot",repos = "http://cran.us.r-project.org")
# and then enable corrplot in your 'Packages' tab if you use RStudio
# We are using R Version 3.5, corrplot 0.84

require(corrplot)

# Set working directory
dir="data/derived/stats"
setwd(dir)

# Load lookup table for matching phenotypes and Wilson&Reeder names
# Use the following line to generate figure with some important names for reference 
#lut <- read.csv("../../../src/10kTrees_34PrimateSpeciesGroups2.tsv", sep='\t')
# Use the following line to generate figure without name labels
lut <- read.csv("../../../src/10kTrees_34PrimateSpecies.tsv", sep='\t')

print(lut)

# Load neuroanatomical phenotypes
# and append Wilson&Reeder names
pheno <- read.csv("stats.csv", sep="\t")
pheno <- merge(pheno,lut,by="SpecimenID")

# Keep folding length not logarithmic
FL <- pheno$Folding.Length
# Move measurements to logarithms (except AbsGI)
pheno$Surface.Area <- log10(pheno$Surface.Area)
pheno$Volume <- log10(pheno$Volume)
# pheno$AbsGI <- log10(pheno$AbsGI)
pheno$Folding.Length <- log10(pheno$Folding.Length)
pheno$Folding.number <- log10(pheno$Folding.number)
pheno$Delta[pheno$Delta<0] <- 0

# Create output directory
dir.create('1neuroanat', showWarnings=FALSE)

# Plot matrix of scatterplots of the raw data
panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=10)
  breaks <- h$breaks;
  nB <- length(breaks)
  y <- h$counts;
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

panel.cor <- function(x, y) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  test <- cor.test(x,y)
  rlow <- round(test$conf.int[1], digits=2)
  rup <- round(test$conf.int[2], digits=2)
  #p <- format.pval(test$p.value, eps=1e-6)
  txt <- paste0(r, "\n[", rlow, ", ",rup, "]")
  text(0.5, 0.5, txt, cex = 1)
}

# Save figure of scatterplots
pdf('1neuroanat/1.1.scatterplots.pdf')
pairs(~Surface.Area+Volume+AbsGI+Folding.Length+Folding.number+Lambda+Delta,
      data=pheno,
      labels=c("Surface Area","Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth"),
      panel=panel.smooth)
graphics.off()

# Plot & save correlation matrix
M <- cbind(pheno$Surface.Area,pheno$Volume,pheno$AbsGI,pheno$Folding.Length,pheno$Folding.number,pheno$Lambda,pheno$Delta)
R <- cor(M)
rownames(R) <- colnames(R) <- c("Surface Area","Volume","AbsGI","Folding Length","Folding Number","Fold Wavelength","Fold Depth")
pdf('1neuroanat/1.2.correlations.pdf')
corrplot(R, method="ellipse", order="AOE")
dev.off()

V <- pheno$Volume

# Color the dots in the plot according to family
colour <- list()
colour["Galagonidae"] <- "#a5b7b77f"
colour["Loridae"] <- "#efe0ce7f"
colour["Lemuriformes"] <- "#b69b9d7f"
colour["Cebidae"] <- "#f4e9be7f"
colour["Atelidae"] <- "#caddef7f"
colour["Hominoidea"] <- "#e9c3927f"
colour["Hominoidea_Human"] <- "#d9c3927f"
colour["Hominoidea_Chimp"] <- "#ffc3927f"
colour["Colobinae"] <- "#c1cfb67f"
colour["Papionini"] <- "#939db77f"
colour["Cercopithecini"] <- "#ab91817f"
colour["red"] <- "#ff00007f"

# Add family info to each specimen
group <- list()
group["51058"] <- "Hominoidea_Human"
group["51060"] <- "Hominoidea_Human"
group["51061"] <- "Hominoidea_Human"
group["51062"] <- "Hominoidea_Human"
group["51063"] <- "Hominoidea_Human"
group["51114"] <- "Hominoidea_Human"
group["51116"] <- "Hominoidea_Human"
group["51118"] <- "Hominoidea_Human"
group["51152"] <- "Hominoidea_Human"
group["51153"] <- "Hominoidea_Human"
group["Abby_Yerkes"] <- "Hominoidea_Chimp"
group["Amanda_Yerkes"] <- "Hominoidea_Chimp"
group["ApelleCebusApell_f4b9"] <- "Cebidae"
group["Artemus_Yerkes"] <- "Hominoidea_Chimp"
group["Arthur_Yerkes"] <- "Hominoidea_Chimp"
group["Atele_6717"] <- "Atelidae"
group["AteleAtelesAter_7245"] <- "Atelidae"
group["AyeAyeDaubentoni_cb56"] <- "Lemuriformes"
group["Azalea_Yerkes"] <- "Hominoidea_Chimp"
group["BabouinCercopith_a723"] <- "Papionini"
group["Barbara_Yerkes"] <- "Hominoidea_Chimp"
group["Barney_Yerkes"] <- "Hominoidea_Chimp"
group["Bonobo_brian"] <- "Hominoidea"
group["Carl_Yerkes"] <- "Hominoidea_Chimp"
group["ColobeColobusPol_6118"] <- "Colobinae"
group["David_Yerkes"] <- "Hominoidea_Chimp"
group["DouroucouliAotus_091c"] <- "Cebidae"
group["GalagoDemidoviiG_6cf3"] <- "Galagonidae"
group["Gibbon_buddy"] <- "Hominoidea"
group["Gorilla_kinyani"] <- "Hominoidea"
group["GorillaBeringeiG_0854"] <- "Hominoidea"
group["LangurIndienPith_8390"] <- "Colobinae"
group["LemurMongosLemur_b941"] <- "Lemuriformes"
group["LepilemurAQueueR_309f"] <- "Lemuriformes"
group["MacacaFascicular_3a21"] <- "Papionini"
group["MacacaFascicularis_32200"] <- "Papionini"
group["MacacaFascicularis_32202"] <- "Papionini"
group["MacacaFascicularis_32204"] <- "Papionini"
group["MacacaFascicularis_32205"] <- "Papionini"
group["MacacaFascicularis_mf01"] <- "Papionini"
group["MacacaFascicularis_mf02"] <- "Papionini"
group["MacacaMulatta_32198"] <- "Papionini"
group["MacacaMulatta_32199"] <- "Papionini"
group["MacacaMulatta_32201"] <- "Papionini"
group["MacacaMulatta_32203"] <- "Papionini"
group["MacacaMulatta_mm01"] <- "Papionini"
group["MacaqueCrabierMa_8f9e"] <- "Papionini"
group["MacaqueRhesusMac_bdf8"] <- "Papionini"
group["MakiCattaLemurCa_dec8"] <- "Lemuriformes"
group["MangabeyCercoceb_4458"] <- "Papionini"
group["MangabeyCouronne_6b65"] <- "Papionini"
group["MarsousetOedipom_e834"] <- "Cebidae"
group["MicrocebeDeCoque_e3a9"] <- "Lemuriformes"
group["MicrocebeMignonM_ca6c"] <- "Lemuriformes"
group["MoustacLasiopygi_23cc"] <- "Cercopithecini"
group["NycticebusTardig_08cb"] <- "Loridae"
group["Orangoutan_80bb"] <- "Hominoidea"
group["OuistitiAPinceau_b4aa"] <- "Cebidae"
group["PithecusGermaini_6a2c"] <- "Colobinae"
group["SaimiriCassiquia_7176"] <- "Cebidae"
group["SaimiriCassiquia_eb44"] <- "Cebidae"
group["SapajouCapucinCe_d971"] <- "Cebidae"
group["SingeLaineuxLago_0f84"] <- "Atelidae"
group["VariNoirEtBlancL_4d0e"] <- "Lemuriformes"
group["VervetCercopithe_72ae"] <- "Cercopithecini"

col=unlist(colour[unlist(group[pheno$SpecimenID])])

loessCI <- function(x, y, myTitle, expx=FALSE, expy=FALSE) {
  plx <- predict(loess(y~x, span=0.7), se=TRUE)
  a <- order(x)
  if(expx) {
    x <- 10^x
  }
  plot(x,y,pch=19,col=col,cex=3)
  text(x,y,pheno$PlotName,cex=0.5)
  lines(x[a],plx$fit[a])
  lines(x[a],plx$fit[a]-qt(0.975,plx$df)*plx$se[a], lty=2)
  lines(x[a],plx$fit[a]+qt(0.975,plx$df)*plx$se[a], lty=2)
  title(myTitle)
}

# Save figures
pdf('1neuroanat/1.3.wavelength-vs-logVolume.pdf')
loessCI(V, pheno$Lambda, "Fold Wavelength versus log10(Volume)")
graphics.off()
pdf('1neuroanat/1.4.wavelength-vs-Volume.pdf')
loessCI(V, pheno$Lambda, "Fold Wavelength versus Volume", expx=TRUE)
graphics.off()

pdf('1neuroanat/1.5.depth-vs-logVolume.pdf')
loessCI(V, pheno$Delta, "Fold Depth versus log10(Volume)")
graphics.off()
pdf('1neuroanat/1.6.depth-vs-Volume.pdf')
loessCI(V, pheno$Delta, "Fold Depth versus Volume", expx=TRUE)
graphics.off()

pdf('1neuroanat/1.7.length-vs-logVolume.pdf')
loessCI(V, pheno$Folding.Length, "log10(Folding Length) versus log10(Volume)")
graphics.off()
pdf('1neuroanat/1.8.logLength-vs-Volume.pdf')
loessCI(V, pheno$Folding.Length, "log10(Folding Length) versus Volume", expx=TRUE)
graphics.off()
pdf('1neuroanat/1.9.length-vs-Volume.pdf')
loessCI(V, FL, "Folding Length versus Volume", expx=TRUE)
graphics.off()

print(pheno$SpecimenID)
