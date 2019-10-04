# Install required packages
install.packages("corrplot", repos="https://cloud.r-project.org", quiet=TRUE)
install.packages("devtools", repos="https://cloud.r-project.org", quiet=TRUE)

library("ape")
library("phytools")
library("corrplot")
library("devtools")

install_github("r03ert0/Rphylopars",dependencies = TRUE)
# 1
