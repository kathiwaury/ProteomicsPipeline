# packages.R

#install missing packages from CRAN
packagesList <- c("drake", "dplyr", "reader", "DT", "ggplot2", "plotly", "kableExtra", "knitr", 
  "caret", "randomForest", "BiocManager")
newPackages <- packagesList[!(packagesList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0) {
  install.packages(newPackages)
}

#install and load missing packages from BioConductor
packagesListBC <- c("clusterProfiler",  "DOSE", "org.Hs.eg.db",  "ReactomePA")
newPackagesBC <- packagesListBC[!(packagesListBC %in% installed.packages()[,"Package"])]
if(length(newPackagesBC) > 0) {
  BiocManager::install(newPackagesBC)
}

#load packages
library(drake)
library(caret)
library(clusterProfiler)
library(dplyr)
library(DOSE)
library(DT)
library(ggplot2)
library(kableExtra)
library(knitr)
library(plotly)
library(org.Hs.eg.db)
library(randomForest)
library(ReactomePA)
library(reader)





