###############################################################
# Framework to predict functional profiles from amplicon data # 
# Data: Tonia DeBellis Maple Root Bacteria Data               #
# Date: 07/06/2019                                            #
###############################################################

#### Prepare your R environment ####

# Clear your R environment
# rm(list=ls())

# Set Working directory
# setwd("~/Desktop/T4F")
# setwd(choose.dir())

#### Install and load required libraries ####

# Install from CRAN

install.packages("ggplot2")
install.packages("vegan")
install.packages("dplyr")
install.packages("picante")
install.packages("tidyr")
install.packages("qiimer")
install.packages("rpart")
install.packages("Hmisc")
install.packages("devtools")

# Install libraries from BiocManager

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomformat", ask = FALSE)
BiocManager::install("DESeq2", ask = FALSE)
BiocManager::install("phyloseq", ask = FALSE)

# Install packages from Github

devtools::install_github("joey711/biom")

# Install Tax4Fun from source file
install.packages("Tax4Fun_0.3.1.tar.gz", 
                 repos = NULL, 
                 type = "source")

# Load libraries

library(ggplot2)
library(vegan)
library(picante)
library(dplyr)
library(tidyr)
library(rpart)
library(Hmisc)
library(qiimer)
library(DESeq2)
library(phyloseq)
library(biomformat)
library(Tax4Fun)
library(biom)

#### Workflow ####
#### Using Tax4Fun and biomformat ####

# Loading OTU, taxonomic and metadata information 
Q.ROOT.bact.otu <- read.csv("Q.ROOTBACTcomm.8kRARESILVA.T4FB.csv", 
                            row.names = 1) 

Q.ROOT.bact.meta <- read.csv("Q.ROOTBACTmeta.T4F.csv", 
                             row.names = 1)[, -1]

Q.ROOT.bact.taxo <- read.csv("Q.ROOTBACTtaxSILVA.T4F.csv", 
                             row.names = 1, 
                             stringsAsFactors = T)[, 1:7]

# Checking data structure
dim(Q.ROOT.bact.otu) # 76 3608
dim(Q.ROOT.bact.meta) # 76 23
dim(Q.ROOT.bact.taxo) # 3608 7

# No columns have sum zero:
# Q.ROOT.bact.otu[colSums(Q.ROOT.bact.otu) == 0]

# Generate Biom File  (biomformat)
# https://rdrr.io/bioc/biomformat/man/make_biom.html 
# data =   OTUs / species are rows, samples / sites  are columns.
# sample metadata = data.frame with the number of rows equal to the number of
# samples in data
# observational metadata = data.frame with the number of rows equal to the number
# of features / species / OTUs / genes in data

rootBact_biom <- biomformat::make_biom(data = t(Q.ROOT.bact.otu), 
                                       sample_metadata = Q.ROOT.bact.meta,
                                       observation_metadata = Q.ROOT.bact.taxo,
                                       id = NULL)

# You can check the classes (and the content) of each part of biom.tax4fun
class(biom_data(rootBact_biom))
class(sample_metadata(rootBact_biom))
class(observation_metadata(rootBact_biom))

head(rootBact_biom)

# Here, export it and import it back with Tax4Fun::importQIIMEBiomData()

outfile = tempfile() # Create a temporary file / or give a directory for your file, try with getwd()
biomformat::write_biom(rootBact_biom, outfile) # Write the biom file to outfile
rootBact_biom.rewritten = biomformat::read_biom(outfile) # Read it using read_biom()

# They are not identical (it makes no sense hahaha)
identical(rootBact_biom, rootBact_biom.rewritten)

# To use Tax4Fun, you have to import the biom file with Tax4Fun::importQIIMEBiomData()
# which has a different structure in comparison to the outfile one

rootBact_biom.tax4fun <- importQIIMEBiomData(outfile)

## infer function from 16S using Silva123 annotations 
rootBact.Tax4Fun.fctProf <- Tax4Fun(rootBact_biom.tax4fun, 
                            "SILVA123", 
                            fctProfiling = TRUE, 
                            shortReadMode = FALSE, 
                            normCopyNo = TRUE)

str(rootBact.Tax4Fun.fctProf)

rootBact.Tax4Fun.fctProf$FTU

rootBact.Tax4Fun.fctProf.out <- data.frame(rootBact.Tax4Fun.fctProf$Tax4FunProfile)
dim(rootBact.Tax4Fun.fctProf.out) # 6393 orthologs across 76 samples

## exploring metabolic pathways

rootBact.Tax4Fun.metaPath <- Tax4Fun(rootBact_biom.tax4fun, 
                                    "SILVA123",
                                    fctProfiling = FALSE, # look for metabolic pathways
                                    shortReadMode = FALSE, 
                                    normCopyNo = TRUE)

rootBact.Tax4Fun.metaPath$FTU

rootBact.Tax4Fun.metaPath.out <- data.frame(rootBact.Tax4Fun.metaPath$Tax4FunProfile)

dim(rootBact.Tax4Fun.metaPath.out) # 279 pathways across 76 samples

######################################################################

head(rootBact.Tax4Fun.fctProf.out)

######################################################################

Q.ROOT.bact.meta[, 1, drop = FALSE]

# PERMANOVA
m <- adonis2(decostand(rootBact.Tax4Fun.fctProf.out, method="hellinger") ~ StandType, data = Q.ROOT.bact.meta[, 1, drop = FALSE])
m$aov.tab

# PERMANOVA
m <- adonis(decostand(rootBact.Tax4Fun.metaPath.out, method="hellinger") ~ StandType, data = Q.ROOT.bact.meta[, 1, drop = FALSE])
m$aov.tab

# 
div <- diversity(rootBact.Tax4Fun.fctProf.out, "simpson")
boxplot(div ~ Q.ROOT.bact.meta[, 1, drop = FALSE]$StandType)
oneway.test(div ~ Q.ROOT.bact.meta[, 1, drop = FALSE]$StandType, var.equal = T)
aov.out <- aov(div ~  Q.ROOT.bact.meta[, 1, drop = FALSE]$StandType)
summary(aov.out)
TukeyHSD(aov.out)


# How functionally rich the sites are? E: More functional alpha diversity as you go south.
rootBact.Tax4Fun.fctProf.metaMDS <- metaMDS(decostand(rootBact.Tax4Fun.fctProf.out, method = "hellinger"))

stressplot(rootBact.Tax4Fun.fctProf.metaMDS)

ordiplot(rootBact.Tax4Fun.fctProf.metaMDS, 
         display="sites", 
         type="text")

ordihull(rootBact.Tax4Fun.fctProf.metaMDS,
         groups = Q.ROOT.bact.meta[, 1, drop = FALSE]$StandType,
         draw = "polygon",
         col = 1:3,
         label = T)

rootBact.Tax4Fun.metaPath.metaMDS <- metaMDS(decostand(rootBact.Tax4Fun.metaPath.out, method = "hellinger"))
ordiplot(rootBact.Tax4Fun.metaPath.metaMDS, 
         display="sites", 
         type="text")



# How functionally unique sites are? E:
rootBact.Tax4Fun.fctProf.sorensen <- betadiver(decostand(rootBact.Tax4Fun.metaPath.out, method = "hellinger"), "w")

rootBact.Tax4Fun.fctProf.BrayDis <- vegdist(decostand(rootBact.Tax4Fun.fctProf.out, "hellinger"), "bray")


# run without copy number correction
biom.tax4fun.test.uncorr <- Tax4Fun(tax4fun.test.biom, "SILVA123", shortReadMode = FALSE, normCopyNo = FALSE)

## summarize
biom.tax4fun.test.uncorr

## extract gene relative abundances for samples
gene.abunds <- apply(rootBact.Tax4Fun.fctProf$Tax4FunProfile, 2, sum)
gene.abunds.uncorr <- apply(biom.tax4fun.test.uncorrr$Tax4FunProfile,2,sum)

## write results to csv files
write.csv(biom.tax4fun$Tax4FunProfile, "Q.rootBact.Tax4Fun.fctProf.test.gene.abunds.corr.csv")
write.csv(biom.tax4fun.test.uncorrr$Tax4FunProfile, "Q.ROOTbact.test.tax4fun.gene.abunds.uncorr.csv")

save.image(file= "tax4fun.test.RData")

metaMDS.gene.abund <- metaMDS(decostand(gene.abunds, method = "hellinger"))
ordiplot()

########################
# Forthcoming analyses #
########################

# BetaDiv distances

# LCBD: Functional and Taxonomic Composition 
# SCBD: Functional and Taxonomic Composition

# permanova - nestedness and random effect

# Plots

# ~ Latitude

