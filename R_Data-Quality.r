
# R Script for the analysis of the quality of metagenome samples of the DOK rainout-trial 2022
# Elena Kost
# 23.08.24

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, plyr, Nonpareil)
library(pacman)



# Import ontologies data ==============================================================================================================================================================================================================================================

seed.c <- read.delim("SEED_count.tsv",sep="\t")
colnames(seed.c) <- sub("S\\.", "S-", colnames(seed.c))
str(seed.c)
seed.c[,-c(1:5)][is.na(seed.c[,-c(1:5)])] <- 0

eggnog.c <- read.delim("EGGNOG_count.tsv",sep="\t")
colnames(eggnog.c) <- sub("S\\.", "S-", colnames(eggnog.c))
str(eggnog.c)
eggnog.c[,-c(1:4)][is.na(eggnog.c[,-c(1:4)])] <- 0

inter.c <- read.delim("INTERPRO2GO_count.tsv",sep="\t")
colnames(inter.c) <- sub("S\\.", "S-", colnames(inter.c))
str(inter.c)
inter.c[,-c(1:9)][is.na(inter.c[,-c(1:9)])] <- 0

ec.c <- read.delim("EC_count.tsv",sep="\t")
colnames(ec.c) <- sub("S\\.", "S-", colnames(ec.c))
str(ec.c)
ec.c[,-c(1:5)][is.na(ec.c[,-c(1:5)])] <- 0

cazy.c <- read.delim("CAZy_count.tsv",sep="\t")
colnames(cazy.c) <- sub("S\\.", "S-", colnames(cazy.c))
str(cazy.c)
cazy.c[,-c(1:3)][is.na(cazy.c[,-c(1:3)])] <- 0

ncyc.c <- read.delim("NCyc_count.tsv",sep="\t")
colnames(ncyc.c) <- sub("S\\.", "S-", colnames(ncyc.c))
str(ncyc.c)
ncyc.c[,-c(1:3)][is.na(ncyc.c[,-c(1:3)])] <- 0

pcyc.c <- read.delim("PCyc_count.tsv",sep="\t")
colnames(pcyc.c) <- sub("S\\.", "S-", colnames(pcyc.c))
str(pcyc.c)
pcyc.c[,-c(1:3)][is.na(pcyc.c[,-c(1:3)])] <- 0



# Check annotation success ==============================================================================================================================================================================================================================================

nrow(seed.c[,-c(1:5)]) # number of genes: 789 genes
sum(seed.c[,-c(1:5)]) # annotated reads: 1'037'551'356 reads

nrow(eggnog.c[,-c(1:4)]) # number of genes: 2'681 genes
sum(eggnog.c[,-c(1:4)]) # annotated reads: 500'119'176 reads

nrow(inter.c[,-c(1:9)]) # number of genes: 11'994 genes
sum(inter.c[,-c(1:9)]) # annotated reads: 1'189'490'242 reads

nrow(ec.c[,-c(1:5)]) # number of genes: 4'311 genes
sum(ec.c[,-c(1:5)]) # annotated reads: 1'293'766'860 reads

nrow(cazy.c[,-c(1:3)]) # number of genes: 778 genes
sum(cazy.c[,-c(1:3)]) # annotated reads: 515'018'928 reads

nrow(ncyc.c[,-c(1:3)]) # number of genes: 69
sum(ncyc.c[,-c(1:3)]) # annotated reads: 46'928'906

nrow(pcyc.c[,-c(1:3)]) # number of genes: 139
sum(pcyc.c[,-c(1:3)]) # annotated reads: 180'919'530



# Import design file ==============================================================================================================================================================================================================================================

design <- read.table("design_metagenomic.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE,
                     colClasses = "factor")
str(design)



# Analysis of rowsums from data files using anova ==============================================================================================================================================================================================================================================

# SEED

model.seed <- aov(rowSums(as.data.frame(t(seed.c[,-c(1:5)]))) ~ Irrigation * Treatment * Date, data = design)
summary(model.seed) # no effect
shapiro.test(residuals(model.seed))
leveneTest(model.seed)

# eggNOG

model.eggnog <- aov(log(rowSums(as.data.frame(t(eggnog.c[,-c(1:4)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.eggnog) # no effect
shapiro.test(residuals(model.eggnog))
leveneTest(model.eggnog)

# INTERPRO2GO

model.inter <- aov(log(rowSums(as.data.frame(t(inter.c[,-c(1:9)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.inter) # no effect
shapiro.test(residuals(model.inter))
leveneTest(model.inter)

# EC

model.ec <- aov(log(rowSums(as.data.frame(t(ec.c[,-c(1:5)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.ec) # no effect
shapiro.test(residuals(model.ec))
leveneTest(model.ec)

# CAZyme
model.cazy <- aov(log(rowSums(as.data.frame(t(cazy.c[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.cazy) # no effect
shapiro.test(residuals(model.cazy))
leveneTest(model.cazy)

# NCycDB
model.ncyc <- aov(log(rowSums(as.data.frame(t(ncyc.c[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.ncyc) # no effect
shapiro.test(residuals(model.ncyc))
leveneTest(model.ncyc)

# PCycDB
model.pcyc <- aov(log(rowSums(as.data.frame(t(pcyc.c[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.pcyc) # no effect
shapiro.test(residuals(model.pcyc))
leveneTest(model.pcyc)



# Check sequencing coverage ==============================================================================================================================================================================================================================================

import.npo.processed <- list.files(pattern = "\\.npo$", full.names = T)
nonpareil.set <- Nonpareil.set(import.npo.processed, plot = FALSE, star = 0)
nonpareil.curve <- lapply(import.npo.processed, Nonpareil.curve)

# Add dashed line for coverage
plot(nonpareil.curve[[1]], col = 1,lwd = 2, xlim = c(1000000000,500000000000),
     main = "Nonpareil Curves preprocessed", xlab = "Sequencing effort (bp)", ylab = "Estimated Average Coverage")

for (i in 2:length(nonpareil.curve)) {
  plot(nonpareil.curve[[i]], col = i, lwd = 2, xlim = c(1000000000,500000000000), add = TRUE)
  abline(h = nonpareil.curve[[i]]$C, col = i, lty = 2)
}

Nonpreil.processed.plot 



# Check sequencing coverage range ==============================================================================================================================================================================================================================================

coverage_values <- sapply(nonpareil.curve, function(x) x$C)
min(coverage_values)
max(coverage_values)


