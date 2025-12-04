
# R Script for Analysis of the metagenome samples of the DOK rainout-trial 2022
# gene level analysis of CAZyme
# Elena Kost
# 12.09.24

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, ape, ade4)
library(pacman)



# Import data ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy <- read.delim("CAZy_RPKmaxORF.tsv",sep="\t")
colnames(cazy) <- sub("S\\.", "S-", colnames(cazy))
str(cazy)
cazy[,-c(1:3)][is.na(cazy[,-c(1:3)])] <- 0
cazy[,c(1:3)][is.na(cazy[,c(1:3)])] <- "unclassified"

cazy <- cazy %>% relocate(accession) %>% 
  column_to_rownames("accession") 

cazy.genes <- cazy[,c(1:2)]
cazy.genes <- cazy.genes %>% add_column(level_0 = c("carbon_cycling"))



# NCycDB =====================

ncyc <- read.delim("NCyc_RPKmaxORF.tsv",sep="\t")
colnames(ncyc) <- sub("S\\.", "S-", colnames(ncyc))
str(ncyc)
ncyc[,-c(1:3)][is.na(ncyc[,-c(1:3)])] <- 0
ncyc[,c(1:3)][is.na(ncyc[,c(1:3)])] <- "unclassified"

ncyc <- ncyc %>% relocate(accession) %>% 
  column_to_rownames("accession") 

ncyc.genes <- ncyc[,c(1:2)]
ncyc.genes <- ncyc.genes %>% add_column(level_0 = c("nitrogen_cycling"))



# PCycDB =====================

pcyc <- read.delim("PCyc_RPKmaxORF.tsv",sep="\t")
colnames(pcyc) <- sub("S\\.", "S-", colnames(pcyc))
str(pcyc)
pcyc[,-c(1:3)][is.na(pcyc[,-c(1:3)])] <- 0
pcyc[,c(1:3)][is.na(pcyc[,c(1:3)])] <- "unclassified"

pcyc <- pcyc %>% relocate(accession) %>% 
  column_to_rownames("accession") 

pcyc.genes <- pcyc[,c(1:2)]
pcyc.genes <- pcyc.genes %>% add_column(level_0 = c("phosphorus_cycling"))



# Import Design ==============================================================================================================================================================================================================================================

design <- read.table("design_metagenomic.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE, colClasses = "factor")
str(design)
identical(rownames(design), colnames(cazy[,-c(1:2)]))

design$IxT <- with(design, interaction(Irrigation, Treatment, lex.order = TRUE, drop = TRUE))



# Normalization of data ==============================================================================================================================================================================================================================================

# rarefying with 100

iters <- 100

# CAZyme =====================

cazy.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(cazy[,-c(1:2)])),  min(rowSums(as.data.frame(t(cazy[,-c(1:2)]))))),  mc.cores = 10)
cazy.array <- laply(cazy.iters, as.matrix)
cazy.iss <- apply(cazy.array, 2:3, mean)


# NCycDB =====================

ncyc.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(ncyc[,-c(1:2)])),  min(rowSums(as.data.frame(t(ncyc[,-c(1:2)]))))),  mc.cores = 10)
ncyc.array <- laply(ncyc.iters, as.matrix)
ncyc.iss <- apply(ncyc.array, 2:3, mean)


# PCycDB =====================

pcyc.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(pcyc[,-c(1:2)])),  min(rowSums(as.data.frame(t(pcyc[,-c(1:2)]))))),  mc.cores = 10)
pcyc.array <- laply(pcyc.iters, as.matrix)
pcyc.iss <- apply(pcyc.array, 2:3, mean)




# Aggregate data  ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy.prop <- prop.table(cazy.iss, margin = 1) * 100  # get proportions

cazy.level1 <- aggregate(t(cazy.prop) ~ level_1, data = cazy.genes, FUN = sum)  
cazy.level1$sort <- cazy.level1[, 1]
cazy.level1 <- data.frame(cazy.level1[, -1], row.names = paste("l1:", cazy.level1[, 1], sep = ""), check.names = FALSE)

cazy.level2 <- aggregate(t(cazy.prop) ~ level_1 + level_2, data = cazy.genes, FUN = sum)
cazy.level2 <- subset(cazy.level2, cazy.level2$level_2 != "unclassified")
cazy.level2$sort <- paste(cazy.level2[, 1], cazy.level2[, 2], sep = ";")
cazy.level2 <- data.frame(cazy.level2[, -c(1:2)], row.names = paste("_l2:", cazy.level2[, 2], sep = ""), check.names = FALSE)

cazy.gene.summary <- rbind(cazy.level1, cazy.level2)
cazy.gene.summary <- cazy.gene.summary[order(as.character(cazy.gene.summary$sort)), ]
cazy.gene.summary <- cazy.gene.summary[, 1:length(cazy.gene.summary) - 1]
cazy.gene.summary[1:50, 1:3]

cazy.gene.summary <- as.data.frame(t(cazy.gene.summary))
cazy.gene.list <- as.list(data.table(cazy.gene.summary))



# NCycDB =====================

ncyc.prop <- prop.table(ncyc.iss, margin = 1) * 100  # get proportions

ncyc.level1 <- aggregate(t(ncyc.prop) ~ level_1, data = ncyc.genes, FUN = sum)  
ncyc.level1$sort <- ncyc.level1[, 1]
ncyc.level1 <- data.frame(ncyc.level1[, -1], row.names = paste("l1:", ncyc.level1[, 1], sep = ""), check.names = FALSE)

ncyc.level2 <- aggregate(t(ncyc.prop) ~ level_1 + level_2, data = ncyc.genes, FUN = sum)
ncyc.level2 <- subset(ncyc.level2, ncyc.level2$level_2 != "unclassified")
ncyc.level2$sort <- paste(ncyc.level2[, 1], ncyc.level2[, 2], sep = ";")
ncyc.level2 <- data.frame(ncyc.level2[, -c(1:2)], row.names = paste("_l2:", ncyc.level2[, 2], sep = ""), check.names = FALSE)

ncyc.gene.summary <- rbind(ncyc.level1, ncyc.level2)
ncyc.gene.summary <- ncyc.gene.summary[order(as.character(ncyc.gene.summary$sort)), ]
ncyc.gene.summary <- ncyc.gene.summary[, 1:length(ncyc.gene.summary) - 1]
ncyc.gene.summary[1:50, 1:3]

ncyc.gene.summary <- as.data.frame(t(ncyc.gene.summary))
ncyc.gene.list <- as.list(data.table(ncyc.gene.summary))



# PCycDB =====================

pcyc.prop <- prop.table(pcyc.iss, margin = 1) * 100  # get proportions

pcyc.level1 <- aggregate(t(pcyc.prop) ~ level_1, data = pcyc.genes, FUN = sum)  
pcyc.level1$sort <- pcyc.level1[, 1]
pcyc.level1 <- data.frame(pcyc.level1[, -1], row.names = paste("l1:", pcyc.level1[, 1], sep = ""), check.names = FALSE)

pcyc.level2 <- aggregate(t(pcyc.prop) ~ level_1 + level_2, data = pcyc.genes, FUN = sum)
pcyc.level2 <- subset(pcyc.level2, pcyc.level2$level_2 != "unclassified")
pcyc.level2$sort <- paste(pcyc.level2[, 1], pcyc.level2[, 2], sep = ";")
pcyc.level2 <- data.frame(pcyc.level2[, -c(1:2)], row.names = paste("_l2:", pcyc.level2[, 2], sep = ""), check.names = FALSE)

pcyc.gene.summary <- rbind(pcyc.level1, pcyc.level2)
pcyc.gene.summary <- pcyc.gene.summary[order(as.character(pcyc.gene.summary$sort)), ]
pcyc.gene.summary <- pcyc.gene.summary[, 1:length(pcyc.gene.summary) - 1]
pcyc.gene.summary[1:50, 1:3]

pcyc.gene.summary <- as.data.frame(t(pcyc.gene.summary))
pcyc.gene.list <- as.list(data.table(pcyc.gene.summary))



# Permanova  ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy.gene.list.permanova <- mclapply(cazy.gene.list, function(x) adonis2(vegdist(x, method = "euclidean") ~ Irrigation * Treatment * Date, data = design, permutations = 9999), mc.cores = 20)

cazy.level.permanova.W <- do.call(rbind, lapply(cazy.gene.list.permanova, function(x) x[1, c(4, 5)]))
colnames(cazy.level.permanova.W) <- c("permanova.W.F", "permanova.W.P")
cazy.level.permanova.C <- do.call(rbind, lapply(cazy.gene.list.permanova, function(x) x[2, c(4, 5)]))
colnames(cazy.level.permanova.C) <- c("permanova.C.F", "permanova.C.P")
cazy.level.permanova.WxC <- do.call(rbind, lapply(cazy.gene.list.permanova, function(x) x[4, c(4, 5)]))
colnames(cazy.level.permanova.WxC) <- c("permanova.WxC.F", "permanova.WxC.P")

cazy.level.permanova <- cbind(cazy.level.permanova.W, cazy.level.permanova.C, cazy.level.permanova.WxC)
str(cazy.level.permanova)



# NCycDB =====================

ncyc.gene.list.permanova <- mclapply(ncyc.gene.list, function(x) adonis2(vegdist(x, method = "euclidean") ~ Irrigation * Treatment * Date, data = design, permutations = 9999), mc.cores = 20)

ncyc.level.permanova.W <- do.call(rbind, lapply(ncyc.gene.list.permanova, function(x) x[1, c(4, 5)]))
colnames(ncyc.level.permanova.W) <- c("permanova.W.F", "permanova.W.P")
ncyc.level.permanova.C <- do.call(rbind, lapply(ncyc.gene.list.permanova, function(x) x[2, c(4, 5)]))
colnames(ncyc.level.permanova.C) <- c("permanova.C.F", "permanova.C.P")
ncyc.level.permanova.WxC <- do.call(rbind, lapply(ncyc.gene.list.permanova, function(x) x[4, c(4, 5)]))
colnames(ncyc.level.permanova.WxC) <- c("permanova.WxC.F", "permanova.WxC.P")

ncyc.level.permanova <- cbind(ncyc.level.permanova.W, ncyc.level.permanova.C, ncyc.level.permanova.WxC)
str(ncyc.level.permanova)



# PCycDB =====================

pcyc.gene.list.permanova <- mclapply(pcyc.gene.list, function(x) adonis2(vegdist(x, method = "euclidean") ~ Irrigation * Treatment * Date, data = design, permutations = 9999), mc.cores = 20)

pcyc.level.permanova.W <- do.call(rbind, lapply(pcyc.gene.list.permanova, function(x) x[1, c(4, 5)]))
colnames(pcyc.level.permanova.W) <- c("permanova.W.F", "permanova.W.P")
pcyc.level.permanova.C <- do.call(rbind, lapply(pcyc.gene.list.permanova, function(x) x[2, c(4, 5)]))
colnames(pcyc.level.permanova.C) <- c("permanova.C.F", "permanova.C.P")
pcyc.level.permanova.WxC <- do.call(rbind, lapply(pcyc.gene.list.permanova, function(x) x[4, c(4, 5)]))
colnames(pcyc.level.permanova.WxC) <- c("permanova.WxC.F", "permanova.WxC.P")

pcyc.level.permanova <- cbind(pcyc.level.permanova.W, pcyc.level.permanova.C, pcyc.level.permanova.WxC)
str(pcyc.level.permanova)



# Permdisp ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy.list.permdisp.W <- mclapply(cazy.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Irrigation), permutations = how(nperm = 9999)), mc.cores = 20)
cazy.list.permdisp.C <- mclapply(cazy.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Treatment), permutations = how(nperm = 9999)), mc.cores = 20)
cazy.list.permdisp.WxC <- mclapply(cazy.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$IxT), permutations = how(nperm = 9999)), mc.cores = 20)

cazy.permdisp.W <- do.call(rbind, lapply(cazy.list.permdisp.W, function(x) x$tab[1, c(4, 6)]))
colnames(cazy.permdisp.W) <- c("permdisp.W.F", "permdisp.W.P")
cazy.permdisp.C <- do.call(rbind, lapply(cazy.list.permdisp.C, function(x) x$tab[1, c(4, 6)]))
colnames(cazy.permdisp.C) <- c("permdisp.C.F", "permdisp.C.P")
cazy.permdisp.WxC <- do.call(rbind, lapply(cazy.list.permdisp.WxC, function(x) x$tab[1, c(4, 6)]))
colnames(cazy.permdisp.WxC) <- c("permdisp.WxC.F", "permdisp.WxC.P")

cazy.permdisp <- cbind(cazy.permdisp.W, cazy.permdisp.C, cazy.permdisp.WxC)
str(cazy.permdisp)



# NCycDB =====================

ncyc.list.permdisp.W <- mclapply(ncyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Irrigation), permutations = how(nperm = 9999)), mc.cores = 20)
ncyc.list.permdisp.C <- mclapply(ncyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Treatment), permutations = how(nperm = 9999)), mc.cores = 20)
ncyc.list.permdisp.WxC <- mclapply(ncyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$IxT), permutations = how(nperm = 9999)), mc.cores = 20)

ncyc.permdisp.W <- do.call(rbind, lapply(ncyc.list.permdisp.W, function(x) x$tab[1, c(4, 6)]))
colnames(ncyc.permdisp.W) <- c("permdisp.W.F", "permdisp.W.P")
ncyc.permdisp.C <- do.call(rbind, lapply(ncyc.list.permdisp.C, function(x) x$tab[1, c(4, 6)]))
colnames(ncyc.permdisp.C) <- c("permdisp.C.F", "permdisp.C.P")
ncyc.permdisp.WxC <- do.call(rbind, lapply(ncyc.list.permdisp.WxC, function(x) x$tab[1, c(4, 6)]))
colnames(ncyc.permdisp.WxC) <- c("permdisp.WxC.F", "permdisp.WxC.P")

ncyc.permdisp <- cbind(ncyc.permdisp.W, ncyc.permdisp.C, ncyc.permdisp.WxC)
str(ncyc.permdisp)



# PCycDB =====================

pcyc.list.permdisp.W <- mclapply(pcyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Irrigation), permutations = how(nperm = 9999)), mc.cores = 20)
pcyc.list.permdisp.C <- mclapply(pcyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Treatment), permutations = how(nperm = 9999)), mc.cores = 20)
pcyc.list.permdisp.WxC <- mclapply(pcyc.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$IxT), permutations = how(nperm = 9999)), mc.cores = 20)

pcyc.permdisp.W <- do.call(rbind, lapply(pcyc.list.permdisp.W, function(x) x$tab[1, c(4, 6)]))
colnames(pcyc.permdisp.W) <- c("permdisp.W.F", "permdisp.W.P")
pcyc.permdisp.C <- do.call(rbind, lapply(pcyc.list.permdisp.C, function(x) x$tab[1, c(4, 6)]))
colnames(pcyc.permdisp.C) <- c("permdisp.C.F", "permdisp.C.P")
pcyc.permdisp.WxC <- do.call(rbind, lapply(pcyc.list.permdisp.WxC, function(x) x$tab[1, c(4, 6)]))
colnames(pcyc.permdisp.WxC) <- c("permdisp.WxC.F", "permdisp.WxC.P")

pcyc.permdisp <- cbind(pcyc.permdisp.W, pcyc.permdisp.C, pcyc.permdisp.WxC)
str(pcyc.permdisp)



# Scale - relative abundances ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy.scale <- scale(cazy.gene.summary)
round(colMeans(cazy.scale[, 1:5]), digits = 5)
colSds(cazy.scale[, 1:5])

# means, stdev, sterr by water regime

cazy.scale.W.mean <- aggregate(cazy.scale, list(design$Irrigation), mean)
cazy.scale.W.mean <- t(as.matrix(data.frame(cazy.scale.W.mean[, -1], row.names = paste(cazy.scale.W.mean[, 1], "W.mean", sep = "."), check.names = FALSE)))

cazy.scale.W.sd <- aggregate(cazy.scale, list(design$Irrigation), FUN = sd)
cazy.scale.W.sd <- t(as.matrix(data.frame(cazy.scale.W.sd[, -1], row.names = paste(cazy.scale.W.sd[, 1], "W.sd", sep = "."), check.names = FALSE)))

cazy.scale.W.se <- aggregate(cazy.scale, list(design$Irrigation), FUN = function(x) sd(x)/length(x))
cazy.scale.W.se <- t(as.matrix(data.frame(cazy.scale.W.se[, -1], row.names = paste(cazy.scale.W.se[, 1], "W.se", sep = "."), check.names = FALSE)))

cazy.scale.W <- cbind(cazy.scale.W.mean, cazy.scale.W.sd, cazy.scale.W.se)
str(cazy.scale.W)

# means, stdev, sterr by cropping system

cazy.scale.C.mean <- aggregate(cazy.scale, list(design$Treatment), mean)
cazy.scale.C.mean <- t(as.matrix(data.frame(cazy.scale.C.mean[, -1], row.names = paste(cazy.scale.C.mean[, 1], "C.mean", sep = "."), check.names = FALSE)))

cazy.scale.C.sd <- aggregate(cazy.scale, list(design$Treatment), FUN = sd)
cazy.scale.C.sd <- t(as.matrix(data.frame(cazy.scale.C.sd[, -1], row.names = paste(cazy.scale.C.sd[, 1], "C.sd", sep = "."), check.names = FALSE)))

cazy.scale.C.se <- aggregate(cazy.scale, list(design$Treatment), FUN = function(x) sd(x)/length(x))
cazy.scale.C.se <- t(as.matrix(data.frame(cazy.scale.C.se[, -1], row.names = paste(cazy.scale.C.se[, 1], "C.se", sep = "."), check.names = FALSE)))

cazy.scale.C <- cbind(cazy.scale.C.mean, cazy.scale.C.sd, cazy.scale.C.se)
str(cazy.scale.C)

# means, stdev, sterr by interaction

cazy.scale.WxC.mean <- aggregate(cazy.scale, list(design$IxT), mean)
cazy.scale.WxC.mean <- t(as.matrix(data.frame(cazy.scale.WxC.mean[, -1], row.names = paste(cazy.scale.WxC.mean[, 1], "WxC.mean", sep = "."), check.names = FALSE)))

cazy.scale.WxC.sd <- aggregate(cazy.scale, list(design$IxT), FUN = sd)
cazy.scale.WxC.sd <- t(as.matrix(data.frame(cazy.scale.WxC.sd[, -1], row.names = paste(cazy.scale.WxC.sd[, 1], "WxC.sd", sep = "."), check.names = FALSE)))

cazy.scale.WxC.se <- aggregate(cazy.scale, list(design$IxT), FUN = function(x) sd(x)/length(x))
cazy.scale.WxC.se <- t(as.matrix(data.frame(cazy.scale.WxC.se[, -1], row.names = paste(cazy.scale.WxC.se[, 1], "WxC.se", sep = "."), check.names = FALSE)))

cazy.scale.WxC <- cbind(cazy.scale.WxC.mean, cazy.scale.WxC.sd, cazy.scale.WxC.se)
str(cazy.scale.WxC)

cazy.scale <- cbind(cazy.scale.W, cazy.scale.C, cazy.scale.WxC)
str(cazy.scale)



# NCycDB =====================

ncyc.scale <- scale(ncyc.gene.summary)
round(colMeans(ncyc.scale[, 1:5]), digits = 5)
colSds(ncyc.scale[, 1:5])

# means, stdev, sterr by water regime

ncyc.scale.W.mean <- aggregate(ncyc.scale, list(design$Irrigation), mean)
ncyc.scale.W.mean <- t(as.matrix(data.frame(ncyc.scale.W.mean[, -1], row.names = paste(ncyc.scale.W.mean[, 1], "W.mean", sep = "."), check.names = FALSE)))

ncyc.scale.W.sd <- aggregate(ncyc.scale, list(design$Irrigation), FUN = sd)
ncyc.scale.W.sd <- t(as.matrix(data.frame(ncyc.scale.W.sd[, -1], row.names = paste(ncyc.scale.W.sd[, 1], "W.sd", sep = "."), check.names = FALSE)))

ncyc.scale.W.se <- aggregate(ncyc.scale, list(design$Irrigation), FUN = function(x) sd(x)/length(x))
ncyc.scale.W.se <- t(as.matrix(data.frame(ncyc.scale.W.se[, -1], row.names = paste(ncyc.scale.W.se[, 1], "W.se", sep = "."), check.names = FALSE)))

ncyc.scale.W <- cbind(ncyc.scale.W.mean, ncyc.scale.W.sd, ncyc.scale.W.se)
str(ncyc.scale.W)

# means, stdev, sterr by cropping system

ncyc.scale.C.mean <- aggregate(ncyc.scale, list(design$Treatment), mean)
ncyc.scale.C.mean <- t(as.matrix(data.frame(ncyc.scale.C.mean[, -1], row.names = paste(ncyc.scale.C.mean[, 1], "C.mean", sep = "."), check.names = FALSE)))

ncyc.scale.C.sd <- aggregate(ncyc.scale, list(design$Treatment), FUN = sd)
ncyc.scale.C.sd <- t(as.matrix(data.frame(ncyc.scale.C.sd[, -1], row.names = paste(ncyc.scale.C.sd[, 1], "C.sd", sep = "."), check.names = FALSE)))

ncyc.scale.C.se <- aggregate(ncyc.scale, list(design$Treatment), FUN = function(x) sd(x)/length(x))
ncyc.scale.C.se <- t(as.matrix(data.frame(ncyc.scale.C.se[, -1], row.names = paste(ncyc.scale.C.se[, 1], "C.se", sep = "."), check.names = FALSE)))

ncyc.scale.C <- cbind(ncyc.scale.C.mean, ncyc.scale.C.sd, ncyc.scale.C.se)
str(ncyc.scale.C)

# means, stdev, sterr by interaction

ncyc.scale.WxC.mean <- aggregate(ncyc.scale, list(design$IxT), mean)
ncyc.scale.WxC.mean <- t(as.matrix(data.frame(ncyc.scale.WxC.mean[, -1], row.names = paste(ncyc.scale.WxC.mean[, 1], "WxC.mean", sep = "."), check.names = FALSE)))

ncyc.scale.WxC.sd <- aggregate(ncyc.scale, list(design$IxT), FUN = sd)
ncyc.scale.WxC.sd <- t(as.matrix(data.frame(ncyc.scale.WxC.sd[, -1], row.names = paste(ncyc.scale.WxC.sd[, 1], "WxC.sd", sep = "."), check.names = FALSE)))

ncyc.scale.WxC.se <- aggregate(ncyc.scale, list(design$IxT), FUN = function(x) sd(x)/length(x))
ncyc.scale.WxC.se <- t(as.matrix(data.frame(ncyc.scale.WxC.se[, -1], row.names = paste(ncyc.scale.WxC.se[, 1], "WxC.se", sep = "."), check.names = FALSE)))

ncyc.scale.WxC <- cbind(ncyc.scale.WxC.mean, ncyc.scale.WxC.sd, ncyc.scale.WxC.se)
str(ncyc.scale.WxC)

ncyc.scale <- cbind(ncyc.scale.W, ncyc.scale.C, ncyc.scale.WxC)
str(ncyc.scale)



# PCycDB =====================

pcyc.scale <- scale(pcyc.gene.summary)
round(colMeans(pcyc.scale[, 1:5]), digits = 5)
colSds(pcyc.scale[, 1:5])

# means, stdev, sterr by water regime

pcyc.scale.W.mean <- aggregate(pcyc.scale, list(design$Irrigation), mean)
pcyc.scale.W.mean <- t(as.matrix(data.frame(pcyc.scale.W.mean[, -1], row.names = paste(pcyc.scale.W.mean[, 1], "W.mean", sep = "."), check.names = FALSE)))

pcyc.scale.W.sd <- aggregate(pcyc.scale, list(design$Irrigation), FUN = sd)
pcyc.scale.W.sd <- t(as.matrix(data.frame(pcyc.scale.W.sd[, -1], row.names = paste(pcyc.scale.W.sd[, 1], "W.sd", sep = "."), check.names = FALSE)))

pcyc.scale.W.se <- aggregate(pcyc.scale, list(design$Irrigation), FUN = function(x) sd(x)/length(x))
pcyc.scale.W.se <- t(as.matrix(data.frame(pcyc.scale.W.se[, -1], row.names = paste(pcyc.scale.W.se[, 1], "W.se", sep = "."), check.names = FALSE)))

pcyc.scale.W <- cbind(pcyc.scale.W.mean, pcyc.scale.W.sd, pcyc.scale.W.se)
str(pcyc.scale.W)

# means, stdev, sterr by cropping system

pcyc.scale.C.mean <- aggregate(pcyc.scale, list(design$Treatment), mean)
pcyc.scale.C.mean <- t(as.matrix(data.frame(pcyc.scale.C.mean[, -1], row.names = paste(pcyc.scale.C.mean[, 1], "C.mean", sep = "."), check.names = FALSE)))

pcyc.scale.C.sd <- aggregate(pcyc.scale, list(design$Treatment), FUN = sd)
pcyc.scale.C.sd <- t(as.matrix(data.frame(pcyc.scale.C.sd[, -1], row.names = paste(pcyc.scale.C.sd[, 1], "C.sd", sep = "."), check.names = FALSE)))

pcyc.scale.C.se <- aggregate(pcyc.scale, list(design$Treatment), FUN = function(x) sd(x)/length(x))
pcyc.scale.C.se <- t(as.matrix(data.frame(pcyc.scale.C.se[, -1], row.names = paste(pcyc.scale.C.se[, 1], "C.se", sep = "."), check.names = FALSE)))

pcyc.scale.C <- cbind(pcyc.scale.C.mean, pcyc.scale.C.sd, pcyc.scale.C.se)
str(pcyc.scale.C)

# means, stdev, sterr by interaction

pcyc.scale.WxC.mean <- aggregate(pcyc.scale, list(design$IxT), mean)
pcyc.scale.WxC.mean <- t(as.matrix(data.frame(pcyc.scale.WxC.mean[, -1], row.names = paste(pcyc.scale.WxC.mean[, 1], "WxC.mean", sep = "."), check.names = FALSE)))

pcyc.scale.WxC.sd <- aggregate(pcyc.scale, list(design$IxT), FUN = sd)
pcyc.scale.WxC.sd <- t(as.matrix(data.frame(pcyc.scale.WxC.sd[, -1], row.names = paste(pcyc.scale.WxC.sd[, 1], "WxC.sd", sep = "."), check.names = FALSE)))

pcyc.scale.WxC.se <- aggregate(pcyc.scale, list(design$IxT), FUN = function(x) sd(x)/length(x))
pcyc.scale.WxC.se <- t(as.matrix(data.frame(pcyc.scale.WxC.se[, -1], row.names = paste(pcyc.scale.WxC.se[, 1], "WxC.se", sep = "."), check.names = FALSE)))

pcyc.scale.WxC <- cbind(pcyc.scale.WxC.mean, pcyc.scale.WxC.sd, pcyc.scale.WxC.se)
str(pcyc.scale.WxC)

pcyc.scale <- cbind(pcyc.scale.W, pcyc.scale.C, pcyc.scale.WxC)
str(pcyc.scale)



# Q-value - correction for multiple testing ==============================================================================================================================================================================================================================================


# CAZyme =====================

# check p-value distribution

par(mfrow = c(1, 3))
hist(cazy.permanova$permanova.W.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Water regime", xlim = c(0, 1))
hist(cazy.permanova$permanova.C.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Cropping system", xlim = c(0, 1))
hist(cazy.permanova$permanova.WxC.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Treatment x System", xlim = c(0, 1))

# apply FDR

FDR <- qvalue(cazy.permanova$permanova.W.P)  # check fit with 'hist(FDR)'
cazy.permanova$permanova.W.q <- FDR$qvalues
FDR <- qvalue(cazy.permanova$permanova.C.P)  # check fit with 'hist(FDR)'
cazy.permanova$permanova.C.q <- FDR$qvalues
FDR <- qvalue(cazy.permanova$permanova.WxC.P, pi0 = 0.6)  # check fit with 'hist(FDR)'
cazy.permanova$permanova.WxC.q <- FDR$qvalues



# NCycDB =====================

# check p-value distribution

par(mfrow = c(1, 3))
hist(ncyc.permanova$permanova.W.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Water regime", xlim = c(0, 1))
hist(ncyc.permanova$permanova.C.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Cropping system", xlim = c(0, 1))
hist(ncyc.permanova$permanova.WxC.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Treatment x System", xlim = c(0, 1))

# apply FDR

FDR <- qvalue(ncyc.permanova$permanova.W.P)  # check fit with 'hist(FDR)'
ncyc.permanova$permanova.W.q <- FDR$qvalues
FDR <- qvalue(ncyc.permanova$permanova.C.P)  # check fit with 'hist(FDR)'
ncyc.permanova$permanova.C.q <- FDR$qvalues
FDR <- qvalue(ncyc.permanova$permanova.WxC.P)  # check fit with 'hist(FDR)'
ncyc.permanova$permanova.WxC.q <- FDR$qvalues



# PCycDB =====================

# check p-value distribution

par(mfrow = c(1, 3))
hist(pcyc.permanova$permanova.W.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Water regime", xlim = c(0, 1))
hist(pcyc.permanova$permanova.C.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Cropping system", xlim = c(0, 1))
hist(pcyc.permanova$permanova.WxC.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Treatment x System", xlim = c(0, 1))

# apply FDR

FDR <- qvalue(pcyc.permanova$permanova.W.P)  # check fit with 'hist(FDR)'
pcyc.permanova$permanova.W.q <- FDR$qvalues
FDR <- qvalue(pcyc.permanova$permanova.C.P, pi0 = 0.3)  # check fit with 'hist(FDR)'
pcyc.permanova$permanova.C.q <- FDR$qvalues
FDR <- qvalue(pcyc.permanova$permanova.WxC.P)  # check fit with 'hist(FDR)'
pcyc.permanova$permanova.WxC.q <- FDR$qvalues



# Combine files ==============================================================================================================================================================================================================================================

# CAZyme =====================

cazy.gene.level.analysis <- cbind(cazy.scale, cazy.permanova, cazy.permdisp)
cazy.gene.level.analysis$count <- colSums(cazy.gene.summary)
cazy.gene.level.analysis$relabund <- 100/sum(colSums(cazy.gene.summary)) * colSums(cazy.gene.summary)
str(cazy.gene.level.analysis)



# NCycDB =====================

ncyc.gene.level.analysis <- cbind(ncyc.scale, ncyc.permanova, ncyc.permdisp)
ncyc.gene.level.analysis$count <- colSums(ncyc.gene.summary)
ncyc.gene.level.analysis$relabund <- 100/sum(colSums(ncyc.gene.summary)) * colSums(ncyc.gene.summary)
str(ncyc.gene.level.analysis)



# PCycDB =====================

pcyc.gene.level.analysis <- cbind(pcyc.scale, pcyc.permanova, pcyc.permdisp)
pcyc.gene.level.analysis$count <- colSums(pcyc.gene.summary)
pcyc.gene.level.analysis$relabund <- 100/sum(colSums(pcyc.gene.summary)) * colSums(pcyc.gene.summary)
str(pcyc.gene.level.analysis)



# Check results ==============================================================================================================================================================================================================================================

# CAZyme ===================== 

level1.cazy <- cazy.gene.level.analysis %>% filter(grepl("^l[0-9]:", rownames(cazy.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) 

level2.cazy <- cazy.gene.level.analysis %>% filter(grepl("^_l[0-9]:", rownames(cazy.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) 



# NCycDB ===================== 

level1.ncyc <- ncyc.gene.level.analysis %>% filter(grepl("^l[0-9]:", rownames(ncyc.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) 

level2.ncyc <- ncyc.gene.level.analysis %>% filter(grepl("^_l[0-9]:", rownames(ncyc.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) 



# PCycDB ===================== 
 
level1.pcyc <- pcyc.gene.level.analysis %>% filter(grepl("^l[0-9]:", rownames(pcyc.gene.level.analysis))) %>%
   filter(permanova.W.q < 0.05) 
 
level2.pcyc <- pcyc.gene.level.analysis %>% filter(grepl("^_l[0-9]:", rownames(pcyc.gene.level.analysis))) %>%
   filter(permanova.W.q < 0.05) 



# Visualization Level 1 ==============================================================================================================================================================================================================================================

# rename 

rownames(level1.cazy) <- sub("l1:*", "", rownames(level1.cazy))
rownames(level1.ncyc) <- sub("l1:*", "", rownames(level1.ncyc))
rownames(level1.pcyc) <- sub("l1:*", "", rownames(level1.pcyc))


# combine files

all.levels.1 <- rbind(level1.cazy, level1.ncyc, level1.pcyc)

all.genes.1 <- rbind(cazy.genes, ncyc.genes, pcyc.genes)
row.names(all.genes.1) <- NULL
all.genes.1 <- all.genes.1[,-2]
all.genes.1 <- all.genes.1[!duplicated(all.genes.1), ]
row.names(all.genes.1) <- NULL
all.genes.1 <- tibble::column_to_rownames(all.genes.1, "level_1")


# get gene levels which were significant

sig.genes.1 <- all.genes.1 %>%
  filter(rownames(.) %in% rownames(all.levels.1))
rownames(sig.genes.1) <- gsub(" ", "_", rownames(sig.genes.1))


# Create files for iTol tree 

# prepare clusters

taxdis.iToL.1 <- taxa2dist(sig.genes.1, varstep = TRUE)
hc.1 <- hclust(taxdis.iToL.1)
hcp.1 <- as.phylo(hc.1)

write.tree(hcp.1, file = "iTol.nc.l1.newick", append = FALSE, digits = 10, tree.names = FALSE)

# get enrichment/depletion under drought

rownames(all.levels.1) <- gsub(" ", "_", rownames(all.levels.1))

iTol_relab.1 <- all.levels.1 %>%
  dplyr::select(control.W.mean, rainout.W.mean) %>%
  mutate(relab.pos = rainout.W.mean) %>%
  mutate(relab.neg = control.W.mean * -1) %>%
  dplyr::select(relab.pos, relab.neg)

iTol_relab.1$relab.neg[iTol_relab.1$relab.neg > 0] <- 0
iTol_relab.1$relab.pos[iTol_relab.1$relab.pos < 0] <- 0

write.table(iTol_relab.1, file = "barplots.relab.nc.l1.txt", sep = ",") # write txt file

# get relative abundance

iTol_ab.1 <- all.levels.1 %>% dplyr::select(relabund) 

write.table(iTol_ab.1, file = "barplots.ab.nc.l1.txt", sep = ",") # write txt file


 
# Visualization Level 2 ==============================================================================================================================================================================================================================================

# rename 

rownames(level2.cazy) <- sub("_l2:*", "", rownames(level2.cazy))
rownames(level2.ncyc) <- sub("_l2:*", "", rownames(level2.ncyc))
rownames(level2.pcyc) <- sub("_l2:*", "", rownames(level2.pcyc))


# combine files

all.levels.2 <- rbind(level2.cazy, level2.ncyc, level2.pcyc)

all.genes.2 <- rbind(cazy.genes, ncyc.genes, pcyc.genes)
row.names(all.genes.2) <- NULL
all.genes.2 <- tibble::column_to_rownames(all.genes.2, "level_2")


# get gene levels which were significant

sig.genes.2 <- all.genes.2 %>%
  filter(rownames(.) %in% rownames(all.levels.2))
sig.genes$level_1 <- gsub(" ", "_", sig.genes$level_1)


# Create files for iTol tree 

# prepare clusters

taxdis.iToL.2 <- taxa2dist(sig.genes.2, varstep = TRUE)
hc.2 <- hclust(taxdis.iToL.2)
hcp.2 <- as.phylo(hc.2)

write.tree(hcp.2, file = "iTol.nc.l2.newick", append = FALSE, digits = 10, tree.names = FALSE)

# get enrichment/depletion under drought

iTol_relab.2 <- all.levels.2 %>%
  dplyr::select(control.W.mean, rainout.W.mean) %>%
  mutate(relab.pos = rainout.W.mean) %>%
  mutate(relab.neg = control.W.mean * -1) %>%
  dplyr::select(relab.pos, relab.neg)

iTol_relab.2$relab.neg[iTol_relab$relab.neg > 0] <- 0
iTol_relab.2$relab.pos[iTol_relab$relab.pos < 0] <- 0

write.table(iTol_relab.2, file = "barplots.relab.nc.l2.txt", sep = ",") # write txt file

# get relative abundance

iTol_ab.2 <- all.levels.2 %>% dplyr::select(relabund) 

write.table(iTol_ab.2, file = "barplots.ab.nc.l2.txt", sep = ",") # write txt file


