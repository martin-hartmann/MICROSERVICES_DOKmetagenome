
# R Script for Analysis of the metagenome samples of the DOK rainout-trial 2022 - indicator analysis
# Elena Kost
# 16.04.25

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, plyr, Nonpareil, eulerr, 
               ggvenn, VennDiagram, UpSetR)
library(pacman)



# Import data SEED ==============================================================================================================================================================================================================================================

seed <- read.delim("SEED_RPKmaxORF.tsv",sep="\t")
colnames(seed) <- sub("S\\.", "S-", colnames(seed))
str(seed)
seed[,-c(1:5)][is.na(seed[,-c(1:5)])] <- 0
seed[,c(1:5)][is.na(seed[,c(1:5)])] <- "unclassified"

seed.genes <- seed[,c(1:5)]

seed.genes <- seed.genes %>% relocate(accession) %>% 
  column_to_rownames("accession") 



# Import Design ==============================================================================================================================================================================================================================================

design <- read.table("design_metagenomic.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE,
                     colClasses = "factor")
str(design)

design$IxT <- with(design, interaction(Irrigation, Treatment, lex.order = TRUE, drop = TRUE))



# Normalization data ==============================================================================================================================================================================================================================================

iters <- 100
seed.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(seed[,-c(1:5)])),  min(rowSums(as.data.frame(t(seed[,-c(1:5)]))))), mc.cores = detectCores())
seed.array <- laply(seed.iters, as.matrix)
seed.iss <- apply(seed.array, 2:3, mean)



# aggregate data SEED ==============================================================================================================================================================================================================================================

seed.prop <- prop.table(seed.iss, margin = 1) * 100  # get proportions

seed.level1 <- aggregate(t(seed.prop) ~ level_1, data = seed.genes, FUN = sum)  
seed.level1$sort <- seed.level1[, 1]
seed.level1 <- data.frame(seed.level1[, -1], row.names = paste("l1:", seed.level1[, 1], sep = ""), check.names = FALSE)

seed.level2 <- aggregate(t(seed.prop) ~ level_1 + level_2, data = seed.genes, FUN = sum)
seed.level2 <- subset(seed.level2, seed.level2$level_2 != "unclassified")
seed.level2$sort <- paste(seed.level2[, 1], seed.level2[, 2], sep = ";")
seed.level2 <- data.frame(seed.level2[, -c(1:2)], row.names = paste("_l2:", seed.level2[, 2], sep = ""), check.names = FALSE)

seed.level3 <- aggregate(t(seed.prop) ~ level_1 + level_2 + level_3, data = seed.genes, FUN = sum)
seed.level3 <- subset(seed.level3, seed.level3$level_3 != "unclassified")
seed.level3$sort <- paste(seed.level3[, 1], seed.level3[, 2], seed.level3[, 3], sep = ";")
seed.level3 <- data.frame(seed.level3[, -c(1:3)], row.names = paste("__l3:", seed.level3[, 3], sep = ""), check.names = FALSE)

seed.level4 <- aggregate(t(seed.prop) ~ level_1 + level_2 + level_3 + level_4, data = seed.genes, FUN = sum)
seed.level4 <- subset(seed.level4, seed.level4$level_4 != "unclassified")
seed.level4$sort <- paste(seed.level4[, 1], seed.level4[, 2], seed.level4[, 3], seed.level4[, 4], sep = ";")
seed.level4 <- data.frame(seed.level4[, -c(1:4)], row.names = paste("___l4:", seed.level4[, 4], sep = ""), check.names = FALSE)

seed.gene.summary <- rbind(seed.level1, seed.level2, seed.level3, seed.level4)
seed.gene.summary <- seed.gene.summary[order(as.character(seed.gene.summary$sort)), ]
seed.gene.summary <- seed.gene.summary[, 1:length(seed.gene.summary) - 1]
seed.gene.summary[1:50, 1:3]

seed.gene.summary <- t(seed.gene.summary)



# Indicator analysis SEED ==============================================================================================================================================================================================================================================

# Correlation value ===========

indR.seed.cor <- multipatt(seed.gene.summary, design$IxT, func = "r.g", control = how(nperm = 9999))
indR.seed.cor <- cbind(indR.seed.cor$str, indR.seed.cor$sign)

hist(indR.seed.cor$p.value, breaks = seq(0, 1, 0.05)) # adjust for multiple testing
FDR <- qvalue(indR.seed.cor$p.value)
hist(FDR)

indR.seed.cor$q.value <- FDR$qvalues
indR.seed.cor$count <- colSums(seed.gene.summary)

rownames(subset(indR.seed.cor, q.value < 0.05 & index == 4)) # get indicator gene functions for enriched under drought in BIODYN
rownames(subset(indR.seed.cor, q.value < 0.05 & index == 5)) # get indicator gene functions for enriched under drought in CONFYM
rownames(subset(indR.seed.cor, q.value < 0.05 & index == 6)) # get indicator gene functions for enriched under drought in CONMIN

rownames(subset(indR.seed.cor, q.value < 0.05 & index == 19)) # get indicator gene functions for enriched under drought in BIODYN and CONFYM
rownames(subset(indR.seed.cor, q.value < 0.05 & index == 20)) # get indicator gene functions for enriched under drought in BIODYN and CONMIN
rownames(subset(indR.seed.cor, q.value < 0.05 & index == 21)) # get indicator gene functions for enriched under drought in CONFYM and CONMIN

rownames(subset(indR.seed.cor, q.value < 0.05 & index == 41)) # get indicator gene functions for enriched under drought in all systems combined



# Visualization using upset plot ===========

indR.seed.final <- indR.seed.cor %>% filter(q.value < 0.05)
indR.seed.final <- indR.seed.final %>% filter(index == 4 | index == 5 | index == 6 | index == 19 | index == 20 | index == 21 | index == 41)  

mat.seed.final <- indR.seed.final[, 66:68]

mat.seed.final <- mat.seed.final %>% 
  rename("BIODYN drought-induced" = "s.rainout.d",
         "CONFYM drought-induced" = "s.rainout.k",
         "CONMIN drought-induced" = "s.rainout.m")

upset.plot <- upset(mat.seed.final, sets.bar.color = "grey23", 
                    mainbar.y.label = "# gene functions",
                    sets.x.label = "Total # gene functions",
                    main.bar.color = c("#E69F00", "#FF618C", "#009E73", "#F58359", "#A38D81", "#9BA14A", "#C09262"))

upset.plot



