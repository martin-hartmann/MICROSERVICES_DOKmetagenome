
# R Script for Analysis of the metagenome samples of the DOK rainout-trial 2022
# Gene level analysis of SEED
# Elena Kost
# 02.09.24

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, stats, ape)
library(pacman)



# Import data ==============================================================================================================================================================================================================================================

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



# Aggregate data  ==============================================================================================================================================================================================================================================

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

seed.gene.summary <- as.data.frame(t(seed.gene.summary))
seed.gene.list <- as.list(data.table(seed.gene.summary))



# Permanova  ==============================================================================================================================================================================================================================================

seed.gene.list.permanova <- mclapply(seed.gene.list, function(x) adonis2(vegdist(x, method = "euclidean") ~ Irrigation * Treatment * Date, data = design, permutations = 9999), mc.cores = 20)

seed.level.permanova.W <- do.call(rbind, lapply(seed.gene.list.permanova, function(x) x[1, c(4, 5)]))
colnames(seed.level.permanova.W) <- c("permanova.W.F", "permanova.W.P")
seed.level.permanova.C <- do.call(rbind, lapply(seed.gene.list.permanova, function(x) x[2, c(4, 5)]))
colnames(seed.level.permanova.C) <- c("permanova.C.F", "permanova.C.P")
seed.level.permanova.WxC <- do.call(rbind, lapply(seed.gene.list.permanova, function(x) x[4, c(4, 5)]))
colnames(seed.level.permanova.WxC) <- c("permanova.WxC.F", "permanova.WxC.P")

seed.level.permanova <- cbind(seed.level.permanova.W, seed.level.permanova.C, seed.level.permanova.WxC)
str(seed.level.permanova)



# Permdisp ==============================================================================================================================================================================================================================================

seed.list.permdisp.W <- mclapply(seed.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Irrigation), permutations = how(nperm = 9999)), mc.cores = 20)
seed.list.permdisp.C <- mclapply(seed.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$Treatment), permutations = how(nperm = 9999)), mc.cores = 20)
seed.list.permdisp.WxC <- mclapply(seed.gene.list, function(x) permutest(betadisper(vegdist(x, method = "euclidean"), design$IxT), permutations = how(nperm = 9999)), mc.cores = 20)

# Extract info
seed.permdisp.W <- do.call(rbind, lapply(seed.list.permdisp.W, function(x) x$tab[1, c(4, 6)]))
colnames(seed.permdisp.W) <- c("permdisp.W.F", "permdisp.W.P")
seed.permdisp.C <- do.call(rbind, lapply(seed.list.permdisp.C, function(x) x$tab[1, c(4, 6)]))
colnames(seed.permdisp.C) <- c("permdisp.C.F", "permdisp.C.P")
seed.permdisp.WxC <- do.call(rbind, lapply(seed.list.permdisp.WxC, function(x) x$tab[1, c(4, 6)]))
colnames(seed.permdisp.WxC) <- c("permdisp.WxC.F", "permdisp.WxC.P")

seed.permdisp <- cbind(seed.permdisp.W, seed.permdisp.C, seed.permdisp.WxC)
str(seed.permdisp)



# Scale - relative abundances ==============================================================================================================================================================================================================================================

seed.scale <- scale(seed.gene.summary)
round(colMeans(seed.scale[, 1:5]), digits = 5)
colSds(seed.scale[, 1:5])

# means, stdev, sterr by water regime

seed.scale.W.mean <- aggregate(seed.scale, list(design$Irrigation), mean)
seed.scale.W.mean <- t(as.matrix(data.frame(seed.scale.W.mean[, -1], row.names = paste(seed.scale.W.mean[, 1], "W.mean", sep = "."), check.names = FALSE)))

seed.scale.W.sd <- aggregate(seed.scale, list(design$Irrigation), FUN = sd)
seed.scale.W.sd <- t(as.matrix(data.frame(seed.scale.W.sd[, -1], row.names = paste(seed.scale.W.sd[, 1], "W.sd", sep = "."), check.names = FALSE)))

seed.scale.W.se <- aggregate(seed.scale, list(design$Irrigation), FUN = function(x) sd(x)/length(x))
seed.scale.W.se <- t(as.matrix(data.frame(seed.scale.W.se[, -1], row.names = paste(seed.scale.W.se[, 1], "W.se", sep = "."), check.names = FALSE)))

seed.scale.W <- cbind(seed.scale.W.mean, seed.scale.W.sd, seed.scale.W.se)
str(seed.scale.W)

# means, stdev, sterr by cropping system

seed.scale.C.mean <- aggregate(seed.scale, list(design$Treatment), mean)
seed.scale.C.mean <- t(as.matrix(data.frame(seed.scale.C.mean[, -1], row.names = paste(seed.scale.C.mean[, 1], "C.mean", sep = "."), check.names = FALSE)))

seed.scale.C.sd <- aggregate(seed.scale, list(design$Treatment), FUN = sd)
seed.scale.C.sd <- t(as.matrix(data.frame(seed.scale.C.sd[, -1], row.names = paste(seed.scale.C.sd[, 1], "C.sd", sep = "."), check.names = FALSE)))

seed.scale.C.se <- aggregate(seed.scale, list(design$Treatment), FUN = function(x) sd(x)/length(x))
seed.scale.C.se <- t(as.matrix(data.frame(seed.scale.C.se[, -1], row.names = paste(seed.scale.C.se[, 1], "C.se", sep = "."), check.names = FALSE)))

seed.scale.C <- cbind(seed.scale.C.mean, seed.scale.C.sd, seed.scale.C.se)
str(seed.scale.C)

# means, stdev, sterr by interaction

seed.scale.WxC.mean <- aggregate(seed.scale, list(design$IxT), mean)
seed.scale.WxC.mean <- t(as.matrix(data.frame(seed.scale.WxC.mean[, -1], row.names = paste(seed.scale.WxC.mean[, 1], "WxC.mean", sep = "."), check.names = FALSE)))

seed.scale.WxC.sd <- aggregate(seed.scale, list(design$IxT), FUN = sd)
seed.scale.WxC.sd <- t(as.matrix(data.frame(seed.scale.WxC.sd[, -1], row.names = paste(seed.scale.WxC.sd[, 1], "WxC.sd", sep = "."), check.names = FALSE)))

seed.scale.WxC.se <- aggregate(seed.scale, list(design$IxT), FUN = function(x) sd(x)/length(x))
seed.scale.WxC.se <- t(as.matrix(data.frame(seed.scale.WxC.se[, -1], row.names = paste(seed.scale.WxC.se[, 1], "WxC.se", sep = "."), check.names = FALSE)))

seed.scale.WxC <- cbind(seed.scale.WxC.mean, seed.scale.WxC.sd, seed.scale.WxC.se)
str(seed.scale.WxC)

seed.scale <- cbind(seed.scale.W, seed.scale.C, seed.scale.WxC)
str(seed.scale)



# Q-value - correction for multiple testing ==============================================================================================================================================================================================================================================

# check p-value distribution

par(mfrow = c(1, 3))
hist(seed.permanova$permanova.W.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Water regime", xlim = c(0, 1))
hist(seed.permanova$permanova.C.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Cropping system", xlim = c(0, 1))
hist(seed.permanova$permanova.WxC.P, breaks = seq(0, 1, 0.05), xlab = "p-value", main = "Treatment x System", xlim = c(0, 1))

# apply FDR

FDR <- qvalue(seed.permanova$permanova.W.P)  # check fit with 'hist(FDR)'
seed.permanova$permanova.W.q <- FDR$qvalues
FDR <- qvalue(seed.permanova$permanova.C.P)  # check fit with 'hist(FDR)'
seed.permanova$permanova.C.q <- FDR$qvalues
FDR <- qvalue(seed.permanova$permanova.WxC.P, pi0 = 0.65)  # check fit with 'hist(FDR)'
seed.permanova$permanova.WxC.q <- FDR$qvalues



# Combine files ==============================================================================================================================================================================================================================================

seed.gene.level.analysis <- cbind(seed.scale, seed.permanova, seed.permdisp)
seed.gene.level.analysis$count <- colSums(seed.gene.summary)
seed.gene.level.analysis$relabund <- 100/sum(colSums(seed.gene.summary)) * colSums(seed.gene.summary)
str(seed.gene.level.analysis)



# Check results ==============================================================================================================================================================================================================================================

# Check results - water regime ==============

seed.gene.level.analysis %>% filter(grepl("l1:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) # first level

seed.gene.level.analysis %>% filter(grepl("_l2:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) # second level

seed.gene.level.analysis %>% filter(grepl("__l3:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) # third level

seed.gene.level.analysis %>% filter(grepl("___l4:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05) # fourth level



# Check results - water regime x cropping system ==============

seed.gene.level.analysis %>% filter(grepl("l1:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.WxC.q < 0.05) # first level

seed.gene.level.analysis %>% filter(grepl("_l2:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.WxC.q < 0.05) # second level

seed.gene.level.analysis %>% filter(grepl("__l3:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.WxC.q < 0.05) # third level

seed.gene.level.analysis %>% filter(grepl("___l4:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.WxC.q < 0.05) # fourth level



# Visualization ==============================================================================================================================================================================================================================================

# prepare files with levels ==============

seed.genes.filtered <- seed.genes

seed.genes.filtered$level_2 <- gsub(",", "", seed.genes.filtered$level_2) # remove commas
seed.genes.filtered$level_3 <- gsub(",", "", seed.genes.filtered$level_3) # remove commas
seed.genes.filtered$level_4 <- gsub(",", "", seed.genes.filtered$level_4) # remove commas

seed.genes.filtered$level_2 <- gsub(" ", "_", seed.genes.filtered$level_2) # remove line
seed.genes.filtered$level_3 <- gsub(" ", "_", seed.genes.filtered$level_3) # remove line
seed.genes.filtered$level_4 <- gsub(" ", "_", seed.genes.filtered$level_4) # remove line

seed.genes.filtered$level_2 <- gsub("\\(", "", seed.genes.filtered$level_2) # remove brackets
seed.genes.filtered$level_3 <- gsub("\\(", "", seed.genes.filtered$level_3) # remove brackets
seed.genes.filtered$level_4 <- gsub("\\(", "", seed.genes.filtered$level_4) # remove brackets

seed.genes.filtered$level_2 <- gsub(")", "", seed.genes.filtered$level_2) # remove brackets
seed.genes.filtered$level_3 <- gsub(")", "", seed.genes.filtered$level_3) # remove brackets
seed.genes.filtered$level_4 <- gsub(")", "", seed.genes.filtered$level_4) # remove brackets

seed.genes.filtered$level_2 <- gsub(":", "", seed.genes.filtered$level_2) # remove colon
seed.genes.filtered$level_3 <- gsub(":", "", seed.genes.filtered$level_3) # remove colon
seed.genes.filtered$level_4 <- gsub(":", "", seed.genes.filtered$level_4) # remove colon

seed.genes.filtered$level_2 <- gsub("\\[", "", seed.genes.filtered$level_2) # remove brackets
seed.genes.filtered$level_3 <- gsub("\\[", "", seed.genes.filtered$level_3) # remove brackets
seed.genes.filtered$level_4 <- gsub("\\[", "", seed.genes.filtered$level_4) # remove brackets

seed.genes.filtered$level_2 <- gsub("]", "", seed.genes.filtered$level_2) # remove brackets
seed.genes.filtered$level_3 <- gsub("]", "", seed.genes.filtered$level_3) # remove brackets
seed.genes.filtered$level_4 <- gsub("]", "", seed.genes.filtered$level_4) # remove brackets

# filter for specific branches

seed.df1 <- seed.genes.filtered %>%
  filter(level_3 == "Capsule_and_Slime_layer" |
           level_3 == "Stress_Response_Osmotic_stress" |
           level_2 == "Betaine_crotonobetaine_L-carnitine_trimethyl_ammonium_compounds" |
           level_2 == "DNA_repair" |
           level_2 == "Uni-_Sym-_and_Antiporters" |
           level_2 == "ABC_transporters" |
           level_2 == "Protein_secretion_system" |
           level_3 == "Sporulation" |
           level_3 == "Programmed_Cell_Death_and_Toxin-antitoxin_Systems" |
           level_3 == "Stationary_phase_Dormancy_Persistence" |
           level_3 == "Cell_division_initiation_related_cluster" |
           level_3 == "Control_of_cell_elongation_-_division_cycle_in_Bacilli"  |
           level_2 == "Secondary_Metabolism" |
           level_2 == "Phosphate_Metabolism" |
           level_2 == "Iron_acquisition_and_metabolism" |
           level_3 == "Phospholipids" |
           level_3 == "Isoprenoids" |
           level_3 == "Unsaturated_fatty_acids" |
           level_3 == "Thiamin" |
           level_3 == "Tetrapyrroles" |
           level_3 == "Riboflavin_FMN_FAD" |
           level_3 == "NAD_and_NADP" |
           level_3 == "Mycofactocin" |
           level_3 == "Fe-S_clusters" |
           level_3 == "Coenzyme_F420" |
           level_3 == "Proline_and_4-hydroxyproline" )

seed.df2 <- seed.df1 %>%
  mutate(name = if_else(level_4 == "unclassified", level_3, level_4)) %>%
  mutate(name = if_else(name == "unclassified", level_2, name)) %>%
  mutate(level_3_new = if_else(level_3 == "unclassified", level_2, level_3)) %>%
  select(-level_3) %>%
  rename("level_3" = "level_3_new")

seed.df3 <- seed.df2 %>% relocate(name) %>% column_to_rownames("name") %>%
  select(-level_4) 

# prepare tree

seed.df4 <- taxa2dist(seed.df3, varstep = TRUE)
hc <- hclust(seed.df4)
hcp <- as.phylo(hc)
write.tree(hcp, file = "iTol.SEED.figure2.newick", append = FALSE, digits = 10, tree.names = FALSE)



# prepare files with values ==============

# filter for specific gene functions for 2nd level

seed.l2.df1 <- seed.gene.level.analysis %>% filter(grepl("_l2:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05)

rownames(seed.l2.df1) <- sub("_l2:*", "", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub(",", "", rownames(seed.l2.df1)) 
rownames(seed.l2.df1) <- gsub(" ", "_", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub("\\(", "", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub(")", "", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub(":", "", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub("\\[", "", rownames(seed.l2.df1))
rownames(seed.l2.df1) <- gsub("]", "", rownames(seed.l2.df1))

seed.l2.df2 <- seed.l2.df1 %>% mutate(rainout.mean = rainout.W.mean - control.W.mean, control.mean = control.W.mean - control.W.mean)
seed.l2.df2 <- tibble::rownames_to_column(seed.l2.df2, "name")

seed.l2.df3 <- seed.l2.df2 %>%
  dplyr::select(name, control.W.mean, rainout.W.mean, permanova.W.q) %>%
  column_to_rownames("name") %>%
  mutate(value_1 = 0.5) %>%
  mutate(value_2 = if_else(control.W.mean > 0 & permanova.W.q < 0.05,  control.W.mean, if_else(rainout.W.mean > 0 & permanova.W.q < 0.05, rainout.W.mean, 0))) %>%
  mutate(value_3 = if_else(control.W.mean > 0 & permanova.W.q < 0.05,  control.W.mean, 0)) %>%
  mutate(value_4 = if_else(rainout.W.mean > 0 & permanova.W.q < 0.05,  rainout.W.mean, 0)) %>%
  select(value_1, value_2, value_3, value_4)

# filter for specific gene functions for 3rd level

seed.l3.df1 <- seed.gene.level.analysis %>% filter(grepl("__l3:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05)

rownames(seed.l3.df1) <- sub("__l3:*", "", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub(",", "", rownames(seed.l3.df1)) 
rownames(seed.l3.df1) <- gsub(" ", "_", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub("\\(", "", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub(")", "", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub(":", "", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub("\\[", "", rownames(seed.l3.df1))
rownames(seed.l3.df1) <- gsub("]", "", rownames(seed.l3.df1))

seed.l3.df2 <- seed.l3.df1 %>% mutate(rainout.mean = rainout.W.mean - control.W.mean, control.mean = control.W.mean - control.W.mean)
seed.l3.df2 <- tibble::rownames_to_column(seed.l3.df2, "name")

seed.l3.df3 <- seed.l3.df2 %>%
  dplyr::select(name, control.W.mean, rainout.W.mean, permanova.W.q) %>%
  column_to_rownames("name") %>%
  mutate(value_1 = 0.5) %>%
  mutate(value_2 = if_else(control.W.mean > 0 & permanova.W.q < 0.05, control.W.mean, if_else(rainout.W.mean > 0 & permanova.W.q < 0.05, rainout.W.mean, 0))) %>%
  mutate(value_3 = if_else(control.W.mean > 0 & permanova.W.q < 0.05, control.W.mean, 0)) %>%
  mutate(value_4 = if_else(rainout.W.mean > 0 & permanova.W.q < 0.05, rainout.W.mean, 0)) %>%
  select(value_1, value_2, value_3, value_4)


# filter for specific gene functions for 4th level

seed.l4.df1 <- seed.gene.level.analysis %>% filter(grepl("___l4:", rownames(seed.gene.level.analysis))) %>%
  filter(permanova.W.q < 0.05)

rownames(seed.l4.df1) <- sub("___l4:*", "", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub(",", "", rownames(seed.l4.df1)) 
rownames(seed.l4.df1) <- gsub(" ", "_", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub("\\(", "", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub(")", "", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub(":", "", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub("\\[", "", rownames(seed.l4.df1))
rownames(seed.l4.df1) <- gsub("]", "", rownames(seed.l4.df1))

seed.l4.df2 <- seed.l4.df1 %>% mutate(rainout.mean = rainout.W.mean - control.W.mean, control.mean = control.W.mean - control.W.mean)
seed.l4.df2 <- tibble::rownames_to_column(seed.l4.df2, "name")

seed.l4.df3 <- seed.l4.df2 %>%
  dplyr::select(name, control.W.mean, rainout.W.mean, permanova.W.q) %>%
  column_to_rownames("name") %>%
  mutate(value_1 = 0.5) %>%
  mutate(value_2 = if_else(control.W.mean > 0 & permanova.W.q < 0.05, control.W.mean, if_else(rainout.W.mean > 0 & permanova.W.q < 0.05, rainout.W.mean, 0))) %>%
  mutate(value_3 = if_else(control.W.mean > 0 & permanova.W.q < 0.05, control.W.mean, 0)) %>%
  mutate(value_4 = if_else(rainout.W.mean > 0 & permanova.W.q < 0.05, rainout.W.mean, 0)) %>%
  select(value_1, value_2, value_3, value_4)


seed.relab <- rbind(seed.l2.df3, seed.l3.df3, seed.l4.df3)

write.table(seed.relab, file = "relab.SEED.txt", sep = ",") # write txt file



