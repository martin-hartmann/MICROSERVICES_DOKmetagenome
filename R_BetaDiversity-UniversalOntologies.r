
# R Script for analysis of the metagenome samples of the DOK rainout-trial 2022
# data normalized by open reading frame length
# Elena Kost
# 19.04.24

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, plotrix, lmerTest, geomtextpath, 
               lsmeans, multcomp)
library(pacman)



# Import data ==============================================================================================================================================================================================================================================

seed <- read.delim("SEED_RPKmaxORF.tsv",sep="\t")
colnames(seed) <- sub("S\\.", "S-", colnames(seed))
str(seed)
seed[,-c(1:5)][is.na(seed[,-c(1:5)])] <- 0

eggnog <- read.delim("EGGNOG_RPKmaxORF.tsv",sep="\t")
colnames(eggnog) <- sub("S\\.", "S-", colnames(eggnog))
str(eggnog)
eggnog[,-c(1:4)][is.na(eggnog[,-c(1:4)])] <- 0

inter <- read.delim("INTERPRO2GO_RPKmaxORF.tsv",sep="\t")
colnames(inter) <- sub("S\\.", "S-", colnames(inter))
str(inter)
inter[,-c(1:9)][is.na(inter[,-c(1:9)])] <- 0

ec <- read.delim("EC_RPKmaxORF.tsv",sep="\t")
colnames(ec) <- sub("S\\.", "S-", colnames(ec))
str(ec)
ec[,-c(1:5)][is.na(ec[,-c(1:5)])] <- 0



# Import Design ==============================================================================================================================================================================================================================================

design <- read.table("design_metagenomic.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE,
                     colClasses = "factor")
str(design)

identical(rownames(design), colnames(seed[,-c(1:5)]))
identical(rownames(design), colnames(eggnog[,-c(1:4)]))
identical(rownames(design), colnames(inter[,-c(1:9)]))
identical(rownames(design), colnames(ec[,-c(1:5)]))

design$Pairs <- c("P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", "P7", "P8", "P8", "P9", "P9", "P10", "P10", "P11", "P11", "P12", "P12",
                  "P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", "P7", "P8", "P8", "P9", "P9", "P10", "P10", "P11", "P11", "P12", "P12",
                  "P1", "P1", "P2", "P2", "P3", "P3", "P4", "P4", "P5", "P5", "P6", "P6", "P7", "P7", "P8", "P8", "P9", "P9", "P10", "P10", "P11", "P11", "P12", "P12")

design$Pairs <- as.factor(design$Pairs)

design$IxT <- with(design, interaction(Irrigation, Treatment, lex.order = TRUE, drop = TRUE))

design$Date <- factor(design$Date, levels = c("28.04.22", "01.06.22", "05.07.22"))



# Check rowsums of data files ==============================================================================================================================================================================================================================================

scientific_10 <- function(x) {parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))}

rowsum.seed <- ggplot(data = as.data.frame(t(seed[,-c(1:5)])), aes(x= rownames(as.data.frame(t(seed[,-c(1:5)]))), y = rowSums(as.data.frame(t(seed[,-c(1:5)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

rowsum.eggnog <- ggplot(data = as.data.frame(t(eggnog[,-c(1:4)])), aes(x= rownames(as.data.frame(t(eggnog[,-c(1:4)]))), y = rowSums(as.data.frame(t(eggnog[,-c(1:4)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

rowsum.inter <- ggplot(data = as.data.frame(t(inter[,-c(1:9)])), aes(x= rownames(as.data.frame(t(inter[,-c(1:9)]))), y = rowSums(as.data.frame(t(inter[,-c(1:9)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

rowsum.ec <- ggplot(data = as.data.frame(t(ec[,-c(1:5)])), aes(x= rownames(as.data.frame(t(ec[,-c(1:5)]))), y = rowSums(as.data.frame(t(ec[,-c(1:5)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

ggarrange(rowsum.seed, rowsum.eggnog, rowsum.inter, rowsum.ec)



# Anova ================

model.seed <- aov(log(rowSums(as.data.frame(t(seed[,-c(1:5)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.seed)
shapiro.test(residuals(model.seed))
leveneTest(model.seed)

model.eggnog <- aov(log(rowSums(as.data.frame(t(eggnog[,-c(1:4)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.eggnog)
shapiro.test(residuals(model.eggnog))
leveneTest(model.eggnog)

model.inter <- aov(log(rowSums(as.data.frame(t(inter[,-c(1:9)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.inter)
shapiro.test(residuals(model.inter))
leveneTest(model.inter)

model.ec <- aov(log(rowSums(as.data.frame(t(ec[,-c(1:5)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.ec)
shapiro.test(residuals(model.ec))
leveneTest(model.ec)



# Normalization of data ==============================================================================================================================================================================================================================================

# rarefying with 100 

iters <- 100
seed.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(seed[,-c(1:5)])),  min(rowSums(as.data.frame(t(seed[,-c(1:5)]))))), mc.cores = detectCores())

eggnog.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(eggnog[,-c(1:4)])),  min(rowSums(as.data.frame(t(eggnog[,-c(1:4)]))))),  mc.cores = detectCores())

inter.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(inter[,-c(1:9)])),  min(rowSums(as.data.frame(t(inter[,-c(1:9)]))))),  mc.cores = detectCores())

ec.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(ec[,-c(1:5)])),  min(rowSums(as.data.frame(t(ec[,-c(1:5)]))))),  mc.cores = detectCores())



# Beta diversity ==============================================================================================================================================================================================================================================

# Bray curtis dissimilarity ============

# SEED

ISS.iters.bray.seed <- mclapply(seed.iters, function(x) vegdist(x, method = "bray"))
ISS.array.bray.seed <- laply(ISS.iters.bray.seed, as.matrix)
seed.bray <- as.dist(apply(seed.array.bray, 2:3, median))

# EGGnog

ISS.iters.bray.eggnog <- mclapply(eggnog.iters, function(x) vegdist(x, method = "bray"))
ISS.array.bray.eggnog <- laply(ISS.iters.bray.eggnog, as.matrix)
eggnog.bray <- as.dist(apply(eggnog.array.bray, 2:3, median))


# InterPro2GO

ISS.iters.bray.inter <- mclapply(inter.iters, function(x) vegdist(x, method = "bray"))
ISS.array.bray.inter <- laply(ISS.iters.bray.inter, as.matrix)
inter.bray <- as.dist(apply(inter.array.bray, 2:3, median))

# EC

ISS.iters.bray.ec <- mclapply(ec.iters, function(x) vegdist(x, method = "bray"))
ISS.array.bray.ec <- laply(ISS.iters.bray.ec, as.matrix)
ec.bray <- as.dist(apply(ec.array.bray, 2:3, median))



# Create principal coordinates ============

# SEED

ISS.pco.seed <- cmdscale(seed.bray, eig = TRUE, add = TRUE)
ISS.pco1.seed <- (ISS.pco.seed$eig/sum(ISS.pco.seed$eig))[1]
ISS.pco2.seed <- (ISS.pco.seed$eig/sum(ISS.pco.seed$eig))[2]
ISS.pco3.seed <- (ISS.pco.seed$eig/sum(ISS.pco.seed$eig))[3]
ISS.pco1.seed <- paste("PCO1 (", format(round(ISS.pco1.seed * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco2.seed <- paste("PCO2 (", format(round(ISS.pco2.seed * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco3.seed <- paste("PCO3 (", format(round(ISS.pco3.seed * 100, 1), nsmall = 1), "%)", sep = "")

# EGGnog

ISS.pco.eggnog <- cmdscale(eggnog.bray, eig = TRUE, add = TRUE)
ISS.pco1.eggnog <- (ISS.pco.eggnog$eig/sum(ISS.pco.eggnog$eig))[1]
ISS.pco2.eggnog <- (ISS.pco.eggnog$eig/sum(ISS.pco.eggnog$eig))[2]
ISS.pco3.eggnog <- (ISS.pco.eggnog$eig/sum(ISS.pco.eggnog$eig))[3]
ISS.pco1.eggnog <- paste("PCO1 (", format(round(ISS.pco1.eggnog * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco2.eggnog <- paste("PCO2 (", format(round(ISS.pco2.eggnog * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco3.eggnog <- paste("PCO3 (", format(round(ISS.pco3.eggnog * 100, 1), nsmall = 1), "%)", sep = "")

# Interpro2GO

ISS.pco.inter <- cmdscale(inter.bray, eig = TRUE, add = TRUE)
ISS.pco1.inter <- (ISS.pco.inter$eig/sum(ISS.pco.inter$eig))[1]
ISS.pco2.inter <- (ISS.pco.inter$eig/sum(ISS.pco.inter$eig))[2]
ISS.pco3.inter <- (ISS.pco.inter$eig/sum(ISS.pco.inter$eig))[3]
ISS.pco1.inter <- paste("PCO1 (", format(round(ISS.pco1.inter * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco2.inter <- paste("PCO2 (", format(round(ISS.pco2.inter * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco3.inter <- paste("PCO3 (", format(round(ISS.pco3.inter * 100, 1), nsmall = 1), "%)", sep = "")

# EC

ISS.pco.ec <- cmdscale(ec.bray, eig = TRUE, add = TRUE)
ISS.pco1.ec <- (ISS.pco.ec$eig/sum(ISS.pco.ec$eig))[1]
ISS.pco2.ec <- (ISS.pco.ec$eig/sum(ISS.pco.ec$eig))[2]
ISS.pco3.ec <- (ISS.pco.ec$eig/sum(ISS.pco.ec$eig))[3]
ISS.pco1.ec <- paste("PCO1 (", format(round(ISS.pco1.ec * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco2.ec <- paste("PCO2 (", format(round(ISS.pco2.ec * 100, 1), nsmall = 1), "%)", sep = "")
ISS.pco3.ec <- paste("PCO3 (", format(round(ISS.pco3.ec * 100, 1), nsmall = 1), "%)", sep = "")



# Visualization and analysis of principal coordinates ============

# SEED

ggplot(as.data.frame(ISS.pco.seed$points), aes(x = ISS.pco.seed$points[,1], y = ISS.pco.seed$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.pco1.seed) + ylab(ISS.pco2.seed) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(seed.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999) # Permanova
permutest(betadisper(seed.bray, design$Irrigation), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(seed.bray, design$Treatment), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(seed.bray, design$Date), permutations = how(nperm = 9999)) # Permdisp


# EGGnog

ggplot(as.data.frame(ISS.pco.eggnog$points), aes(x = ISS.pco.eggnog$points[,1], y = ISS.pco.eggnog$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.pco1.eggnog) + ylab(ISS.pco2.eggnog) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(eggnog.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999) # Permanova
permutest(betadisper(eggnog.bray, design$Irrigation), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(eggnog.bray, design$Treatment), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(eggnog.bray, design$Date), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(eggnog.bray, design$IxT), permutations = how(nperm = 9999)) # Permdisp


# Interpro2GO

ggplot(as.data.frame(ISS.pco.inter$points), aes(x = ISS.pco.inter$points[,1], y = ISS.pco.inter$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.pco1.inter) + ylab(ISS.pco2.inter) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(inter.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999) # Permanova
permutest(betadisper(inter.bray, design$Irrigation), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(inter.bray, design$Treatment), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(inter.bray, design$Date), permutations = how(nperm = 9999)) # Permdisp


# EC

ggplot(as.data.frame(ISS.pco.ec$points), aes(x = ISS.pco.ec$points[,1], y = ISS.pco.ec$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.pco1.ec) + ylab(ISS.pco2.ec) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(ec.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999) # Permanova
permutest(betadisper(ec.bray, design$Irrigation), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(ec.bray, design$Treatment), permutations = how(nperm = 9999)) # Permdisp
permutest(betadisper(ec.bray, design$Date), permutations = how(nperm = 9999)) # Permdisp



# Canonical analysis of principal coordinates (CAP) ============

# SEED - irrigation and cropping system

nc <- nrow(as.matrix(seed.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(seed.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

ISS.cap.seed <- CAPdiscrim(seed.bray ~ IxT, data = design, m = 12, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(ISS.cap.seed$group), data.frame(ISS.cap.seed$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.seed$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.seed <- paste("CAP1 (", round((100/sum(ISS.cap.seed$lda.other$svd^2) * ISS.cap.seed$lda.other$svd^2)[1],
                                       digits = 1), "%)", sep = "")
ISS.cap2.seed <- paste("CAP2 (", round((100/sum(ISS.cap.seed$lda.other$svd^2) * ISS.cap.seed$lda.other$svd^2)[2],
                                       digits = 1), "%)", sep = "")


# EGGnog - irrigation and cropping system

nc <- nrow(as.matrix(eggnog.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(eggnog.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

ISS.cap.eggnog <- CAPdiscrim(eggnog.bray ~ IxT, data = design, m = 17, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(ISS.cap.eggnog$group), data.frame(ISS.cap.eggnog$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.eggnog$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.eggnog <- paste("CAP1 (", round((100/sum(ISS.cap.eggnog$lda.other$svd^2) * ISS.cap.eggnog$lda.other$svd^2)[1],
                                         digits = 1), "%)", sep = "")
ISS.cap2.eggnog <- paste("CAP2 (", round((100/sum(ISS.cap.eggnog$lda.other$svd^2) * ISS.cap.eggnog$lda.other$svd^2)[2],
                                         digits = 1), "%)", sep = "")


# Interpro2GO - irrigation and cropping system

nc <- nrow(as.matrix(inter.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(inter.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

ISS.cap.inter <- CAPdiscrim(inter.bray ~ IxT, data = design, m = 34, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(ISS.cap.inter$group), data.frame(ISS.cap.inter$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.inter$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.inter <- paste("CAP1 (", round((100/sum(ISS.cap.inter$lda.other$svd^2) * ISS.cap.inter$lda.other$svd^2)[1],
                                        digits = 1), "%)", sep = "")
ISS.cap2.inter <- paste("CAP2 (", round((100/sum(ISS.cap.inter$lda.other$svd^2) * ISS.cap.inter$lda.other$svd^2)[2],
                                        digits = 1), "%)", sep = "")


# EC - irrigation and cropping system

nc <- nrow(as.matrix(ec.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(ec.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

ISS.cap.ec <- CAPdiscrim(ec.bray ~ IxT, data = design, m = 19, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(ISS.cap.ec$group), data.frame(ISS.cap.ec$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(ISS.cap.ec$PCoA)
success <- success[order(success$source), ]
success

ISS.cap1.ec <- paste("CAP1 (", round((100/sum(ISS.cap.ec$lda.other$svd^2) * ISS.cap.ec$lda.other$svd^2)[1],
                                     digits = 1), "%)", sep = "")
ISS.cap2.ec <- paste("CAP2 (", round((100/sum(ISS.cap.ec$lda.other$svd^2) * ISS.cap.ec$lda.other$svd^2)[2],
                                     digits = 1), "%)", sep = "")



# Visualization of CAP ============

# SEED

plot.cap.seed <- ggplot(as.data.frame(ISS.cap.seed$x), aes(x = ISS.cap.seed$x[,1], y = ISS.cap.seed$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.cap1.seed) + ylab(ISS.cap2.seed) +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") +
  geom_mark_ellipse(aes(fill = design$IxT), 
                    expand = 0, linewidth = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  theme_classic() + theme(legend.position="bottom", legend.box = "horizontal") +
  
  annotate("text", x= -7, y= -7,label= "Overall reclassification rate: 92%", hjust = 0, size = 4) +
  annotate("text", x= -7, y= -8,label= "Pillais test=3.0 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -6, y=  3,label= "100%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -4, y= -2,label= "92%", hjust = 0, size = 4, color = "#009E73") +
  
  annotate("text", x=  4, y= 5,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x= -1, y= 0.5,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= 5, y= -3.5,label= "83%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= 0, y= -4.5,label= "92%", hjust = 0, size = 4, color = "#E69F00") 


# EGGnog

plot.cap.eggnog <- ggplot(as.data.frame(ISS.cap.eggnog$x), aes(x = ISS.cap.eggnog$x[,1], y = ISS.cap.eggnog$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.cap1.eggnog) + ylab(ISS.cap2.eggnog) +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") +
  geom_mark_ellipse(aes(fill = design$IxT), 
                    expand = 0, linewidth = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  theme_classic() + theme(legend.position="bottom", legend.box = "horizontal") +
  
  annotate("text", x= -7, y= -7.5,label= "Overall reclassification rate: 90%", hjust = 0, size = 4) +
  annotate("text", x= -7, y= -8.5,label= "Pillais test=3.5 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -3, y=  0,label= "100%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -4, y= 4,label= "92%", hjust = 0, size = 4, color = "#009E73")  +
  
  annotate("text", x=  2, y= -2,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x= -4, y= -5,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= 3, y= 4,label= "83%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= 3, y= 0,label= "83%", hjust = 0, size = 4, color = "#E69F00") 


# Interpro2GO

plot.cap.inter <- ggplot(as.data.frame(ISS.cap.inter$x), aes(x = ISS.cap.inter$x[,1], y = ISS.cap.inter$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.cap1.inter) + ylab(ISS.cap2.inter) +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") +
  geom_mark_ellipse(aes(fill = design$IxT), 
                    expand = 0, linewidth = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  theme_classic() + theme(legend.position="bottom", legend.box = "horizontal") +
  
  annotate("text", x= -13, y= -11,label= "Overall reclassification rate: 88%", hjust = 0, size = 4) +
  annotate("text", x= -13, y= -12,label= "Pillais test=4.1 ***", hjust = 0, size = 4) +
  
  annotate("text", x= - 5, y= 1,label= "67%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -12, y= 4.5,label= "92%", hjust = 0, size = 4, color = "#009E73")  +
  
  annotate("text", x= -2.5, y= -5,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x=  5.5, y= -6.5,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= 3.5, y= 8,label= "92%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= 3, y= 4.5,label= "92%", hjust = 0, size = 4, color = "#E69F00") 


# EC

plot.cap.ec <- ggplot(as.data.frame(ISS.cap.ec$x), aes(x = ISS.cap.ec$x[,1], y = ISS.cap.ec$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(ISS.cap1.ec) + ylab(ISS.cap2.ec) +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") +
  geom_mark_ellipse(aes(fill = design$IxT), 
                    expand = 0, linewidth = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("#009E73", "#FF618C", "#E69F00", "#009E73", "#FF618C", "#E69F00")) +
  theme_classic() + theme(legend.position="bottom", legend.box = "horizontal") +
  
  annotate("text", x= -10, y= -5,label= "Overall reclassification rate: 97%", hjust = 0, size = 4) +
  annotate("text", x= -10, y= -6,label= "Pillais test=3.5 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -7, y=  1,label= "92%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -10, y= -3,label= "92%", hjust = 0, size = 4, color = "#009E73")  +
  
  annotate("text", x= -2, y= 5,label= "100%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x=  5.5, y= 3,label= "100%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= 3, y= -1,label= "100%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= 0.5, y= -5,label= "100%", hjust = 0, size = 4, color = "#E69F00") 


# combine

ggarrange(plot.cap.seed, plot.cap.ec, plot.cap.inter, plot.cap.eggnog, ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom", labels = c("A", "B", "C", "D"))



