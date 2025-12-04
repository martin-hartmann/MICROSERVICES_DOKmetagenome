
# R Script for Analysis of the metagenome samples of the DOK rainout-trial 2022 
# composition analysis for nutrient cycling (cazy, ncyc, pcyc)
# Elena Kost
# 19.09.24

# Load packages ==============================================================================================================================================================================================================================================

pacman::p_load(BiocManager, bestNormalize, BiodiversityR, car, compositions, data.table, 
               dbstats, devtools, DESeq2, Hmisc, indicspecies, matrixStats, parallel, plyr, 
               qvalue, reshape2, RVAideMemoire, tidyverse, ggpubr, agricolae, interp, lubridate,
               ggforce, multcompView, vegan, RVAideMemoire, scales, plyr, Nonpareil)
library(pacman)



# Import raw files  ==============================================================================================================================================================================================================================================

cazy <- read.delim("CAZy_RPKmaxORF.tsv",sep="\t")
colnames(cazy) <- sub("S\\.", "S-", colnames(cazy))
str(cazy)
cazy[,-c(1:3)][is.na(cazy[,-c(1:3)])] <- 0

ncyc <- read.delim("NCyc_RPKmaxORF.tsv",sep="\t")
colnames(ncyc) <- sub("S\\.", "S-", colnames(ncyc))
str(ncyc)
ncyc[,-c(1:3)][is.na(ncyc[,-c(1:3)])] <- 0

pcyc <- read.delim("PCyc_RPKmaxORF.tsv",sep="\t")
colnames(pcyc) <- sub("S\\.", "S-", colnames(pcyc))
str(pcyc)
pcyc[,-c(1:3)][is.na(pcyc[,-c(1:3)])] <- 0



# Import Design ==============================================================================================================================================================================================================================================

design <- read.table("design_metagenomic.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = TRUE,
                     colClasses = "factor")
str(design)

design$IxT <- with(design, interaction(Irrigation, Treatment, lex.order = TRUE, drop = TRUE))

identical(rownames(design), colnames(cazy[,-c(1:3)]))
identical(rownames(design), colnames(ncyc[,-c(1:3)]))
identical(rownames(design), colnames(pcyc[,-c(1:3)]))



# Check rowsums of data files ==============================================================================================================================================================================================================================================

scientific_10 <- function(x) {parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))}

rowsum.cazy <- ggplot(data = as.data.frame(t(cazy[,-c(1:3)])), aes(x= rownames(as.data.frame(t(cazy[,-c(1:3)]))), y = rowSums(as.data.frame(t(cazy[,-c(1:3)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

rowsum.ncyc <- ggplot(data = as.data.frame(t(ncyc[,-c(1:3)])), aes(x= rownames(as.data.frame(t(ncyc[,-c(1:3)]))), y = rowSums(as.data.frame(t(ncyc[,-c(1:3)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

rowsum.pcyc <- ggplot(data = as.data.frame(t(pcyc[,-c(1:3)])), aes(x= rownames(as.data.frame(t(pcyc[,-c(1:3)]))), y = rowSums(as.data.frame(t(pcyc[,-c(1:3)]))))) +
  geom_bar(stat="identity") + theme_classic() + scale_y_continuous(label=scientific_10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")

ggarrange(rowsum.cazy, rowsum.ncyc, rowsum.pcyc, ncol = 3)



# Anova ================

model.cazy <- aov(log(rowSums(as.data.frame(t(cazy[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.cazy)
shapiro.test(residuals(model.cazy))
leveneTest(model.cazy)

model.ncyc <- aov(log(rowSums(as.data.frame(t(ncyc[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.ncyc)
shapiro.test(residuals(model.ncyc))
leveneTest(model.ncyc)

model.pcyc <- aov(log(rowSums(as.data.frame(t(pcyc[,-c(1:3)])))) ~ Irrigation * Treatment * Date, data = design)
summary(model.pcyc)
shapiro.test(residuals(model.pcyc))
leveneTest(model.pcyc)



# Normalization of data ==============================================================================================================================================================================================================================================

# rarefying with 100 
iters <- 100

cazy.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(cazy[,-c(1:4)])),  min(rowSums(as.data.frame(t(cazy[,-c(1:5)]))))), mc.cores = detectCores())

ncyc.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(ncyc[,-c(1:3)])),  min(rowSums(as.data.frame(t(ncyc[,-c(1:4)]))))),  mc.cores = detectCores())

pcyc.iters <- mclapply(as.list(1:iters), function(x) rrarefy(as.data.frame(t(pcyc[,-c(1:3)])),  min(rowSums(as.data.frame(t(pcyc[,-c(1:9)]))))),  mc.cores = detectCores())



# Beta diversity ==============================================================================================================================================================================================================================================

# Bray curtis dissimilarity ============

# CAZyme

iters.bray.cazy <- mclapply(cazy.iters, function(x) vegdist(x, method = "bray"))
array.bray.cazy <- laply(iters.bray.cazy, as.matrix)
cazy.bray <- as.dist(apply(array.bray.cazy, 2:3, median))

# NCycDB

iters.bray.ncyc <- mclapply(ncyc.iters, function(x) vegdist(x, method = "bray"))
array.bray.ncyc <- laply(iters.bray.ncyc, as.matrix)
ncyc.bray <- as.dist(apply(array.bray.ncyc, 2:3, median))

# PCycDB

iters.bray.pcyc <- mclapply(pcyc.iters, function(x) vegdist(x, method = "bray"))
array.bray.pcyc <- laply(iters.bray.pcyc, as.matrix)
pcyc.bray <- as.dist(apply(array.bray.pcyc, 2:3, median))


# Create principal coordinates ============

# CAZyme

pco.cazy <- cmdscale(cazy.bray, eig = TRUE, add = TRUE)
pco1.cazy <- (pco.cazy$eig/sum(pco.cazy$eig))[1]
pco2.cazy <- (pco.cazy$eig/sum(pco.cazy$eig))[2]
pco1.cazy <- paste("PCO1 (", format(round(pco1.cazy * 100, 1), nsmall = 1), "%)", sep = "")
pco2.cazy <- paste("PCO2 (", format(round(pco2.cazy * 100, 1), nsmall = 1), "%)", sep = "")

# NCycDB

pco.ncyc <- cmdscale(ncyc.bray, eig = TRUE, add = TRUE)
pco1.ncyc <- (pco.ncyc$eig/sum(pco.ncyc$eig))[1]
pco2.ncyc <- (pco.ncyc$eig/sum(pco.ncyc$eig))[2]
pco1.ncyc <- paste("PCO1 (", format(round(pco1.ncyc * 100, 1), nsmall = 1), "%)", sep = "")
pco2.ncyc <- paste("PCO2 (", format(round(pco2.ncyc * 100, 1), nsmall = 1), "%)", sep = "")

# PCycDB

pco.pcyc <- cmdscale(pcyc.bray, eig = TRUE, add = TRUE)
pco1.pcyc <- (pco.pcyc$eig/sum(pco.pcyc$eig))[1]
pco2.pcyc <- (pco.pcyc$eig/sum(pco.pcyc$eig))[2]
pco1.pcyc <- paste("PCO1 (", format(round(pco1.pcyc * 100, 1), nsmall = 1), "%)", sep = "")
pco2.pcyc <- paste("PCO2 (", format(round(pco2.pcyc * 100, 1), nsmall = 1), "%)", sep = "")


# Visualization and analysis of principal coordinates ============

colors <- c("#009E73", "#FF618C", "#E69F00") 

# CAZyme

ggplot(as.data.frame(pco.cazy$points), aes(x = pco.cazy$points[,1], y = pco.cazy$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(pco1.cazy) + ylab(pco2.cazy) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(cazy.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999)
permutest(betadisper(cazy.bray, design$Irrigation), permutations = how(nperm = 9999))
permutest(betadisper(cazy.bray, design$Treatment), permutations = how(nperm = 9999))
permutest(betadisper(cazy.bray, design$Date), permutations = how(nperm = 9999))


# NCycDB

ggplot(as.data.frame(pco.ncyc$points), aes(x = pco.ncyc$points[,1], y = pco.ncyc$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(pco1.ncyc) + ylab(pco2.ncyc) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(ncyc.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999)
permutest(betadisper(ncyc.bray, design$Irrigation), permutations = how(nperm = 9999))
permutest(betadisper(ncyc.bray, design$Treatment), permutations = how(nperm = 9999))


# PCycDB

ggplot(as.data.frame(pco.pcyc$points), aes(x = pco.pcyc$points[,1], y = pco.pcyc$points[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(pco1.pcyc) + ylab(pco2.pcyc) +
  theme_classic() +
  scale_color_manual(values =  colors, 
                     name = "Cropping system",
                     labels = c("BIODYN", "CONFYM", "CONMIN")) +
  scale_shape_manual(values = c(8, 17),
                     name = "Water regime",
                     labels = c("rainfed control", "drought-induced")) +
  theme(legend.position = "bottom") 

adonis2(pcyc.bray ~ Irrigation * Treatment * Date, data = design, permutations = 9999)
permutest(betadisper(pcyc.bray, design$Irrigation), permutations = how(nperm = 9999))
permutest(betadisper(pcyc.bray, design$Treatment), permutations = how(nperm = 9999))
permutest(betadisper(pcyc.bray, design$Date), permutations = how(nperm = 9999))



# Canonical analysis of principal coordinates (CAP) ============

# CAZyme - irrigation and cropping system

nc <- nrow(as.matrix(cazy.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(cazy.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

cap.cazy <- CAPdiscrim(cazy.bray ~ IxT, data = design, m = 27, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(cap.cazy$group), data.frame(cap.cazy$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(cap.cazy$PCoA)
success <- success[order(success$source), ]
success

cap1.cazy <- paste("CAP1 (", round((100/sum(cap.cazy$lda.other$svd^2) * cap.cazy$lda.other$svd^2)[1],
                                   digits = 1), "%)", sep = "")
cap2.cazy <- paste("CAP2 (", round((100/sum(cap.cazy$lda.other$svd^2) * cap.cazy$lda.other$svd^2)[2],
                                   digits = 1), "%)", sep = "")


# NCycDB - irrigation and cropping system

nc <- nrow(as.matrix(ncyc.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(ncyc.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

cap.ncyc <- CAPdiscrim(ncyc.bray ~ IxT, data = design, m = 13, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(cap.ncyc$group), data.frame(cap.ncyc$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(cap.ncyc$PCoA)
success <- success[order(success$source), ]
success

cap1.ncyc <- paste("CAP1 (", round((100/sum(cap.ncyc$lda.other$svd^2) * cap.ncyc$lda.other$svd^2)[1],
                                   digits = 1), "%)", sep = "")
cap2.ncyc <- paste("CAP2 (", round((100/sum(cap.ncyc$lda.other$svd^2) * cap.ncyc$lda.other$svd^2)[2],
                                   digits = 1), "%)", sep = "")



# PCycDB - irrigation and cropping system

nc <- nrow(as.matrix(pcyc.bray))
success <- data.frame(m = numeric(nc), class.success = numeric(nc))
for (i in 1:60) {
  cap <- CAPdiscrim(pcyc.bray ~ IxT, data = design, m = i, add = TRUE)
  success[i, 1] <- cap$m
  success[i, 2] <- 100/length(cap$group) * length(which(cap$group == cap$CV))
}

par(mfrow = c(1,1))

plot(success$m, success$class.success, xlab = "Number of PCoA axes included", ylab = "Reclassification success rate (%)")
text(success$m, success$class.success, labels = success$m, pos = 1, cex = 0.6)

cap.pcyc <- CAPdiscrim(pcyc.bray ~ IxT, data = design, m = 21, permutations = 9999, add = TRUE) #try with different

success <- cbind(data.frame(cap.pcyc$group), data.frame(cap.pcyc$CV))
colnames(success) <- c("source", "classified")
rownames(success) <- rownames(cap.pcyc$PCoA)
success <- success[order(success$source), ]
success

cap1.pcyc <- paste("CAP1 (", round((100/sum(cap.pcyc$lda.other$svd^2) * cap.pcyc$lda.other$svd^2)[1],
                                   digits = 1), "%)", sep = "")
cap2.pcyc <- paste("CAP2 (", round((100/sum(cap.pcyc$lda.other$svd^2) * cap.pcyc$lda.other$svd^2)[2],
                                   digits = 1), "%)", sep = "")



# Visualization of CAP ============

# CAZyme

ggplot(as.data.frame(cap.cazy$x), aes(x = cap.cazy$x[,1], y = cap.cazy$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(cap1.cazy) + ylab(cap2.cazy) +
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
  
  annotate("text", x= -7, y= -4.5,label= "Overall reclassification rate: 81%", hjust = 0, size = 4) +
  annotate("text", x= -7, y= -5,label= "Pillais test=3.7 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -5, y= -3,label= "92%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -5, y=  4,label= "83%", hjust = 0, size = 4, color = "#009E73") +
  
  annotate("text", x=  2.5, y= -1.5, label= "83%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x=  5, y= -4,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x=  2, y= 1,label= "67%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= -1.5, y= 5,label= "67%", hjust = 0, size = 4, color = "#E69F00") 


# NCycDB

ggplot(as.data.frame(cap.ncyc$x), aes(x = cap.ncyc$x[,1], y = cap.ncyc$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(cap1.ncyc) + ylab(cap2.ncyc) +
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
  
  annotate("text", x= -5, y= -3.5,label= "Overall reclassification rate: 74%", hjust = 0, size = 4) +
  annotate("text", x= -5, y= -4,label= "Pillais test=2.7 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -5.5, y= -2.5,label= "100%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -4, y=  3,label= "83%", hjust = 0, size = 4, color = "#009E73") +
  
  annotate("text", x= 0.2, y= -3, label= "75%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x= 4.5, y= -3,label= "50%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= -0.5, y= -0.5,label= "75%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x=  4, y= 4,label= "58%", hjust = 0, size = 4, color = "#E69F00") 


# PCycDB

ggplot(as.data.frame(cap.pcyc$x), aes(x = cap.pcyc$x[,1], y = cap.pcyc$x[,2])) +
  geom_point(aes(shape = design$Irrigation, color = design$Treatment), size = 2) +
  xlab(cap1.pcyc) + ylab(cap2.pcyc) +
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
  
  annotate("text", x= -8, y= -4.5,label= "Overall reclassification rate: 79%", hjust = 0, size = 4) +
  annotate("text", x= -8, y= -5,label= "Pillais test=3.1 ***", hjust = 0, size = 4) +
  
  annotate("text", x= -4.5, y=  3,label= "92%", hjust = 0, size = 4, color = "#009E73") +
  annotate("text", x= -7.5, y= -2,label= "92%", hjust = 0, size = 4, color = "#009E73") +
  
  annotate("text", x= -1, y= -0, label= "67%", hjust = 0, size = 4, color = "#FF618C") +
  annotate("text", x=  5, y= -4,label= "92%", hjust = 0, size = 4, color = "#FF618C") +
  
  annotate("text", x= 1, y= 3,label= "67%", hjust = 0, size = 4, color = "#E69F00") +
  annotate("text", x= 2, y= 5,label= "67%", hjust = 0, size = 4, color = "#E69F00") 


