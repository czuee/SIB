setwd("~/Documents/SIB-internship")
load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-RNASeq2GeneNorm.Rdata")
load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-Clinical.Rdata")

install.packages("tidyr")
library(tidyr)
library(dplyr)

###Check for the duplicates #### to see if there are any normal tissues

X <- RNASeq2GeneNorm$Patient_Id

duplicates <- duplicated(X) | duplicated(X, fromLast = TRUE)
table(duplicates) 
#duplicates
#FALSE  TRUE 
#476    90 

rows <- data.frame(dupli = rownames(RNASeq2GeneNorm)[duplicates])

tcga.code <- tidyr::separate(rows, dupli, c("Project", "Tissue source", "Participant", "Sample/Vial", "Portion/Analyte", "Plate", "Seq center"))

#The duplicates are in the samples as expected
tcga.code <- data.frame(tcga.code)
table(tcga.code$Sample.Vial) # 01 = Primary solid tumor, 06 = Metastatic, 11 = Solid tissue normal
#01A 06A 11A 
#45   2  43

#Since the gene expression data has normal tissue expression, these need to be removed. All normal tissues have Sample_Id =  11. Use this to remove the normal samples. The rest of the duplicates (Metastatic) can stay.

data <- RNASeq2GeneNorm
dat.normal <- data[data$Sample_Id == "11", ]
dat.tumor <- data[data$Sample_Id != "11", ]
nrow(dat.normal) #44 : There is also ONE normal tissue sample without any cancer tissue

dupli.bin <- which(duplicates)
normal <- grep("11", RNASeq2GeneNorm$Sample_Id)
which(!(normal %in% dupli.bin)) %>% normal[.] %>% data[., 1]
data[data$Patient_Id == "TCGA-H7-A6C5", 1:3] 

#X -> filter for SD -> set gene-mean=0 -> double clustering (gene/sample) -> heatmap /gene correlation. DESEQ/LIMMA for differential expression (another progam)

#PLot SD of the genes to filter low SD genes
colnames(dat.tumor)[-(1:2)] <- paste("X", colnames(dat.tumor)[-(1:2)], sep = "")
dat.log <- cbind(dat.tumor$Patient_Id, log(1+ dat.tumor[ ,-(1:2)])) # Data with log(1+e)

sd.g <- sapply( dat.log[ ,-1], FUN = sd )
sd.g <- sort(sd.g, decreasing = TRUE)
plot(sd.g, main = "SD distribution log(1+e)", xlab = "Genes", ylab = "Std Dev")
sd.h <- hist(sd.g, xlab = "SD (gene expression)", main = "SD values log(1+e)")

#Draw a normal curve on the histogram
# xfit<-seq(min(sd.g),max(sd.g),length=40)
# yfit<-dnorm(xfit,mean=mean(sd.g),sd=sd(sd.g))
# yfit <- yfit*diff(sd.h$mids[1:2])*length(sd.g)
# lines(xfit, yfit, col="red", lwd=2)

###Testing for normality
# sd.norm <- rnorm(1000, mean=mean(sd.g), sd=sd(sd.g))
# plot(density(sd.g))
# plot(density(sd.norm))
# shapiro.test(sample(sd.g, size = 5000))
# shapiro.test(sd.norm)
# 
# qqnorm(sd.g);qqline(sd.g, col = 2) #NOT a normal distr
# qqnorm(sd.norm);qqline(sd.norm, col = 2)

#Convert to decreasing and plot the top 1000/2000 cutoff line
plot(sd.g)
abline(v = 2000, col = 2)
text(2000, 0.5, "top 2000", col = 2, pos = 4)
abline(v = 1000, col = "blue")
text(1000, 2.5, "top 1000", col = "blue", pos = 4)

sd1000 <- sd.g[1:1000] #Consider top 5% of genes with highest variance
sd2000 <- sd.g[1:2000] #Consider top 10% of genes with highest variance

# Histogram overlapping shows that the SD values are not very different in the two sets
# hist(sd2000, col=rgb(0.1,0.1,0.1,0.5), main="Overlapping Histogram")
# hist(sd1000, col=rgb(0.8,0.8,0.8,0.5), add=T)
# box()

sd1000 <- data.frame(t(sd.g[1:1000]))
sd2000 <- data.frame(t(sd.g[1:2000]))

g1000 <- dat.log[ ,match(colnames(sd1000), colnames(dat.log))]
g2000 <- dat.log[ ,match(colnames(sd2000), colnames(dat.log))]

### Set expression mean = 0 for each gene ###

g1000.mean0 <- sapply(g1000, function(x) (x - mean(x))) #Subtract mean of each column from each variable
g2000.mean0 <- sapply(g2000, function(x) (x - mean(x))) #Subtract mean of each column from each variable

hist(g2000.mean0)

# Integrate the gender data with the expression data
g1k <- cbind(dat.tumor[ ,1] ,data.frame(g1000.mean0))
g2k <- cbind(dat.tumor[ ,1] ,data.frame(g2000.mean0))
colnames(g1k)[1] <- "Patient_Id"; colnames(g2k)[1] <- "Patient_Id"

gender <- subset(ClinicalData, select = c(Patient_Id, gender))
g1k.mf <- merge(gender, g1k, by = "Patient_Id")
g2k.mf <- merge(gender, g2k, by = "Patient_Id")

### Controls ###
test_expr <- g1000.mean0[ ,1:10] # To see the spread of genes with highest SD
dat.m <- melt(test_expr, id.vars = colnames(test_expr))
ggplot(dat.m) + geom_boxplot(aes(x=Var2, y=value)) + coord_flip() + ggtitle("10 highest SD genes") + ylab("log 1+ gene expression")

x = g1k.mf[ ,1:22]
dat <- melt(x, id.vars = 'gender', measure.vars = colnames(x[-(1:2)]))
ggplot(dat) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for top 20 genes") + xlab("Genes with highest variance") + ylab("log expression, mean 0")

GeneId[grep("7503", GeneId$Entrez_Gene_Id, fixed = TRUE), ] # XIST X-gene inactivation
GeneId[grep("6192", GeneId$Entrez_Gene_Id, fixed = TRUE), ] # RPS4Y1 ribosomal protein
GeneId[grep("8653", GeneId$Entrez_Gene_Id, fixed = TRUE), ] # DDX3Y Y-linked DEAD-box RNA helicase

y = cbind(g1k.mf[ ,1:2], g1k.mf[, 23:50])
dat <- melt(y, id.vars = 'gender', measure.vars = colnames(y[-(1:2)]))
ggplot(dat) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for top 20 genes") + xlab("Genes with highest variance") + ylab("log expression, mean 0")

#Which genes have a significant difference in their means? Abs Mean diff > SD of either M or F

#Kurtosis calculation 
install.packages("e1071")
library(e1071)

g2k.kurt <- sapply(g2k.mf[ ,-(1:2)], kurtosis)
g2k.kurt.m <- sapply(g2k.m[ ,-(1:2)], kurtosis)
g2k.kurt.f <- sapply(g2k.f[ ,-(1:2)], kurtosis)

g2k.k[1250 ,]
g2k.k <- data.frame(total=g2k.kurt, male=g2k.kurt.m, female=g2k.kurt.f)
g2k.k[, 4] = (g2k.k$total > max(g2k.k$male, g2k.k$female))

#The kurtosis might not be a good measure because of the outliers.

g2k.f <- filter(g2k.mf, gender == "female")
g2k.m <- filter(g2k.mf, gender == "male")

mean.m <- sapply(g2k.m[ ,-(1:2)], mean)
sd.m <- sapply(g2k.m[ ,-(1:2)], sd)
mean.f <- sapply(g2k.f[ ,-(1:2)], mean)
sd.f <- sapply(g2k.f[ ,-(1:2)], sd)
mean.mf = data.frame(meanM=mean.m, sdM=sd.m, meanF=mean.f, sdF=sd.f)

#Abs Mean diff > SD of either M or F
mean.mf[ ,5] = abs(mean.mf$meanM - mean.mf$meanF) > mean.mf$sdM & abs(mean.mf$meanM - mean.mf$meanF) > mean.mf$sdF
colnames(mean.mf)[5] = "signif"

# There are 16 genes in 2000 that have a significant difference in means
plot(meanM ~ meanF, mean.mf, col=(signif + 1), main = "Mean differences (16 signif genes)", xlab = "Mean for F", ylab = "mean for M")

Genes.signif <- rownames(mean.mf)[mean.mf$signif == TRUE]

### End Controls ###

#Clustering
#The dist function calculates distances between rows of a matrix. I'm using Euclidean distance here since the data is log tranformed.
g1k.mf.m <- data.matrix(g1k.mf[ ,-(1:2)])
g1k.dist.sam <- dist(g1k.mf.m)

clust.sam <- hclust(g1k.dist.sam, method = "ward.D")
plot(clust.sam)
groups.sam <- cutree(clust.sam, 6)
rect.hclust(clust.sam, k=6, border="red")
table(groups.sam)
#  1   2   3   4   5   6 
#103  47 123  79  69 101 
groups.gender <- cbind(group=groups.sam, gender=g1k.mf$gender)
plot(gender ~ group, groups.gender)

clust.sam1 <- hclust(g1k.dist.sam) #method =complete
plot(clust.sam1)
groups.sam1 <- cutree(clust.sam1, 6)
rect.hclust(clust.sam1, k=6, border="red")
table(groups.sam1)
#  1   2   3   4   5   6 
# 52  49 103  55 145 118 
groups.gender1 <- cbind(group=groups.sam1, gender=g1k.mf$gender)
plot(gender ~ group, groups.gender1)

#The clustering does not show differences between M and F gene expression with both the methods

#Lets look at the gene clustering
g1k.mf.t <- t(g1k.mf.m) #Transpose the matrix
g1k.dist <- dist(g1k.mf.t)

clust.gene <- hclust(g1k.dist, method = "ward.D")
plot(clust.gene)
groups.gene <- cutree(clust.gene, 4)
gene.cluster <- rect.hclust(clust.gene, k=4, border="red")
sapply(gene.cluster, length) #The actual order in the dendrogram
# 134 259 400 207
table(groups.gene)
#   1   2   3   4 
# 134 400 207 259 

clust.gene1 <- hclust(g1k.dist, method = "complete")
plot(clust.gene1)
groups.gene1 <- cutree(clust.gene1, 5)
gene.cluster1 <- rect.hclust(clust.gene1, k=5, border="red")
sapply(gene.cluster1, length) #The actual order in the dendrogram
# 13 522 238 133  94
table(groups.gene1)
#   1   2   3   4   5 
# 133  94  13 238 522 

library(RColorBrewer)
display.brewer.all() # To see all availabel colors
niceCols <- brewer.pal(9, "Set1")

#Plot the distribution for the clustered genes from both the sets
#Matplot plots the columns of one matrix against columns of another
oPar <- par(mfrow=c(2,2))
for (i in 1:4) {
  matplot(t(g1k.mf.t[groups.gene == i,]),
          type="p", col=niceCols[i], 
          xlab="sample",ylab="log expression value")
}
par(oPar)

oPar <- par(mfrow=c(3,2))
for (i in 1:5) {
  matplot(t(g1k.mf.t[groups.gene1 == i,]),
          type="l", col=niceCols[i], 
          xlab="sample",ylab="log expression value")
}
par(oPar)

#Does the group 1 from both clustering methods have the same genes?
match.gene <- (colnames(g1k.mf.m[ ,groups.gene == 1]) %in% colnames(g1k.mf.m[ ,groups.gene1 == 1]))
which(match.gene == FALSE)
#All the genes match except #130 that is extra in groups.gene
colnames(g1k.mf.m[ ,groups.gene == 1])[130] #X6444

