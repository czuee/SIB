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
### End Controls ###

#Which genes have a significant difference in their means? Abs Mean diff > SD of entire distribution

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

#Which genes have a significant difference in their means? Abs Mean diff > SD of entire distribution

plot_signif <- function(x){ #Returns a plot of mean M vs mean F & vector of "signif" genes
  x.f <- dplyr::filter(x, gender == "female")
  x.m <- dplyr::filter(x, gender == "male")
  
  mean.m <- sapply(x.m[ ,-(1:2)], mean)
  sd.all <- sapply(x[ ,-(1:2)], sd)
  mean.f <- sapply(x.f[ ,-(1:2)], mean)
  mean.mf = data.frame(meanM=mean.m, meanF=mean.f, sdAll=sd.all)
  
  #Abs Mean diff > SD of entire distrib
  mean.mf[ ,4] = abs(mean.mf$meanM - mean.mf$meanF) > mean.mf$sdAll
  
  colnames(mean.mf)[4] = c("signif")
  num <- sum(mean.mf$signif == TRUE) #Number of signif genes = 16
  genes.signif <- rownames(mean.mf)[mean.mf$signif == TRUE]
  write(genes.signif, file = "genes.signif.Rdata")
  plot(meanM ~ meanF, mean.mf, col=(signif + 1), main = paste("Mean differences:", num, "signif genes"), xlab = "Mean for F", ylab = "mean for M")
}

plot_signif(g2k.mf) #16 signif genes

dat.log <- cbind(dat.tumor[ ,1:2], log(1+ dat.tumor[ ,-(1:2)])) #log(1+e)
colnames(dat.log)[1] = "Patient_Id"
dat.log.mean0 <- data.frame(sapply(dat.log[ ,-(1:2)], function(x) (x-mean(x)))) #Mean centering
dat.mean0 <- cbind(dat.log[ ,1], dat.log.mean0)
colnames(dat.mean0)[1] = "Patient_Id"

genes.all <- merge(gender, dat.log, by = "Patient_Id")
genes.mean0 <- merge(gender, dat.mean0, by = "Patient_Id")

plot_signif(genes.all) #21 signif genes
plot_signif(genes.mean0)

# g2k.f <- dplyr::filter(g2k.mf, gender == "female")
# g2k.m <- dplyr::filter(g2k.mf, gender == "male")
# 
# mean.m <- sapply(g2k.m[ ,-(1:2)], mean)
# sd.all <- sapply(g2k.mf[ ,-(1:2)], sd)
# mean.f <- sapply(g2k.f[ ,-(1:2)], mean)
# mean.mf = data.frame(meanM=mean.m, meanF=mean.f, sdAll=sd.all)
# 
# #Abs Mean diff > SD of entire distrib
# mean.mf[ ,4] = abs(mean.mf$meanM - mean.mf$meanF) > mean.mf$sdAll
# 
# colnames(mean.mf)[4] = c("signif")
# sum(mean.mf$signif == TRUE) #Number of signif genes = 16
# 
# Genes.signif.2k <- rownames(mean.mf)[mean.mf$signif == TRUE]
# # There are 16 genes in 2000 that have a significant difference in means
# plot(meanM ~ meanF, mean.mf, col=(signif + 1), main = "Mean differences (16 signif genes)", xlab = "Mean for F", ylab = "mean for M")

plot(abs(meanM - meanF) ~ sdAll, mean.mf, col = signif+1, main = "Significant genes selection:top16 genes", xlab = "SD of distribution", ylab = "Diff in means btw M and F")

#Find signif genes in entire distribution [ALL]
# genes.all.m <- dplyr::filter(genes.all, gender == "male")
# genes.all.f <- dplyr::filter(genes.all, gender == "female")
# 
# all.mean.m <- sapply(genes.all.m[ ,-(1:2)], mean)
# all.sd.all <- sapply(genes.all[ ,-(1:2)], sd)
# all.mean.f <- sapply(genes.all.f[ ,-(1:2)], mean)
# all.mean.mf = data.frame(meanM=all.mean.m, meanF=all.mean.f, sdAll=all.sd.all)
# 
# #Abs Mean diff > SD of entire distrib
# all.mean.mf[ ,4] = abs(all.mean.mf$meanM - all.mean.mf$meanF) > all.mean.mf$sdAll
# 
# colnames(all.mean.mf)[4] = c("signif")
# sum(all.mean.mf$signif == TRUE) #Number of signif genes
# # 21 (5 more genes found!)
# 
Genes.signif <- rownames(all.mean.mf)[all.mean.mf$signif == TRUE]
# # There are 21 genes from 20531 that have a significant difference in means
# plot(meanM ~ meanF, all.mean.mf, col = (signif +1), main = "Mean differences (21 signif genes)", xlab = "Mean for F", ylab = "mean for M")

plot(abs(meanM - meanF) ~ sdAll, all.mean.mf, col = signif+1, main = "Significant genes selection:21 genes", xlab = "SD of distribution", ylab = "Diff in means btw M and F")

# Get Protein coding names for the signif genes
signif21 <- data.frame(Genes.signif)
colnames(signif21)[1] = "Entrez_Gene_Id"
GeneId$Entrez_Gene_Id <- paste("X", GeneId$Entrez_Gene_Id, sep = "")
genes.signif21 <- merge(signif21, GeneId, by = "Entrez_Gene_Id", all.x = TRUE) ### Signif gene names

signif21.exprn = data.frame(t(signif21))
colnames(signif21.exprn) <- genes.signif21$Entrez_Gene_Id
signif21.exprn <- genes.mean0[ ,match(colnames(signif21.exprn), colnames(genes.mean0))]
colnames(signif21.exprn) <- genes.signif21$Protein_Coding
signif21.exprn <- cbind(genes.mean0[ ,1:2], signif21.exprn)

library(reshape2)
library(ggplot2)
dat.g = melt(signif21.exprn, id.vars = 'gender', measure.vars = colnames(signif21.exprn[-(1:2)]))
ggplot(dat.g) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Significant diff btw M & F") + ylab("log 1+gene expression") + xlab("Gene id")

#Getting gene names from HGNC database - Hugo Gene Nomenclature committee. First download the genes required from the download tool at http://www.genenames.org/cgi-bin/download. Here I downloaded the genes for X and Y chromosomes.
genes_XY <- read.table("Data/genedata_XY.txt", header = TRUE, sep = "\t", na.strings = c("", " ", "NA", "NaN", "<NA>"))
names(genes_XY)

#Reading in all the gene names
gene.names <- read.table("Data/genenames_all.txt", header = TRUE, sep = "\t", na.strings = c("", " ", "NA", "NaN", "<NA>"), skip = 29000, nrows = 10000)
classes <- sapply(gene.names, class) #To set column classes in order to speed up data reading
gene.names <- read.table("Data/genenames_all.txt", 
                         header = TRUE, 
                         sep = "\t", 
                         na.strings = c("", " ", "NA", "NaN", "<NA>"), 
                         colClasses = classes, 
                         quote = "")

gene.names$Entrez.Gene.ID <- paste("X", gene.names$Entrez.Gene.ID, sep = "")
colnames(genes.signif21) = c("Entrez.Gene.Id" , "Approved.Symbol")
genenames.signif21 <- merge(genes.signif21, gene.names, by = "Approved.Symbol")
genenames.signif21 <- genenames.signif21[ , -(2)]

prev_symb <- gene.names[unique(grep(paste(genes.signif21$Approved.Symbol, collapse="|"), gene.names$Previous.Symbols)), ]
prev_symb <- rbind(prev_symb[1, ], prev_symb[7:8, ])
genenames.signif21.all <- unique(rbind(genenames.signif21, prev_symb))
nrow(genenames.signif21.all)
#19 descriptions - NCRNA00185 and TTTY14 are the same gene. Also CYorf15A and B have the same description. ALl the genes are X and Y linked only.