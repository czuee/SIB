#Read in clinical data and expression data and plot the expression of genes in males and females as a box plot and heatplot

#Extract patient id and gender data from the Clinical data
setwd("~/Documents/SIB-internship")
load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-Clinical.Rdata")

str(ClinicalData) #Check the data
View(ClinicalData)

gender1 <- subset(ClinicalData, select = c(Patient_Id, gender))

#Load the expression data and match the patient ids to extract the gender data
load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-RNASeq2GeneNorm.Rdata")
RNASeq2Gene_g <- merge(gender1, RNASeq2GeneNorm, by="Patient_Id")

#Save colnames as x+Entrez_id IF NECESSARY
#colnames(RNASeq2GeneNorm)[-(1:2)] = paste("x", colnames(RNASeq2GeneNorm)[-(1:2)], sep="")
#colnames(RNASeq2Gene_g)[-(1:3)] = paste("x", colnames(RNASeq2Gene_g)[-(1:3)], sep="")

table(RNASeq2Gene_g$gender) #The number of male & female entries
#female   male 
#   151    415 

#Convert expression data to log scale
log_RNASeq2Gene <- cbind(RNASeq2Gene_g[,1:3], log(1+(RNASeq2Gene_g[,-(1:3)]))) #log(1+e)

### CONTROLS ###
library(reshape2)
library(ggplot2)

#To check the expression levels for different genes
test_expr <- log_RNASeq2Gene[,1:13]

dat.m <- melt(test_expr, id.vars = 'gender', measure.vars = colnames(test_expr)[-(1:3)])
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for first 10 genes") + ylab("log 1+gene expression")

set.seed(1)
test_genes2 <-  cbind(log_RNASeq2Gene[,1:3], sample(log_RNASeq2Gene[ ,-(1:3)], 10))

dat.m2 = melt(test_genes2, id.vars = 'gender', measure.vars = colnames(test_genes2[-(1:3)]))
ggplot(dat.m2) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for random 10 genes") + ylab("log gene expression")

#To ensure that the expression levels for all genes are in the same range for different patients
test_pat1 <- log_RNASeq2Gene[1:10,] #First 10 patients

dat.p <- melt(test_pat1, id.vars = "Patient_Id", measure.vars = colnames(test_pat1)[-(1:3)])
ggplot(dat.p, aes(Patient_Id, value)) + geom_boxplot() + coord_flip() + ggtitle("Gene expression for first 10 individuals") + ylab("log 1+ gene expression")

set.seed(1)
test_pat2 <- log_RNASeq2Gene[sample(nrow(log_RNASeq2Gene), 10), ]
dat.p2 <- melt(test_pat2, id.vars = "Patient_Id", measure.vars = colnames(test_pat2)[-(1:3)])
ggplot(dat.p2, aes(Patient_Id, value)) + geom_boxplot() + coord_flip() + ggtitle("Gene expression for 10 random individuals") + ylab("log 1+gene expression")

### Test genes from papers ###

#Genes that are significantly mutated in HNSC cancers; Nature 2015 517:576 TCGA paper
mutsig.d <- data.frame(c("CDKN2A", "FAT1","TP53","CASP8", "AJUBA", "PIK3CA","NOTCH1","KMT2D","NSD1", "HLA-A","TGFBR2","HRAS", "FBXW7","RB1","PIK3R1", "TRAF3", "NFE2L2","CUL3", "PTEN"))
names(mutsig.d) <- "Protein_Coding"

load("~/Documents/SIB-internship/Data/reference/GeneId.Rdata")

mutsig <- merge(mutsig.d, GeneId, by="Protein_Coding", all.x = TRUE) #Match for entrez gene id and use NA where no matching id is found

## The genes for which no entrez ids were found were AJUBA and KMT2D. A search in the HGNC data set returns that the entrez gene id for these are present in our data set (AJUBA = 84962/JUB; KMT2D = 8085/MLL2) but with different names that are probably the earlier names used for these genes.

mutsig$Entrez_Gene_Id[(mutsig$Protein_Coding == "AJUBA")] = "84962" #Find a direct way to do this through Bioconductor
mutsig$Entrez_Gene_Id[(mutsig$Protein_Coding == "KMT2D")] = "8085"

mutsig_genes <- data.frame(t(mutsig)) #Transpose the matrix
colnames(mutsig_genes) <- as.vector(mutsig$Entrez_Gene_Id)

#Match the column names in the two data frames to get the expression values for the mutated genes. I also changed the names of the genes to the protein-coding names mentioned in the paper.
mutsig_genes1 <- cbind(log_RNASeq2Gene[ ,1:3], log_RNASeq2Gene[ , match(colnames(mutsig_genes), colnames(log_RNASeq2Gene))])
colnames(mutsig_genes1)[-(1:3)] <- as.vector(mutsig$Protein_Coding)

dat.g = melt(mutsig_genes1, id.vars = 'gender', measure.vars = colnames(mutsig_genes1[-(1:3)]))
ggplot(dat.g) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Significantly mutated genes(MutSig)") + ylab("log gene expression") + xlab("Gene id")

#Genes from PLOS1 2013 paper implicated in HNSC cancers
hnsc.d <- data.frame(c("RPA2", "E2F2", "FGFR3", "SOX2", "NFE2L2", "KEAP1", "PIK3CA", "AKR1C1", "PDGFRA", "PDGFRB", "TWIST1", "TP63", "TGFA", "EGFR"))
names(hnsc.d) <- "Protein_Coding"

hnsc.g <- merge(hnsc.d, GeneId, by="Protein_Coding", all.x = TRUE)
hnsc_gene <- data.frame(t(hnsc.g))
colnames(hnsc_gene) <- as.vector(hnsc.g$Entrez_Gene_Id)

hnsc_gene1 <- cbind(log_RNASeq2Gene[,1:3], log_RNASeq2Gene[ , match(colnames(hnsc_gene), colnames(RNASeq2Gene_g))])
colnames(hnsc_gene1)[-(1:3)] <- as.vector(hnsc.g$Protein_Coding)

dat.g2 = melt(hnsc_gene1, id.vars = 'gender', measure.vars = colnames(hnsc_gene1[-(1:3)]))
ggplot(dat.g2) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Genes assoc with HNSC, PlosOne 2013") + ylab("log gene expression") + xlab("Gene id")

#### Heat maps ###

RNASeq2Gene_mf <- log_RNASeq2Gene[order(log_RNASeq2Gene$gender),] # To separate males and females

#The data needs to be in a matrix format with only numbers. The row headers also need to be saved in another vector
RNASeq2Gene_mf.m <- data.matrix(RNASeq2Gene_mf[-(1:3)]) #Convert to matrix

distr <- as.vector(RNASeq2Gene_mf.m[1:151,]) #Female expression values
hist(distr, main = "Gene expression values for females")
mean(distr <= 1) #20.8%
mean(distr == 0) #14.5%
mean(distr) #4.42
median(distr) #5.16
quantile(distr, c(0.25, 0.75)) # 1.693687 6.710150  
max(distr) #14.88346

distr2 <- as.vector(RNASeq2Gene_mf.m)
hist(distr2, main = "Gene expression values for all")

dist3 <- as.vector(RNASeq2Gene_mf.m[152:nrow(RNASeq2Gene_mf.m), ]) #Male expression values
hist(dist3, main = "Gene expression values for males")

#Plotting the heatmap using heatmap.2 from gplots
install.packages("gplots")
library(gplots)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

fem <- RNASeq2Gene_mf.m[1:151,]
mal <- RNASeq2Gene_mf.m[152:nrow(RNASeq2Gene_mf.m),]
set.seed(1)
test_mf <- rbind(fem[sample(nrow(fem), 100), ], mal[sample(nrow(mal), 100), ]) #100 individuals

mf_heatmap <- heatmap.2(test_mf[ , 100:500], trace="none", Rowv = NA, col = my_palette, margins=c(9,12), labRow = NA, labCol = NA, dendrogram = "column")
