#Read in clinical data and expression data and plot the expression of genes in males and females as a box plot and heatplot

setwd("~/Documents/SIB-internship")
load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-Clinical.Rdata")

str(ClinicalData) #Check the data
View(ClinicalData)

#Extract patient id and gender data
gender1 <- subset(ClinicalData, select = c(Patient_Id, gender))

load("~/Documents/SIB-internship/Data/HNSC/20151101-HNSC-RNASeq2GeneNorm.Rdata")

#Get gender values for entire dataset
RNASeq2Gene_g <- merge(gender1, RNASeq2GeneNorm, by="Patient_Id")

#Save colnames as x+Entrez_id
colnames(RNASeq2GeneNorm)[-(1:2)] = paste("x", colnames(RNASeq2GeneNorm)[-(1:2)], sep="")
colnames(RNASeq2Gene_g)[-(1:3)] = paste("x", colnames(RNASeq2Gene_g)[-(1:3)], sep="")

min(RNASeq2GeneNorm[-(1:2)]); max(RNASeq2GeneNorm[-(1:2)])
#[1] 0
#[1] 6925542

table(RNASeq2Gene_g$gender)
#female   male 
#   151    415 

############### TESTING ###################
#Testing on a subset of the data
test_g <- head(gender1, 20)
test_expr <- RNASeq2GeneNorm[1:20,1:20]
test_tot <- merge(test_g, test_expr, by="Patient_Id") # To get gender in the beginning
test_tot[,1:7]

load("~/Documents/SIB-internship/Data/reference/GeneId.Rdata")

#Work on the entire dataset for few genes
test_expr2 <- RNASeq2GeneNorm[,1:7]
test_tot2 <- merge(gender1, test_expr2, by="Patient_Id") # Also merges duplicates from expression data


library(reshape2)
dat.m <- melt(test_tot2, id.vars = 'gender', measure.vars = c('x100130426', 'x100133144', 'x100134869', 'x10357', 'x10431'))

ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, color=gender))


# Test set with value of all genes for 10 patients. No gender here! See distribution of gene expression across individuals
test_pat1 <- RNASeq2GeneNorm[1:10, ] #First 10 patients

log_test_pat1 <- cbind(test_pat1[,1], 1+log(1+(test_pat1[,-(1:2)]))) #### Take log of the gene expression values
names(log_test_pat1)[1] = "Patient_Id"

dat.p <- melt(log_test_pat1, id.vars = "Patient_Id", measure.vars = colnames(log_test_pat1)[-(1)])
ggplot(dat.p, aes(Patient_Id, value)) + geom_boxplot() + coord_flip() + ggtitle("Gene expression for 10 individuals") + ylab("1+log gene expression")

test_pat2 <- RNASeq2GeneNorm[sample(nrow(RNASeq2GeneNorm), 10), ] #10 random individuals

log_test_pat2 <- cbind(test_pat2[,1], 1+log(1+(test_pat2[,-(1:2)]))) #### Take log of the gene expression values
names(log_test_pat2)[1] = "Patient_Id"

dat.p2 <- melt(log_test_pat2, id.vars = "Patient_Id", measure.vars = colnames(log_test_pat2)[-(1)])
ggplot(dat.p2, aes(Patient_Id, value)) + geom_boxplot() + coord_flip() + ggtitle("Gene expression for 10 random individuals") + ylab("1+log gene expression")

################# Testing over ##########

library(reshape2)
library(ggplot2)

plot_mf <- function(data, title){
  log_data <- cbind(data[,1:2], log(1+(data[,-(1:3)]))) #log(1+e)
  dat.g = melt(log_data, id.vars = 'gender', measure.vars = colnames(log_data[-(1:2)]))
  ggplot(dat.g) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle(title) + ylab("log gene expression") + xlab("Gene id")
}

#Boxplot for first 10 genes
test_genes1 <- RNASeq2Gene_g[, 1:13]
plot_mf(test_genes1, "Gene expression for first 10 genes")

# log_test_genes1 <- cbind(test_genes1[,1:2], log(1+(test_genes1[,-(1:3)]))) #log(1+e) 
# 
# library(reshape2)
# dat.g = melt(log_test_genes1, id.vars = 'gender', measure.vars = colnames(log_test_genes1[-(1:2)]))
# library(ggplot2)
# ggplot(dat.g) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for first 10 genes") + ylab("log gene expression")

#boxplot for random 10 genes
set.seed(1)
test_genes2 <- cbind(RNASeq2Gene_g[ ,1:3] ,sample(RNASeq2Gene_g[, -(1:3)], 10))
plot_mf(test_genes2, "Gene expression for random 10 genes")

# log_test_genes2 <- cbind(test_genes2[,1:2], log(1+(test_genes2[,-(1:3)])))
# 
# dat.g2 = melt(log_test_genes2, id.vars = 'gender', measure.vars = colnames(log_test_genes2[-(1:2)]))
# library(ggplot2)
# ggplot(dat.g2) + geom_boxplot(aes(x=variable, y=value, color=gender)) + coord_flip() + ggtitle("Gene expression for random 10 genes") + ylab("log gene expression")

##Significantly mutated genes mentioned in Nature2015 paper as analyzed by the MutSigCv algorithm. The analysis could be repeated with our dataset to get a better picture.

mutsig.d <- data.frame(c("CDKN2A", "FAT1","TP53","CASP8", "AJUBA", "PIK3CA","NOTCH1","KMT2D","NSD1", "HLA-A","TGFBR2","HRAS", "FBXW7","RB1","PIK3R1", "TRAF3", "NFE2L2","CUL3", "PTEN"))
names(mutsig.d) <- "Protein_Coding"

load("~/Documents/SIB-internship/Data/reference/GeneId.Rdata")

mutsig <- merge(mutsig.d, GeneId, by="Protein_Coding", all.x = TRUE) #Match for entrez gene id and use NA where no matching id is found

## The genes for which no entrez ids were found were AJUBA and KMT2D. A search in the HGNC data set returns that the entrez gene id for these are present in our data set (AJUBA = 84962/JUB; KMT2D = 8085/MLL2) but with different names that are probably the earlier names used for these genes.

mutsig$Entrez_Gene_Id[(mutsig$Protein_Coding == "AJUBA")] = "84962" #Find a direct way to do this through Bioconductor
mutsig$Entrez_Gene_Id[(mutsig$Protein_Coding == "KMT2D")] = "8085"

mutsig_genes <- as.data.frame(t(mutsig))
colnames(mutsig_genes) <- as.vector(mutsig$Entrez_Gene_Id)

#Match the column names in the two data frames to get the expression values for the mutated genes. I also changed the names of the genes to the protein-coding names mentioned in the paper.

mutsig_genes1 <- cbind(RNASeq2Gene_g[,1:3], RNASeq2Gene_g[ , match(colnames(mutsig_genes), colnames(RNASeq2Gene_g))])
colnames(mutsig_genes1)[-(1:3)] <- as.vector(mutsig$Protein_Coding)


plot_mf(mutsig_genes1, "Significantly mutated genes(MutSig)")

#Heatmap
mutsig_mf <- mutsig_genes1[order(mutsig_genes1$gender), ]
log_mutsig_mf <- log(1+mutsig_mf[-(1:3)])

log_mutsig_mf.m <- data.matrix(log_mutsig_mf)
mutsig_heatmap = heatmap.2(log_mutsig_mf.m, trace="none", col = my_palette, margins=c(9,12), labRow = NA, scale = "column")

#Genes from PLOS1 2013 paper implicated in HNSC cancers
hnsc.d <- data.frame(c("RPA2", "E2F2", "FGFR3", "SOX2", "NFE2L2", "KEAP1", "PIK3CA", "AKR1C1", "PDGFRA", "PDGFRB", "TWIST1", "TP63", "TGFA", "EGFR"))
names(hnsc.d) <- "Protein_Coding"

hnsc.g <- merge(hnsc.d, GeneId, by="Protein_Coding", all.x = TRUE)
hnsc_gene <- as.data.frame(t(hnsc.g))
colnames(hnsc_gene) <- as.vector(hnsc.g$Entrez_Gene_Id)

hnsc_gene1 <- cbind(RNASeq2Gene_g[,1:3], RNASeq2Gene_g[ , match(colnames(hnsc_gene), colnames(RNASeq2Gene_g))])
colnames(hnsc_gene1)[-(1:3)] <- as.vector(hnsc.g$Protein_Coding)

plot_mf(hnsc_gene1, "Genes assoc with HNSC, PlosOne 2013")

#### Heat map ###

#The data needs to be in a matrix format with only numbers. The row headers also need to be saved in another vector

RNASeq2Gene_mf <- RNASeq2Gene_g[order(RNASeq2Gene_g$gender),] # To separate males and females

log_mf <- log(1+(RNASeq2Gene_mf[,-(1:3)])) #log(1+e)

gender.v <- RNASeq2Gene_mf$gender #Store the gender values
table(gender.v) #151 female and 415 male

RNASeq2Gene_log_mf <- data.matrix(log_mf) #Convert to matrix

distr <- as.vector(RNASeq2Gene_log_mf[1:151,1:100])
hist(distr)
mean(distr <= 1)
mean(distr == 0)
mean(distr)
median(distr)
quantile(distr, c(0.25, 0.75))
max(distr)

#Plotting the heatmap using heatmap.2 from gplots
install.packages("gplots")
library(gplots)

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

col_breaks = c(seq(0,1,length=100), # for red #Most of the points are 0 or between 0-1
               seq(1, 6,length=100), # for yellow
               seq(6, max(distr),length=100)) # for green

mf_heatmap <- heatmap.2(RNASeq2Gene_log_mf[152:400,100:200], trace="none", Rowv = NA, col = my_palette, margins=c(9,12), labRow = NA, labCol = NA, dendrogram = "column")

