#To find the difference between tumors and normal samples from the same patients. I will go back to the original file RNASeq2GeneNorm. In the file Mean0_signif_genes.R, I removed duplicates that were normal tissues. Let's first see how many females and males have the data required.

# IMPORTANT :: See the revised version below using adjusted p-values.

library(dplyr)
library(ggplot2)
library(lazyeval)
library(tidyr)
library(reshape2)

data <- RNASeq2GeneNorm
gender <- subset(ClinicalData, select = c(Patient_Id, gender))
data.gend <- merge(gender, data, by = "Patient_Id")

#Get rows for patients with both tumor and normal samples. I first filtered to find normal samples (sample_Id = 11) and then used the patient id to get the rows required. The duplicated command could be used but it also gives metastatic tissue.

dat.tn.pat <- dplyr::filter(data.gend, Sample_Id == "11")
dat.tn <-  data.gend[(grep(paste(dat.tn.pat$Patient_Id, collapse = "|"), data.gend$Patient_Id)), ]
table(dat.tn$gender)/2 #Half of these are the numbers of F & M, there is one extra normal tissue from M
#female   male 
#28     59 
#14     29

colnames(dat.tn)[-(1:3)] <- paste("X", colnames(dat.tn)[-(1:3)], sep = "")

#The data set is not large. To find out the difference between tumor & normal

dat.tn <- dplyr::arrange(dat.tn, Patient_Id)

test <- dat.tn[ , (1:10)]
test <- test[-79,] #Remove extra normal tissue

val <- colnames(test)[7]
test.tn <- test %>% group_by(Patient_Id) %>% summarise_(interp(~diff(var), var = as.name(val)))

test.tn <- merge(gender, test.tn, by = "Patient_Id")
t.X100133144 <- t.test(col4~gender, data = test.tn)
femt.X100133144 <- filter(test.tn, gender == "female") %>% select(col4) %>% unlist
malt.X100133144 <- filter(test.tn, gender == "male") %>% select(col4) %>% unlist
qqnorm(test.tn$col4, col = gender)
qqline(femt.X100133144)
qqnorm(malt.X100133144)
qqline(malt.X100133144)
abline(0,1)

#Putting a loop for genes
test.tn <- data.frame(Patient_Id = unique( sort(test$Patient_Id) ) )
test.tn <- merge(test.tn, gender, by = "Patient_Id")
for(i in 4:10){
  val <- colnames(test)[i]
  test.tn[i-1] <- test %>% group_by(Patient_Id) %>% summarise_(interp(~diff(var), var = as.name(val))) %>% select(-(Patient_Id))
}

#Applying to all the genes: In order to get relevant results, I will consider only the genes with maximum variance
sd.tn <- sapply(dat.tn[ ,-(1:3)], FUN = sd)
sd.tn <- sort(sd.tn, decreasing = TRUE)
plot(sd.tn, main = "SD distribution", xlab = "Genes", ylab = "Std Dev")
abline(v = 2000, col = 2)
text(2000, 0.5, "top 2000", col = 2, pos = 4)

#Consider top 10% of genes with highest variance
sd2000 <- data.frame(t(sd.g[1:2000]))
tn2000 <- dat.tn[ ,match(colnames(sd2000), colnames(dat.tn))]
g2k <- cbind(dat.tn[ ,1:3], tn2000)

g2k <- g2k[-79, ]
diff.tn <- data.frame(Patient_Id = unique(g2k$Patient_Id))
diff.tn <- merge(diff.tn, gender, by = "Patient_Id")

for(i in 4:ncol(g2k)){
  val <- colnames(g2k)[i]
  diff.tn[i-1] <- g2k %>% group_by(Patient_Id) %>% summarise_(interp(~diff(var), var = as.name(val))) %>% select(-(Patient_Id))
}

colnames(diff.tn) = colnames(g2k)[-3]
#colnames(diff.tn) = sapply(strsplit(colnames(diff.tn), split = "\\(|\\)"), "[", 2)
Pval = Tstat = vector()
rows = nrow(diff.tn)

for(i in 3:ncol(diff.tn)){
  tstat <- t.test(diff.tn[ ,i] ~ diff.tn[ ,2])
  Pval[i-2] = tstat$p.value
  Tstat[i-2] = tstat$statistic
}

signif_genes = data.frame(colnames(diff.tn)[which(Pval<0.01)+2])
colnames(signif_genes)[1] = "Entrez_Gene_Id"
signif_genes.tn <- merge(signif_genes, GeneId, by = "Entrez_Gene_Id", all.x = TRUE) ### Signif gene names
colnames(signif_genes.tn) = c("Entrez.Gene.Id" , "Approved.Symbol")
genenames.tn <- merge(signif_genes.tn, gene.names, by = "Approved.Symbol")

numeric_mf = as.numeric(as.factor(diff.tn$gender))
qqnorm(diff.tn$X10351, col = numeric_mf+1) #Male is red, female is black

signif_gene_expr = t(signif_genes)
colnames(signif_gene_expr) = signif_gene_expr[1,]
signif_gene_expr.tn = diff.tn[ , match(colnames(signif_gene_expr), colnames(diff.tn))]
signif_gene_expr.tn = cbind(diff.tn[ ,1:2] ,signif_gene_expr.tn)

genes = signif_gene_expr.tn[ ,-c(3,4,24)]

library(ggplot2)
dat.m = tidyr::gather(genes, "genes", "exprn", 3:23)
#dat.m <- melt(signif_gene_expr.tn, id.vars = gender, measure.vars = colnames(signif_gene_expr.tn)[-(1:2)])
ggplot(dat.m) + geom_boxplot(aes(x=genes, y=exprn, col = gender)) + ggtitle("(Tumor - normal) in M vs F") + ylab("gene expression")

### Revision: Using adjusted P-val #######

datTN.log <- cbind(dat.tn[ ,1:3], log(1+ dat.tn[ ,-(1:3)])) #log(1+e)
sd.tn <- sapply(datTN.log[ ,-(1:3)], FUN = sd)
qqnorm(sd.tn)
qqline(sd.tn) 
quants <- qqnorm(sd.tn)
which(quants$x == 0)
#10106 is the point at x=0
abline(0,1)
abline(v=0, col = 2)

#Take 10,000 genes with highest expression variance
sd.tn <- sort(sd.tn, decreasing = TRUE)
plot(sd.tn)
sd10k <- data.frame(t(sd.tn[1:10000]))
tn2000 <- datTN.log[ ,match(colnames(sd10k), colnames(datTN.log))]
g10k <- cbind(datTN.log[ ,1:3], tn2000)
g10k <- dplyr::arrange(g10k, Patient_Id)
g10k <- g10k[-79, ]

#Top 10,000 genes with highest variance for the sample dataset

#Using diff on a column will give the difference between consecutive rows. However, I want the difference only between row1 & row2, row3 & row4, etc. Using sapply on diff and then selecting for the odd number of rows. Appending these to the Patitent_Id and gender info for the samples. Much easier and faster!!

test.tn = sapply(g10k[ ,-(1:3)], diff)
odd_num = seq(1, nrow(g10k), 2)

diff.tn <- data.frame(Patient_Id = unique(g10k$Patient_Id))
diff.tn <- merge(diff.tn, gender, by = "Patient_Id")
diff.tn = cbind(diff.tn, test.tn[odd_num, ])

diff.tn = arrange(diff.tn, gender)
diff.tn[ , 3:nrow(diff.tn)] %>% sapply(mean) %>% summary

#For computingt t-test for all genes over gender (column 2), I made a funciton that can be fed to sapply

Tstat = function(col){
  tstat = t.test(col ~ diff.tn[ ,2])
}
Tstat.tn = sapply(diff.tn[ , -(1:2)], Tstat)
row.names(Tstat.tn)
#[1] "statistic"   "parameter"   "p.value"     "conf.int"    "estimate"    "null.value"  "alternative" "method"     
#[9] "data.name" 

Pval = Tstat.tn[3, ]

#The p-value can be adjusted using the FDR method.
Padj = p.adjust(Pval, method = "fdr")

Pval.tn = unlist(Pval)
Pval.tn1 = -log(Pval.tn)
Padj.tn = unlist(Padj)
Padj.tn1 = -log(Padj.tn)

Padj.tn = unlist(Padj)
summary(Padj.tn)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000976 0.5687000 0.7753000 0.7331000 0.9126000 0.9999000 

#The genes called significant through adjusted Pval = Padj
sum(Padj.tn < 0.1) #FDR < 1
#7
sum(Padj.tn < 0.2)
#36
which(Padj.tn < 0.1)
#X11240 X55410 X90665 X83869  X4171   X993 X54414 
#   521   1176   1886   2256   7435   9035   9658 
summary(Pval.tn)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1422  0.3877  0.4238  0.6845  0.9999 

#The q value is an alternate method for the FDR (note that there are some differences from p.adjust)
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
browseVignettes("qvalue")

qobj = qvalue::qvalue(p = Pval.tn)
summary(qobj)
#pi0:	0.7878925	

#Cumulative number of significant calls:
#          <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value        7     46   291    631  1140 1924 10000
#q-value        2      2     2      2     3   10 10000
#local FDR      0      2     2      2     2    6  7043

#estimate of the proportion of true alternative tests = 1-pi0
1-0.7878925
#0.2121075 

#calculates the maximum estimated q-value among all p-values ≤ 0.01, which is equivalent to estimating the false discovery rate when calling all p-values ≤ 0.01 significant.
max(qval[qobj$pvalues <= 0.01]) #0.2676485

hist(qobj)
plot(qobj)
qval = qobj$qvalues
which(qval < 0.01)
# X55410 X83869 X54414 
#   1176   2256   9658 

which(qval < 0.1) # 1 gene out of 10 is expected to be false positive
#X11240   X9087  X55410  X90665  X83869 X122773   X4171  X79895    X993  X54414 
#   521     630    1176    1886    2256    5105    7435    7853    9035    9658 

qval[which(qval < 0.1)] #Indiv q-values denote proportion of genes that are as or more extreme are FPs 
#X11240        X9087       X55410       X90665       X83869      X122773        X4171       X79895         X993 
#7.713199e-02 9.819673e-02 7.685609e-05 6.349179e-02 7.685609e-05 9.819673e-02 6.349179e-02 9.819673e-02 6.349179e-02 
#X54414 
#4.831142e-02 

#Significance of t-test
Tstat.tn[1, ] %>% unlist %>% summary
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-4.6130 -0.6686  0.1861  0.2040  1.0800  7.8010 

tstat = Tstat.tn[1, ] %>% unlist
hist(tstat, prob=TRUE)
lines(density(tstat), col = "red", lwd = 2)
which(tstat > 5) #13 genes
#X9086.t  X11240.t   X9087.t   X4137.t  X55410.t   X9472.t  X90665.t  X10580.t  X83869.t X122773.t  X79895.t  X57010.t 
#   122       521       630       836      1176      1762      1886      1988      2256      5105      7853      8435 
#X54414.t 
#   9658

# t-statistic for significant genes identified by Padj. The negative value only indicates it is significant in opposite direction.
sig = which(Padj.tn < 0.1) %>% unlist
sigTstat = Tstat.tn[ , sig]
sigTstat[1, ] %>% unlist
# X11240.t  X55410.t  X90665.t  X83869.t   X4171.t    X993.t  X54414.t 
# 4.469653  7.744435  4.796699  7.800996 -4.592771 -4.613302  4.954153 

Tstat.tn[1, ] %>% unlist %>% summary
#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   -4.6130 -0.6686  0.1861  0.2040  1.0800  7.8010 

#Plots with Pval and Padj
par(mfrow = c(2, 2))
qqnorm(Pval.tn1, main = "QQPlot of -log(P) for T-N in F vs M")
qqline(Pval.tn1)
hist(Pval.tn, main =  "Histogram of P-val")
qqnorm(Padj.tn1, main = "QQPlot of -log(Padj) for T-N in F vs M")
qqline(Padj.tn1)
hist(Padj.tn, main =  "Histogram of P-adj")

plot(Pval.tn, Padj.tn) #Show how p-values have been adjustment

#Plot showing distribution of Pval and Padj as function of t-statistic
df = data.frame(x = Tstat.tn[1, ], y1 = Pval.tn, y2 = Padj.tn)
plot(Tstat.tn[1, ], Padj.tn, type = 'p', main = "P-val and Adjusted P-val vs. t-statistic")
par(new = TRUE)
plot(Tstat.tn[1, ], 
     Pval.tn, col = "blue", type = 'p', 
     axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, 
     at = pretty(range(Pval.tn)), col = "blue")

#Plot Significant genes as a function of Padj cutoff. 
interval1 = seq(0,0.5, by = 0.001) #Define axis intervals for plotting
interval2 = seq(0,0.1, by = 0.001)

signifval = function(num){ #Function to calculate cum sum for signif genes
  sum(Padj.tn < num)
}
CumSignif1 = sapply(interval1, signifval)
CumSignif2 = sapply(interval2, signifval)

par(mfrow = c(2,1))
plot(interval1, CumSignif1, ylab = "Number of signif genes", xlab = "P-adjusted", type = 'l', main = "Number of signif genes vs adjusted p value cutoffs")
plot(interval2, CumSignif2, ylab = "Number of signif genes", xlab = "P-adjusted", type = 'l')

#The same genes pop up using different measures for significance. I will use the 10 genes identified with q-value < 0.1. 
library(lazyeval)
library(plotly)
library(ggplot2)
source("/Users/admin/multiplot.r")

signif_tstat = Tstat.tn[ , which(qval < 0.1)]
signif_genes = diff.tn[ , (which(qval < 0.1) + 2)]

diff.tn$gender = as.factor(diff.tn$gender)
diff.tn$gender = relevel(diff.tn$gender, ref = "male")
numeric_mf = as.numeric(diff.tn$gender)
par(mfrow = c(3, 4))
qqnorm(signif_genes[ ,1], col = numeric_mf, main = "X11240") #Male is black, female is red
qqnorm(signif_genes[ ,2], col = numeric_mf, main = "X9087") 
qqnorm(signif_genes[ ,3], col = numeric_mf, main = "X55410") 
qqnorm(signif_genes[ ,4], col = numeric_mf, main = "X90665") 
qqnorm(signif_genes[ ,5], col = numeric_mf, main = "X83869") 
qqnorm(signif_genes[ ,6], col = numeric_mf, main = "X122773") 
qqnorm(signif_genes[ ,7], col = numeric_mf, main = "X4171") 
qqnorm(signif_genes[ ,8], col = numeric_mf, main = "X79895") 
qqnorm(signif_genes[ ,9], col = numeric_mf, main = "X993") 
qqnorm(signif_genes[ ,10], col = numeric_mf, main = "X54414") 

plotq = function(i){
  column = colnames(signif_genes)[i]
  qplot(sample = signif_genes[ , i], data = signif_genes, color = gender) + labs(title= paste0("Gene", column))
}

p1 <- qplot(sample = signif_genes[ , 1], color = gender) + labs(title= "X11240") + theme(legend.position="none")
p2 <- qplot(sample = signif_genes[ , 2], color = gender) + labs(title= "X9087") + theme(legend.position="none")
p3 <- qplot(sample = signif_genes[ , 3], color = gender) + labs(title= "X55410") + theme(legend.position="none")
p4 <- qplot(sample = signif_genes[ , 4], color = gender) + labs(title= "X90665")+ theme(legend.position="none")
p5 <- qplot(sample = signif_genes[ , 5], color = gender) + labs(title= "X83869")+ theme(legend.position="none")
p6 <- qplot(sample = signif_genes[ , 6], color = gender) + labs(title= "X122773") + theme(legend.position="none")
p7 <- qplot(sample = signif_genes[ , 7], color = gender) + labs(title= "X4171") + theme(legend.position="none")
p8 <- qplot(sample = signif_genes[ , 8], color = gender) + labs(title= "X79895") + theme(legend.position="none")
p9 <- qplot(sample = signif_genes[ , 9], color = gender) + labs(title= "X993") + theme(legend.position="none")
p10 <- qplot(sample = signif_genes[ , 10], color = gender) + labs(title= "X54414") + theme(legend.position="none")

multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, cols = 3)

#Extracting gene information for the genes called significant
signif_gene_names = data.frame(colnames(signif_genes))
colnames(signif_gene_names)[1] = "Entrez_Gene_Id"
signif_genes.tn <- merge(signif_gene_names, GeneId, by = "Entrez_Gene_Id", all.x = TRUE) ### Signif gene names
colnames(signif_genes.tn) = c("Entrez.Gene.Id" , "Approved.Symbol")
genenames.tn <- merge(signif_genes.tn, gene.names, by = "Approved.Symbol")

#Only 9 genes extracted since X55410 or NCRNA00185 is the same as X83869(TTTY14).
write.table(genenames.tn, file = "Tumor-Normal_diff_geneexp_MvsF.txt", sep = "\t")

signif_genes = cbind(diff.tn[ ,1:2], signif_genes)
dat.m = tidyr::gather(signif_genes, "genes", "exprn", 3:12)
#dat.m <- melt(signif_gene_expr.tn, id.vars = gender, measure.vars = colnames(signif_gene_expr.tn)[-(1:2)])
ggplot(dat.m) + geom_boxplot(aes(x=genes, y=exprn, col = gender)) + ggtitle("(Tumor - normal) in M vs F") + ylab("gene expression") + coord_flip()

signif_genes = diff.tn[ , (which(qval < 0.2) + 2)]

### Wilcoxon test #####
wilcox.test(diff.tn[  , 523] ~  diff.tn[ , 2])

Wstat = function(col){
  wstat = wilcox.test(col ~ diff.tn[ ,2])
}
Wstat.tn = sapply(diff.tn[ , -(1:2)], Wstat)
Pval.W = Wstat.tn[3, ] %>% unlist
summary(Pval.W)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000023 0.175700 0.433900 0.452400 0.710400 1.000000 
sum(Pval < 0.001)
#46

#q value doesn't seem to be a good adjustment for ranked test.
qobj.W = qvalue::qvalue(p = Pval.W)
summary(qobj.W)
#pi0:	0.9081073	

#Cumulative number of significant calls:
  
#          <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value        2     23   201    453   874 1548  9893
#q-value        0      0     0      0     0    0 10000
#local FDR      0      0     0      0     0    0  4708

hist(qobj.W)
plot(qobj.W)

which(Pval.W < 0.001)
#X3239  X11240   X8838   X4137  X55410 X387890  X90665  X10580   X1524  X83869  X22998 X220965 X131890 X122773 X168537 
# 382     521     835     836    1176    1520    1886    1988    2223    2256    2400    2577    4470    5105    5484 
#X114907   X1869 X474344   X4171  X54741  X57010    X993  X54414 
#   6780    7215    7223    7435    8274    8435    9035    9658 
signif_genesW = diff.tn[ , (which(Pval.W < 0.001) + 2)]
signif_genesW = cbind(diff.tn[ ,1:2], signif_genesW)
dat.w = tidyr::gather(signif_genesW, "genes", "exprn", 3:25)
#dat.m <- melt(signif_gene_expr.tn, id.vars = gender, measure.vars = colnames(signif_gene_expr.tn)[-(1:2)])
ggplot(dat.w) + geom_boxplot(aes(x=genes, y=exprn, col = gender)) + ggtitle("(Tumor - normal) in M vs F Wilcoxon") + ylab("gene expression") + coord_flip()
