# This script runs DESeq (specifically starting with htseq count files)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install('EnhancedVolcano')
# BiocManager::install("apeglm")
# BiocManager::install("GOplot")

#Load Libraries
#We either just installed these or installed these in the last script (countreads)
library(DESeq2); library(ggplot2); library(tidyverse); library(EnhancedVolcano); library(GOplot); library(pheatmap)

# Set working directory to source file location (scripts folder)

#Choose directory with htseq-count data (relative to scripts folder)
directory<-("../Data")

#Create the sample table (this could alternatively be made externally and read in)
#grabbing all the files in the directory with NB and setting that as sampleFiles
sampleFiles <- grep("NB", list.files(directory), value=T) #this is pulling all files from the directory with the NB (no blood meal) because we want to sperate comparisons by blood meal and no blood meal
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions (long days/short day i.e. non diapause/diapause) I pull everything from before the first underscore (which I know is the diapause status)

#a factor=levels for a certain variable; telling deseq it is a variable. for example, 1=LD and 2=ND and then ask R to compare 1 and 2
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable) #str allows you to view the internal structure of an R object
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds

#getting rid of rows with less than 10 reads 
#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically, in this case, it happens to be the control that we want: LD or nondiapause)
dds$condition <- relevel(dds$condition, ref = "LD")

# A bit of quality control
# Look at the distribution of the count data across the samples, will make a boxplot, check to see if the sizes are comparable accross samples. Stats things to make sure our data looks normal. can export and save as image
# PCA longest axis, explains most variation of the data and seperates it based on that
librarySizes <- colSums(counts(dds))

#saved image as "barplot of library sizes"; it gives you an iea of count distribution across samples
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

#difference between per gene counts for each of the sample groups
#saved image as logcounts
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

#Looking at PCA of the data - Do treatments cluster together?
rlogcounts <- rlog(counts(dds)) #transforming data to make it approximately homoskedastic, n < 30 so rlog is better

data_for_PCA <- t(rlogcounts)
dim(data_for_PCA)

#run PCA
pcDat <- prcomp(data_for_PCA, center = T)
pcDat_df <- data.frame('Condition' = dds$condition, pcDat$x[,1:2])

# basic plot
ggplot(pcDat_df, aes(x = PC1, y = PC2, col = Condition)) +
  geom_point() + 
  theme_bw()

# Transform normalized counts using the rlog function ("RNASEQ20_Day3_HandsOn.pdf")
# see what percentage of the variance can be explained by PC1 (LD vs SD i.e. nondiapause vs diapause)
#saved image as PCA counts
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")

# Hierarchial Clustering
### Extract the rlog matrix from the object
# Heat map: look out for a row or block that are either all red or all blue
rld_mat <- assay(rld) #retrieve matrix from the rld object
# compute pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# plot the correlation c=values as a heatmap
heatmap(rld_cor)

# # Look at the difference in total reads from the different samples
# dds <- estimateSizeFactors(dds)
# 
# # Estimate the dispersion (or variation) or the data
# dds <- estimateDispersions(dds)
# 
# plotDispEsts(dds)

# Differential Expression Analysis
#Run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds)
res_LFC <- lfcShrink(dds, coef="condition_SD_vs_LD", type="apeglm")
res_LFC

# Order the table by smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)


#Basic MA-plot
plotMA(res, ylim=c(-3,3))
plotMA(res_LFC, ylim=c(-3,3)) # See lots of shrinkage towards the beginning

#If you want to interactivly find genes associated with certain points
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC,
                lab = rownames(res_LFC),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(select(sig_res, gene, log2FoldChange, padj), 
          file="../misc/NB_LDvSD_LFCshrink_padj.txt", col.names = F, row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res %>% select(gene), 
            file="../misc/NB_LDvSD_test.txt", col.names = F, row.names = F, quote = F)


####################
# Running through another dataset (BM)
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- grep("BM", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds_BM <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                     directory = directory,
                                     design = ~ condition)
dds_BM

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds_21d)) >= 10
dds_BM <- dds_BM [keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds_BM$condition <- relevel(dds_BM$condition, ref = "LD")


#Run DESeq
dds_BM <- DESeq(dds_BM)
res_BM <- results(dds_BM)
res_BM

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this.
resultsNames(dds_BM)
res_LFC_BM <- lfcShrink(dds_BM, coef="condition_SD_vs_LD", type="apeglm")
res_LFC_BM


# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC_BM,
                lab = rownames(res_LFC_BM),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NA,
                #drawConnectors = TRUE,
                xlim = c(-1.5, 1.5),
                ylim = c(0,30),
                pCutoff = 10e-6,
                FCcutoff = 0.58,
                pointSize = 2.0,
                labSize = 5.0)

# Can look at which of the result are significant and have a high enough log2 fold change
sig_res_BM <- res_LFC_BM %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(dplyr::select(sig_res_BM, gene, log2FoldChange, padj), 
          file="../misc/BM_SDvLD_LFCshrink_padj.txt", row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res_BM %>% dplyr::select(gene), 
            file="../misc/BM_SDvLD_test.txt", col.names = F, row.names = F, quote = F)
##################
# Making a venn diagram of the 3 different datasets DEGs
dat_NB <- read.csv("../misc/NB_SDvLD_LFCshrink_padj.txt")
dat_BM <- read.csv("../misc/BM_SDvLD_LFCshrink_padj.txt")



l1 <- dat_NB[, c(1,2)]
l2 <- dat_BM[, c(1,2)]

VennDiag <- GOVenn(l1,l2, label=c('NB','BM'), plot = F)
print(VennDiag$plot)
# The colors inside the labels: Red (up), Blue (down), and Yellow (contra-regulated). These colors can be changed but the order of the regulation (up, down, contra) starts at the top and goes clockwise.


# Make output excel sheet to match M.F.P.'s old excel sheets
#remember that res_LFC is the name used for res_LFC_NB
res_NB_print <- as.data.frame(res_LFC)
res_BM_print <- as.data.frame(res_LFC_BM)

#set res_NB_print as a vector with row names
res_NB_print$Row.names <- row.names(res_NB_print)

#merge row names from res_print of NB and BM
res_merged <- merge(res_NB_print,res_BM_print, by="row.names")

#colnames: retrieve or set a row of column names of a matrix-like object. In this case, setting the column names to following variables inside the vector 
colnames(res_merged) <- c("geneID", "NB_baseMean", "NB_Log2FoldChange", "NB_lfcSE", "NB_pvalue", "NB_padj", "BM_baseMean", "BM_Log2FoldChange", "BM_lfcSE", "BM_pvalue", "BM_padj")

#made sure that the merged table looks ok, i.e. making sure that the table has both NB and BM data
head(res_merged)

BiocManager::install("mygene")
library(mygene)
# Find gene names with mygene package
dat <- queryMany(res_merged$geneID, scopes="symbol", fields="name")
res_merged <- merge(res_merged, dat, by.x="geneID", by.y="query")
keep_cols <- c("geneID", "name", "NB_Log2FoldChange", "NB_padj", "BM_Log2FoldChange", "BM_padj")

# Write out a csv with these data# Write out a csv with these data
write.csv(res_merged[keep_cols], 
          file="../misc/DESeq_results_adultlarvae.csv", row.names = F)
