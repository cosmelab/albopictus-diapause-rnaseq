# This script runs DESeq (specifically starting with htseq count files)

##Install DESeq2
# if (!requireNamespace("BiocManager", quietly = TRUE)) 
install.packages("BiocManager")
 BiocManager::install("DESeq2")
 BiocManager::install('EnhancedVolcano')
 BiocManager::install("apeglm")
 BiocManager::install("GOplot")

#Load Libraries
library(DESeq2); library(ggplot2); library(tidyverse); library(EnhancedVolcano); library(GOplot); library(pheatmap); library(reshape2); library(viridis); library(mygene)

# Set working directory to source file location (going to session, set working directory as source file, same functionally as coding it. source file is the script youre working on)

#Choose directory with htseq-count data (this is directory relative to the script)
directory<-("../HTSeqfiles")

#Create the sample table (this could alternatively be made externally and read in) change to 72h and 35 h
sampleFiles <- grep("72h", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically) change the ND to NDI
dds$condition <- relevel(dds$condition, ref = "NDI")

# A bit of quality control
# Look at the distribution of the count data across the samples
librarySizes <- colSums(counts(dds))

barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, ylim = c(0,1.6e+07),
        main="Barplot of library sizes")

logcounts <- log2(counts(dds) + 1)
head(logcounts)

#Is there any difference between per gene counts for each of the sample groups?
statusCol <- as.numeric(factor(dds$condition)) + 1  # make a colour vector

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)

#Adding median log counts
abline(h=median(as.matrix(logcounts)), col="blue")

#Looking at PCA of the data - Do treatments cluster together?
rlogcounts <- rlog(counts(dds))#transforming data to make it approximately homoskedastic, n < 30 so rlog is better

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
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition")

# Hierarchial Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld) #retrieve matrix from the rld object
# compute pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# plot the correlation c=values as a heatmap
pheatmap(rld_cor)

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
res_LFC <- lfcShrink(dds, coef="condition_DI_vs_NDI", type="apeglm")
res_LFC

# Order the table by smallest p value
resOrdered <- res[order(res$pvalue),]
summary(res)


#Basic MA-plot, shinrks log fold change estimates of some replicates compared to other replicates
plotMA(res, ylim=c(-3,3))
plotMA(res_LFC, ylim=c(-3,3)) # See lots of shrinkage towards the beginning
resOrdered <- res[order(res$pvalue),]
#If you want to interactivly find genes associated with certain points
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`. Did not work first time but did second time through
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

# Can look at which of the result are significant and have a high enough log2 fold change . This command did not produce a result for me (?)
sig_res <- res_LFC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes. was told that attempt to set col.names was ignored?
write.csv(select(sig_res, gene, log2FoldChange, padj), 
            file="../misc/72h_DvND_LFCshrink_padj.txt", col.names = F, row.names = F)

# Write out just the gene names for later analysis in KEGG, need to rename this file to day number of this set)
write.table(sig_res %>% select(gene), 
            file="../misc/72h_DvND_test.txt", col.names = F, row.names = F, quote = F)


####################
# Running through another dataset (135hdays)
#Create the sample table (this could alternatively be made externally and read in)
sampleFiles <- grep("135h", list.files(directory), value=T) #this is pulling all files from the directory with the age specified
sampleNames <- sub("_htseqCount","",sampleFiles) #this is removing the ending of the files to better represent the sample names
sampleConditions <- sub("_.*","", sampleFiles) #to get conditions I pull everything from before the first underscore (which I know is the diapause status)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleConditions)
str(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)

# Make the DESeq dataset
dds_135h <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds_135h

#DESeq recommends a pre-filtering step to reduce memory size and increase speed. They suggest keeping only rows which have 10 reads total
keep <- rowSums(counts(dds_135h)) >= 10
dds_135h <- dds_135h[keep,]

#Relevel the condition for what to compare to. Likely not important in our case with two conditions, but you would want everything compared to the control. (Default is first condition alphabetically)
dds_135h$condition <- relevel(dds_135h$condition, ref = "NDI")


#Run DESeq
dds_135h <- DESeq(dds_135h)
res_135h <- results(dds_135h)
res_135h

# Apparently, it is common to shrink the log fold change estimate for better visualization and ranking of genes. Often, lowly expressed genes tend to have relatively high levels of variability so shrinking can reduce this. 
resultsNames(dds_135h)
res_LFC_135h <- lfcShrink(dds_135h, coef="condition_DI_vs_NDI", type="apeglm")
res_LFC_135h


# Volcano plot of the lfcShrink data. For various setting try: `browseVignettes("EnhancedVolcano")`
EnhancedVolcano(res_LFC_135h,
                lab = rownames(res_LFC_135h),
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
sig_res_135h <- res_LFC_135h %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as.data.frame() %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 0.58)

# Write out a table of these significant differentially expressed genes
write.csv(select(sig_res_135h, gene, log2FoldChange, padj), 
          file="../misc/135h_DIvNDI_LFCshrink_padj.txt", row.names = F)

# Write out just the gene names for later analysis in KEGG
write.table(sig_res_135h %>% select(gene), 
            file="../misc/135h_DvND_test.txt", col.names = F, row.names = F, quote = F)

##################
# Making a venn diagram of the 2 different datasets DEGs
dat_72 <- read.csv("../misc/72h_DvND_LFCshrink_padj.txt")
dat_135 <- read.csv("../misc/135h_DIvNDI_LFCshrink_padj.txt")



l1 <- dat_72[, c(1,2)]
l2 <- dat_135[, c(1,2)]


VennDiag <- GOVenn(l1,l2, label=c('72h','135h'), plot = F)
print(VennDiag$plot)
# The colors inside the labels: Red (up), Blue (down), and Yellow (contra-regulated). These colors can be changed but the order of the regulation (up, down, contra) starts at the top and goes clockwise.

NEW Additions: 
res_print <- as.data.frame(res_LFC) # first for 72h
res_135h_print <- as.data.frame(res_LFC_135h)

NEW: 
# Make output excel sheet to match M.F.P.'s old excel sheets
res_print <- as.data.frame(res_LFC)
res_135h_print <- as.data.frame(res_LFC_135h)


res_merged <- merge(res_print,res_135h_print, by="row.names") 
#renamed with 72h and 135h- just want to make sure i did it in the correct order 
colnames(res_merged) <- c("geneID", "72h_baseMean", "72h_Log2FoldChange", "72h_lfcSE", "72h_pvalue", "72h_padj", "135h_baseMean", "135h_Log2FoldChange", "135h_lfcSE", "135h_pvalue", "135h_padj")

head(res_merged)

# Find gene names with mygene package
dat <- queryMany(res_merged$geneID, scopes="symbol", fields="name")

res_merged <- merge(res_merged, dat, by.x="geneID", by.y="query")

keep_cols <- c("geneID", "name", "72h_Log2FoldChange", "72h_padj", "135h_Log2FoldChange", "135h_padj")

# Write out a csv with these data
write.csv(res_merged[keep_cols], 
          file="../misc/DESeq_results_embryo.csv", row.names = F)
