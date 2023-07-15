library(biomaRt)
library(clusterProfiler)
library(org.Dr.eg.db)

ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
head(datasets)

setwd("~/Documentos/RNAseq_zebrafish/enrichment")
gene_id <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                 mart = useDataset("drerio_gene_ensembl", useMart("ensembl")))

# -- AHcut0_CTcut0 all samples-- ###############################################
data.genes <- read.csv("../differential_expression/allSamples/AHcut0_CTcut0_res05_sig_fc0.csv")
genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)

# -- AHcut0_CTcut0 filt samples-- ##############################################

data.genes <- read.csv("../differential_expression/filtSamples/AHcut0_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)

# -- AHuncut0_CTuncut0 all samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)

browseKEGG(kk, 'dre00190')
browseKEGG(kk, 'dre04260')
browseKEGG(kk, 'dre00983')
browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00120')
browseKEGG(kk, 'dre03320')
browseKEGG(kk, 'dre00980')
browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre00982')
browseKEGG(kk, 'dre00040')
browseKEGG(kk, 'dre00591')
browseKEGG(kk, 'dre00100')
browseKEGG(kk, 'dre04744')
browseKEGG(kk, 'dre00590')
browseKEGG(kk, 'dre03010')

# -- AHuncut0_CTuncut0 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)

browseKEGG(kk, 'dre00190')
browseKEGG(kk, 'dre04260')
browseKEGG(kk, 'dre00983')
browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00120')
browseKEGG(kk, 'dre03320')
browseKEGG(kk, 'dre00980')
browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre00982')
browseKEGG(kk, 'dre00040')
browseKEGG(kk, 'dre00591')
browseKEGG(kk, 'dre00100')
browseKEGG(kk, 'dre04744')
browseKEGG(kk, 'dre00590')
browseKEGG(kk, 'dre03010')


# -- AHuncut0_CTuncut0 filt samples-- ##########################################
data.genes <- read.csv("../differential_expression/filtSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)

browseKEGG(kk, 'dre00190')



# -- CTuncut0_CTcut0 all samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/CTuncut0_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre04920')
browseKEGG(kk, 'dre04060')
browseKEGG(kk, 'dre04217')
browseKEGG(kk, 'dre04625')
browseKEGG(kk, 'dre04216')
browseKEGG(kk, 'dre05168')
browseKEGG(kk, 'dre04621')

# -- CTuncut0_CTcut0 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/CTuncut0_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre04920')
browseKEGG(kk, 'dre04060')
browseKEGG(kk, 'dre04217')
browseKEGG(kk, 'dre04625')
browseKEGG(kk, 'dre04216')
browseKEGG(kk, 'dre05168')
browseKEGG(kk, 'dre04621')

# -- AHcut5_CTcut5 All samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/AHcut5_CTcut5_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre04814')

# -- AHcut5_CTcut5 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/AHcut5_CTcut5_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00190')
browseKEGG(kk, 'dre04260')
browseKEGG(kk, 'dre03010')

# -- AHcut100_CTcut100 all samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/AHcut100_CTcut100_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00982')
browseKEGG(kk, 'dre00980')
browseKEGG(kk, 'dre00591')
browseKEGG(kk, 'dre00983')
browseKEGG(kk, 'dre03320')
browseKEGG(kk, 'dre00040')
browseKEGG(kk, 'dre00053')
browseKEGG(kk, 'dre00590')
browseKEGG(kk, 'dre01240')
browseKEGG(kk, 'dre00860')
browseKEGG(kk, 'dre00360')
browseKEGG(kk, 'dre00120')
browseKEGG(kk, 'dre00260')
browseKEGG(kk, 'dre04744')
browseKEGG(kk, 'dre00100')
browseKEGG(kk, 'dre00500')
browseKEGG(kk, 'dre00130')
browseKEGG(kk, 'dre00350')
browseKEGG(kk, 'dre00010')

# -- AHcut100_CTcut100 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/AHcut100_CTcut100_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00190')
browseKEGG(kk, 'dre04260')
browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00120')
browseKEGG(kk, 'dre00980')
browseKEGG(kk, 'dre00982')
browseKEGG(kk, 'dre00983')
browseKEGG(kk, 'dre00591')


# -- CTcut100_CTcut0 all samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/CTcut100_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00100')
browseKEGG(kk, 'dre00830')

# -- CTcut100_CTcut0 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/CTcut100_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00830')
browseKEGG(kk, 'dre00100')
browseKEGG(kk, 'dre00140')
browseKEGG(kk, 'dre00040')
browseKEGG(kk, 'dre00053')
browseKEGG(kk, 'dre00860')
browseKEGG(kk, 'dre00982')
browseKEGG(kk, 'dre00980')
browseKEGG(kk, 'dre00983')
browseKEGG(kk, 'dre00500')

# -- CTcut5_CTcut0 all samples-- ###########################################
data.genes <- read.csv("../differential_expression/allSamples/CTcut5_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre04146')

# -- CTcut5_CTcut0 filt samples-- ###########################################
data.genes <- read.csv("../differential_expression/filtSamples/CTcut5_CTcut0_res05_sig_fc0.csv")

genes <- data.genes[data.genes$log2FoldChange > 0.5 | data.genes$log2FoldChange < -0.5,]
genes <- genes[!is.na(genes$padj),]
genes <- genes[genes$padj < 0.05,]
genes <- genes$X
idx <- match(genes,gene_id$ensembl_gene_id)
entrez <- gene_id$entrezgene_id[idx]
entrez <- entrez[!is.na(entrez)]
entrez <- as.character(entrez)

ggo <- groupGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               ont = "CC",
               #level = 3,
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="BP",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="MF",
               readable = T)

kk <- enrichGO(gene = entrez,
               OrgDb = org.Dr.eg.db,
               pAdjustMethod = "BH",
               ont="CC",
               readable = T)

head(kk)
barplot(kk)
goplot(kk)

# -- KEEG
kk <- enrichKEGG(gene = entrez,
                 organism = "dre",
                 pvalueCutoff = 0.05)

head(kk)
kk$ID

browseKEGG(kk, 'dre00330')
