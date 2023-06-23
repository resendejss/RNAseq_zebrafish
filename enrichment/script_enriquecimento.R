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
