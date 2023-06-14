################################################################################
## RNA-Seq - zebrafish
## GO enrichment
## date of creation: 2023/06/12
## date of last update: 2023/06/14
## jean resende
################################################################################
library(clusterProfiler)
library(org.Dr.eg.db)

# -- AHcut0_CTcut0 all samples--
data.genes <- read.csv("../differential_expression/allSamples/AHcut0_CTcut0_res05_sig_fc0.csv")

data <- data.genes[data.genes$lfcSE > 1 | data.genes$log2FoldChange < -1,]
data <- data[!is.na(data$padj),]

nrow(data.genes[data.genes$log2FoldChange < 1,])


gene <- data$X

entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
entrez_genes <- entrez_genes[!is.na(entrez_genes)]

ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH"
                )
goplot(ego)
head(ego)

#ggo <- groupGO(gene     = gene.ahctuncut,
#               OrgDb    = org.Dr.eg.db,
#               keyType = "ENSEMBL",
#               ont      = "CC",
#               #level    = 3,
#               readable = TRUE)
#head(ggo)

# -- KEEG
entrez_genes <- mapIds(org.Dr.eg.db, gene.ahctuncut, 'ENTREZID', 'ENSEMBL')

kegg <- enrichKEGG(gene = entrez_genes,
                   organism = "dre",
                   pvalueCutoff = 0.05)

head(kegg)


browseKEGG(kegg, "dre00100")