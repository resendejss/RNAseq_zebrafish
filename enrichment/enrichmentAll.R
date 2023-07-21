################################################################################
## RNA-Seq - zebrafish
## GO enrichment
## date of creation: 2023/06/12
## date of last update: 2023/06/16
## jean resende
################################################################################
library(clusterProfiler)
library(org.Dr.eg.db)

# -- Down -- ###################################################################
## -- GO
data.genes <- read.csv("../differential_expression/allSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")

data <- data.genes[data.genes$log2FoldChange < -0.5,]
data <- data[!is.na(data$padj),]
data <- data[data$padj < 0.05,]
nrow(data)

gene <- data$X

entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
entrez_genes <- entrez_genes[!is.na(entrez_genes)]

### -- BP
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH")
table(ego@result$ONTOLOGY)
data_egoo <- ego@result


head(ego)
barplot(ego)
goplot(ego)

### -- MF
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH")

head(ego)
barplot(ego)
goplot(ego)

### -- CC
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH")

head(ego)
barplot(ego)
goplot(ego)

## -- KEGG
kegg <- enrichKEGG(gene = entrez_genes,
                   organism = "dre",
                   pvalueCutoff = 0.05)

head(kegg)
head(kegg@result)
browseKEGG(kegg, "dre00100")
dotplot(kegg)
dotplot(ego,   showCategory = 20)
?dotplot

# -- Up -- #####################################################################
data.genes <- read.csv("../differential_expression/allSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")

data <- data.genes[data.genes$log2FoldChange > 0.5,]
data <- data[!is.na(data$padj),]
data <- data[data$padj < 0.05,]
nrow(data)

gene <- data$X

entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
entrez_genes <- entrez_genes[!is.na(entrez_genes)]

## -- BP
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH")

head(ego)
barplot(ego)
goplot(ego)

## -- MF
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH")

head(ego)
barplot(ego)
goplot(ego)

## -- CC
ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH")

head(ego)
barplot(ego)
goplot(ego)

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