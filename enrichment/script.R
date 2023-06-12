library(clusterProfiler)
library(org.Dr.eg.db)

data.genes <- read.csv("../differential_expression/AHuncut0_CTuncut0_res05_sig_fc0.csv")

gene.ahctuncut <- data.genes$X

ggo <- groupGO(gene     = gene.ahctuncut,
               OrgDb    = org.Dr.eg.db,
               keyType = "ENSEMBL",
               ont      = "CC",
               #level    = 3,
               readable = TRUE)
head(ggo)

entrez_genes <- mapIds(org.Dr.eg.db, gene.ahctuncut, 'ENTREZID', 'ENSEMBL')

ego <- enrichGO(gene         = entrez_genes,
                OrgDb         = org.Dr.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH"
                )
goplot(ego)
head(ego)



# -- KEEG
entrez_genes <- mapIds(org.Dr.eg.db, gene.ahctuncut, 'ENTREZID', 'ENSEMBL')

kegg <- enrichKEGG(gene = entrez_genes,
                   organism = "dre",
                   pvalueCutoff = 0.05)

head(kegg)


browseKEGG(kegg, "dre00100")