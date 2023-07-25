################################################################################
## RNA-Seq - zebrafish
## GO enrichment
## date of creation: 2023/06/12
## date of last update: 2023/06/16
## jean resende
################################################################################
library(clusterProfiler)
library(ggplot2)

# -- Down -- ###################################################################
## -- GO

## -- functions
downEnrichGO <- function(data.genes){
  require(clusterProfiler)
  require(org.Dr.eg.db)
  data <- data.genes[data.genes$log2FoldChange < -0.5,]
  data <- data[!is.na(data$padj),]
  data <- data[data$padj < 0.05,]
  nrow(data)
  
  gene <- data$X
  
  entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
  entrez_genes <- entrez_genes[!is.na(entrez_genes)]
  
  ego <- enrichGO(gene         = entrez_genes,
                  OrgDb         = org.Dr.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH")
  return(ego)
}

upEnrichGO <- function(data.genes){
  require(clusterProfiler)
  require(org.Dr.eg.db)
  data <- data.genes[data.genes$log2FoldChange > 0.5,]
  data <- data[!is.na(data$padj),]
  data <- data[data$padj < 0.05,]
  nrow(data)
  
  gene <- data$X
  
  entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
  entrez_genes <- entrez_genes[!is.na(entrez_genes)]
  
  ego <- enrichGO(gene         = entrez_genes,
                  OrgDb         = org.Dr.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH")
  return(ego)
}

downKegg <- function(data.genes){
  require(clusterProfiler)
  require(org.Dr.eg.db)
  data <- data.genes[data.genes$log2FoldChange < -0.5,]
  data <- data[!is.na(data$padj),]
  data <- data[data$padj < 0.05,]
  nrow(data)
  
  gene <- data$X
  
  entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
  entrez_genes <- entrez_genes[!is.na(entrez_genes)]
  
  kegg <- enrichKEGG(gene = entrez_genes,
                     organism = "dre",
                     pvalueCutoff = 0.05)
  return(kegg)
}

upKegg <- function(data.genes){
  require(clusterProfiler)
  require(org.Dr.eg.db)
  data <- data.genes[data.genes$log2FoldChange > 0.5,]
  data <- data[!is.na(data$padj),]
  data <- data[data$padj < 0.05,]
  nrow(data)
  
  gene <- data$X
  
  entrez_genes <- mapIds(org.Dr.eg.db, gene, 'ENTREZID', 'ENSEMBL')
  entrez_genes <- entrez_genes[!is.na(entrez_genes)]
  
  kegg <- enrichKEGG(gene = entrez_genes,
                     organism = "dre",
                     pvalueCutoff = 0.05)
  return(kegg)
}

################################################################################
# -- down
fileNames <- list.files("../differential_expression/allSamples/")

for (i in 1:length(fileNames)) {
  data.genes <- read.csv(paste("../differential_expression/allSamples",
                               fileNames[i], sep = "/"))
  
  ego <- downEnrichGO(data.genes) # down
  
  graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
    facet_grid(~ONTOLOGY)+
    labs(title = "GO Analysis",
         subtitle = "Top 10 terms for BP, CC and MF")
  
  ggsave(paste("allSamples/down_GO_analysis_",gsub(".csv","",fileNames[i]), ".pdf",
               sep = ""), graph, width = 7, height = 5)
  
  write.csv(ego@result, file=paste("allSamples/down_GO_analysis_",fileNames[i],
                                   sep = ""))
  
  kegg <- downKegg(data.genes) # down

  graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
  ggsave(paste("allSamples/down_KEGG_analysis_",gsub(".csv","",fileNames[i]),
               ".pdf",sep = ""), graph, width = 6, height = 6)
}

# -- up
fileNames <- list.files("../differential_expression/allSamples/")

for (i in 1:length(fileNames)) {
  data.genes <- read.csv(paste("../differential_expression/allSamples",
                               fileNames[i], sep = "/"))
  
  ego <- upEnrichGO(data.genes) # down
  
  graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
    facet_grid(~ONTOLOGY)+
    labs(title = "GO Analysis",
         subtitle = "Top 10 terms for BP, CC and MF")
  
  ggsave(paste("allSamples/up_GO_analysis_",gsub(".csv","",fileNames[i]), ".pdf",
               sep = ""), graph, width = 7, height = 5)
  
  write.csv(ego@result, file=paste("allSamples/up_GO_analysis_",fileNames[i],
                                   sep = ""))
  
  kegg <- upKegg(data.genes) # down
  
  graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
  ggsave(paste("allSamples/up_KEGG_analysis_",gsub(".csv","",fileNames[i]),
               ".pdf",sep = ""), graph, width = 6, height = 6)
}

