################################################################################
## RNA-Seq - zebrafish
## GO enrichment
## date of creation: 2023/06/12
## date of last update: 2023/08/01
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
## All Saamples                                                               ##
################################################################################
# -- down -- ###################################################################
fileNames <- list.files("../differential_expression/allSamples/")[1:7]

# vou precisar fazer destaa forma pois com o for normal estava dando erro nos
# ultimos dois arquivos

#i = fileNames[1]
#i = fileNames[2]
#i = fileNames[3]
#i = fileNames[4]
#i = fileNames[5]
#i = fileNames[6]
#i = fileNames[7]

data.genes <- read.csv(paste("../differential_expression/allSamples",
                             i, sep = "/"))
ego <- downEnrichGO(data.genes) # down
  
graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
  facet_grid(~ONTOLOGY)+
  labs(title = "GO Analysis",
       subtitle = "Top 10 terms for BP, CC and MF")
  
ggsave(paste("allSamples/down_GO_analysis_",gsub(".csv","",i), ".pdf",
             sep = ""), graph, width = 7, height = 5)
  
write.csv(ego@result, file=paste("allSamples/down_GO_analysis_",i,sep = ""))
  
kegg <- downKegg(data.genes) # down

graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
ggsave(paste("allSamples/down_KEGG_analysis_",gsub(".csv","",i),
             ".pdf",sep = ""), graph, width = 6, height = 6)

write.csv(kegg@result, file=paste("allSamples/down_KEGG_analysis_",i,sep = ""))

rm(ego,graph,kegg,fileNames,i,data.genes)

# -- up -- #####################################################################
fileNames <- list.files("../differential_expression/allSamples/")[1:7]

# vou precisar fazer destaa forma pois com o for normal estava dando erro nos
# ultimos dois arquivos

#i = fileNames[1]
#i = fileNames[2]
#i = fileNames[3]
#i = fileNames[4]
#i = fileNames[5]
#i = fileNames[6]
#i = fileNames[7]

data.genes <- read.csv(paste("../differential_expression/allSamples",
                             i, sep = "/"))
ego <- upEnrichGO(data.genes) # up

graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
  facet_grid(~ONTOLOGY)+
  labs(title = "GO Analysis",
       subtitle = "Top 10 terms for BP, CC and MF")

ggsave(paste("allSamples/up_GO_analysis_",gsub(".csv","",i), ".pdf",
             sep = ""), graph, width = 7, height = 5)

write.csv(ego@result, file=paste("allSamples/up_GO_analysis_",i,sep = ""))

kegg <- upKegg(data.genes) # up

graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
ggsave(paste("allSamples/up_KEGG_analysis_",gsub(".csv","",i),
             ".pdf",sep = ""), graph, width = 6, height = 6)

write.csv(kegg@result, file=paste("allSamples/up_KEGG_analysis_",i,sep = ""))

rm(ego,graph,kegg,fileNames,i,data.genes)
################################################################################
## Filt Saamples                                                               ##
################################################################################
# -- down -- ###################################################################
fileNames <- list.files("../differential_expression/filtSamples/")[1:7]

# vou precisar fazer destaa forma pois com o for normal estava dando erro nos
# ultimos dois arquivos

#i = fileNames[1]
#i = fileNames[2]
#i = fileNames[3]
#i = fileNames[4]
#i = fileNames[5]
#i = fileNames[6]
#i = fileNames[7]

data.genes <- read.csv(paste("../differential_expression/filtSamples",
                             i, sep = "/"))
ego <- downEnrichGO(data.genes) # down

graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
  facet_grid(~ONTOLOGY)+
  labs(title = "GO Analysis",
       subtitle = "Top 10 terms for BP, CC and MF")

ggsave(paste("filtSamples/down_GO_analysis_",gsub(".csv","",i), ".pdf",
             sep = ""), graph, width = 7, height = 5)

write.csv(ego@result, file=paste("filtSamples/down_GO_analysis_",i,sep = ""))

kegg <- downKegg(data.genes) # down

graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
ggsave(paste("filtSamples/down_KEGG_analysis_",gsub(".csv","",i),
             ".pdf",sep = ""), graph, width = 6, height = 6)

write.csv(kegg@result, file=paste("filtSamples/down_KEGG_analysis_",i,sep = ""))

rm(ego,graph,kegg,fileNames,i,data.genes)

# -- up -- #####################################################################
fileNames <- list.files("../differential_expression/filtSamples/")[1:7]

# vou precisar fazer destaa forma pois com o for normal estava dando erro nos
# ultimos dois arquivos

#i = fileNames[1]
#i = fileNames[2]
#i = fileNames[3]
#i = fileNames[4]
#i = fileNames[5]
#i = fileNames[6]
#i = fileNames[7]

data.genes <- read.csv(paste("../differential_expression/filtSamples",
                             i, sep = "/"))
ego <- upEnrichGO(data.genes) # up

graph <- dotplot(ego, split="ONTOLOGY", font.size=8, showCategory=5)+
  facet_grid(~ONTOLOGY)+
  labs(title = "GO Analysis",
       subtitle = "Top 10 terms for BP, CC and MF")

ggsave(paste("filtSamples/up_GO_analysis_",gsub(".csv","",i), ".pdf",
             sep = ""), graph, width = 7, height = 5)

write.csv(ego@result, file=paste("filtSamples/up_GO_analysis_",i,sep = ""))

kegg <- upKegg(data.genes) # up

graph <- dotplot(kegg, showCategory=nrow(kegg), font.size=8)
ggsave(paste("filtSamples/up_KEGG_analysis_",gsub(".csv","",i),
             ".pdf",sep = ""), graph, width = 6, height = 6)

write.csv(kegg@result, file=paste("filtSamples/up_KEGG_analysis_",i,sep = ""))

rm(ego,graph,kegg,fileNames,i,data.genes)
################################################################################