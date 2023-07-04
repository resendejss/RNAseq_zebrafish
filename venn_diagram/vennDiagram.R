library(VennDetail)

ahuncut0_ctuncut0 <- read.csv("../differential_expression/allSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")
ahcut0_ctcut0 <- read.csv("../differential_expression/allSamples/AHcut0_CTcut0_res05_sig_fc0.csv")
ctuncut0_ctcut0 <- read.csv("../differential_expression/allSamples/CTuncut0_CTcut0_res05_sig_fc0.csv")

gene_ahuncut0_ctuncut0 <- ahuncut0_ctuncut0[
  ahuncut0_ctuncut0$log2FoldChange < (-0.5) | ahuncut0_ctuncut0$log2FoldChange > 0.5,]
gene_ahuncut0_ctuncut0 <- gene_ahuncut0_ctuncut0[gene_ahuncut0_ctuncut0$padj < 0.05,]
gene_ahuncut0_ctuncut0 <- gene_ahuncut0_ctuncut0[!is.na(gene_ahuncut0_ctuncut0$padj),]
gene_ahuncut0_ctuncut0 <- gene_ahuncut0_ctuncut0$X

gene_ahcut0_ctcut0 <- ahcut0_ctcut0[ahcut0_ctcut0$log2FoldChange < (-0.5) | ahcut0_ctcut0$log2FoldChange > 0.5,]
gene_ahcut0_ctcut0 <- gene_ahcut0_ctcut0[gene_ahcut0_ctcut0$padj < 0.05,]
gene_ahcut0_ctcut0 <- gene_ahcut0_ctcut0[!is.na(gene_ahcut0_ctcut0$padj),]
gene_ahcut0_ctcut0 <- gene_ahcut0_ctcut0$X

gene_ctuncut0_ctcut0 <- ctuncut0_ctcut0[
  ctuncut0_ctcut0$log2FoldChange < (-0.5) | ctuncut0_ctcut0$log2FoldChange > 0.5,]
gene_ctuncut0_ctcut0 <- gene_ctuncut0_ctcut0[gene_ctuncut0_ctcut0$padj < 0.05,]
gene_ctuncut0_ctcut0 <- gene_ctuncut0_ctcut0[!is.na(gene_ctuncut0_ctcut0$padj),]
gene_ctuncut0_ctcut0 <- gene_ctuncut0_ctcut0$X

ven <- venndetail(list(ahuncut0_ctuncut0 = gene_ahuncut0_ctuncut0,
                       ahcut0_ctcut0 = gene_ahcut0_ctcut0,
                       ctuncut0_ctcut0 = gene_ctuncut0_ctcut0))
plot(ven)
plot(ven, type = "vennpie")
plot(ven, type = "upset")
