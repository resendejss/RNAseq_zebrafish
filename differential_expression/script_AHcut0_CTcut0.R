################################################################################
## EXPRESSAO DIFERENCIAL 
## Bioinformatica Para Novatos
## Jean Resende
################################################################################
library(DESeq2)
library(readxl)

# -- genes codificantes -- 
## -- preparacao do objeto
load("../PreProcSEQ-main/5-expressionMatrix/tximport/mat_gse_filt_20230519.RData")
coldata <- read_xlsx("../PreProcSEQ-main/coldata.xlsx")
load("../PreProcSEQ-main/6-annotationTranscripts/EnsDbAnnotation_20230519_atual.RData")

coldata$Trat_01 <- as.factor(coldata$Trat_01) 
coldata$Trat_02 <- as.factor(coldata$Trat_02)
coldata$Trat_03 <- as.factor(coldata$Trat_03)

data.protein.coding <- EnsDbAnnotation[
  EnsDbAnnotation$gene_biotype == "protein_coding",]

idx <- match(rownames(data.protein.coding),rownames(mat_gse_filt$counts))
idx <- idx[!is.na(idx)]

mat_gse.ahct <- mat_gse_filt

colnames(mat_gse.ahct$abundance)

mat_gse.ahct$abundance <- mat_gse.ahct$abundance[idx,c(1:3,16:19)]
mat_gse.ahct$counts <- mat_gse.ahct$counts[idx,c(1:3,16:19)]
mat_gse.ahct$length <- mat_gse.ahct$length[idx,c(1:3,16:19)]

coldata.ahct <- coldata[c(1:3,17:20),1:2]

ddsTxi <- DESeqDataSetFromTximport(mat_gse.ahct,
                                   colData = coldata.ahct,
                                   design = ~ Trat_01)

ddsTxi$names

## -- pre-filtragem
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

## -- definindo a referencia
dds$Trat_01 <- relevel(ddsTxi$Trat_01, ref = "AH")

## -- expressao diferencial
dds <- DESeq(dds)
res <- results(dds, contrast = c("Trat_01","AH","CT"))

res_padj05_lfc0 <- results(dds, contrast = c("Trat_01","AH","CT"), alpha = 0.05)
#res_padj05_lfc1 <- results(dds,
#                           contrast = c("Trat_01","AH","CT"),
#                           alpha = 0.05,
#                           lfcThreshold = 1)

summary(res_padj05_lfc0)
summary(res)
#summary(res_padj05_lfc1)

res_padj05_lfc0.df <- res_padj05_lfc0[!is.na(res_padj05_lfc0$padj) &
                                        res_padj05_lfc0$padj < 0.05 & 
                                        abs(res_padj05_lfc0$log2FoldChange) > 1,]
summary(res)
summary(res_padj05_lfc1)

write.csv(res_padj05_lfc0.df, file = "AHcut0_CTcut0_res05.sig.fc1.csv")

# -- PLOTS
library(ggplot2)
library(ggrepel)

data <- as.data.frame(res_padj05_lfc0)

names.diff.up <- rownames(res_padj05_lfc0.df[res_padj05_lfc0.df$log2FoldChange > 1,])
names.diff.down <- rownames(res_padj05_lfc0.df[res_padj05_lfc0.df$log2FoldChange < -1,])

idx.up <- match(names.diff.up,rownames(data)) 
idx.down <- match(names.diff.down,rownames(data))

#idx.up <- match(names.diff.up,rownames(data))
#idx.down <- match(names.diff.down, rownames(data))

data$diffexpressed <- "NO"
data$diffexpressed[idx.up] <- "UP"
data$diffexpressed[idx.down] <- "DOWN"

data$delabel <- NA
data$delabel[data$diffexpressed!="NO"] <- 
  rownames(data[data$diffexpressed != "NO",])

mycolors <- c("blue","red","gray")
names(mycolors) <- c("DOWN","UP","NO")

pdf("vulcanoplot_AHcut0_CTcut0_fc1.pdf", width = 8, height = 6)
p <- ggplot(data = data,
             aes(x=log2FoldChange, y= -log10(padj), col=diffexpressed, label=delabel))+
  geom_point()+
  theme_minimal()+
  geom_text_repel()+
  #geom_vline(xintercept = c(-1,1), col="red")+
  #geom_hline(yintercept = -log10(0.05), col="red")+
  scale_color_manual(values = mycolors)+
  xlim(-10,10)
p
dev.off()

library(Glimma)
glimmaMDS(dds)
#glimmaMA(dds)

library(pcaExplorer)
pcaExplorer(dds)

################################################################################
# -- enriquecimento
library("AnnotationDbi")
library(org.Dr.eg.db)

geneUniverse <- rownames(res_padj05_lfc0.df)

ids <- mapIds(org.Dr.eg.db,
              keys=row.names(res_padj05_lfc0.df),
              column="ENSEMBL",
              keytype="ENSEMBL",
              multiVals="first")

#ids <- ids[!is.na(ids)]

#idx <- match(names(ids), rownames(res_padj05_lfc1.df))

#deGenes <- res_padj05_lfc1.df[idx,]
#rownames(deGenes) <- ids

deGenes <- res_padj05_lfc0.df

library(clusterProfiler)
ans.go <- enrichGO(gene = rownames(deGenes), ont = "ALL",
                   OrgDb ="org.Dr.eg.db",
                   universe = ids,
                   readable=TRUE,
                   pvalueCutoff = 0.05)


ans.kegg <- enrichKEGG(gene = rownames(deGenes),
                       organism = 'dre',
                       #universe = ids,
                       pvalueCutoff = 0.05)

tab.go <- as.data.frame(ans.go)
tab.go<- subset(tab.go, Count>5)
tab.go[1:5, 1:6]


