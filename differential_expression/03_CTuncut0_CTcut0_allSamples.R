################################################################################
# Name: 03_CTuncut0_CTcut0_allSamples                                          #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/06/27                                                    #
# Last update date: 2023/06/27                                                 #
# Description: Differential expression                                         #
################################################################################
library(DESeq2)
library(readxl)

# building object
load("matrix_salmon_tximport_20230519.RData")
load("EnsDbAnnotation_20230519_atual.RData")
coldata <- read_xlsx("coldata.xlsx")

coldata$Trat_01 <- as.factor(coldata$Trat_01)
coldata$Trat_02 <- as.factor(coldata$Trat_02)
coldata$Trat_03 <- as.factor(coldata$Trat_03)

data.prot.cod <- EnsDbAnnotation[
  EnsDbAnnotation$gene_biotype == "protein_coding",]

head(rownames(mat_gse$counts))
head(rownames(data.prot.cod))

idx <- match(rownames(data.prot.cod), rownames(mat_gse$counts))
idx <- idx[!is.na(idx)]

mat.gse.ctCutUncut0 <- mat_gse

colnames(mat.gse.ctCutUncut0$counts)

mat.gse.ctCutUncut0$abundance <- mat.gse.ctCutUncut0$abundance[idx,c(17:20,29:32)]
mat.gse.ctCutUncut0$counts <- mat.gse.ctCutUncut0$counts[idx,c(17:20,29:32)]
mat.gse.ctCutUncut0$length <- mat.gse.ctCutUncut0$length[idx,c(17:20,29:32)]

colnames(mat.gse.ctCutUncut0$counts)

coldata.ct <- coldata[coldata$names %in% colnames(mat.gse.ctCutUncut0$counts),c(1,3)]

ddsTxi <- DESeqDataSetFromTximport(mat.gse.ctCutUncut0,
                                   colData = coldata.ct,
                                   design = ~ Trat_02)

# pre-filtration
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

# setting the reference
dds$Trat_02 <- relevel(ddsTxi$Trat_02, ref = "uncut")

# differential expression
dds <- DESeq(dds)
res <- results(dds, contrast = c("Trat_02","uncut","cut"))
res.padj05.lfc0 <- results(dds, contrast = c("Trat_02","uncut","cut"), alpha = 0.05)

#summary(res)
summary(res.padj05.lfc0)

write.csv(res.padj05.lfc0, file = "allSamples/CTuncut0_CTcut0_res05_sig_fc0.csv")
################################################################################
# plot
library(ggplot2)
library(ggrepel)

data <- as.data.frame(res.padj05.lfc0)

data2 <- data[data$padj < 0.05,]
data2 <- data2[!is.na(data2$padj),]

data2 <- data2[data2$log2FoldChange > 1 | data2$log2FoldChange < -1,]

names.diff.up <- rownames(data2[data2$log2FoldChange > 1,])
names.diff.down <- rownames(data2[data2$log2FoldChange < -1,])

idx.up <- match(names.diff.up, rownames(data))
idx.down <- match(names.diff.down, rownames(data))

data$diffexpressed <- "NO"
data$diffexpressed[idx.up] <- "UP"
data$diffexpressed[idx.down] <- "DOWN"

data$delabel <- NA
data$delabel[data$diffexpressed!="NO"] <- rownames(
  data[data$diffexpressed != "NO",])

mycolors <- c("#00AFBB","#bb0c00","grey")
names(mycolors) <- c("DOWN","UP","NO")

pdf("filtSamples/vulcanoplot_AHcut0_CTcut0_fc1.pdf", width = 8, height = 6)
ggplot(data = data,
       aes(x=log2FoldChange, y= -log10(padj), col=diffexpressed, label=delabel))+
  geom_point(size=0.5)+
  theme_minimal()+
  geom_text_repel()+
  #geom_vline(xintercept = c(-1,1), col="red")+
  #geom_hline(yintercept = -log10(0.05), col="red")+
  scale_color_manual(values = mycolors)+
  xlim(-10,10)
dev.off()
