################################################################################
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(ggplot2)
library(readxl)

setwd("/media/jean/hd_biosys_1tb/bioinformatica_para_novatos/mentorias/RodrigoDisner/PreProcSEQ-main/5-expressionMatrix/tximport")

load("matrix_salmon_tximport_20230519.RData")
load("../../6-annotationTranscripts/EnsDbAnnotation_20230519_atual.RData")
coldata <- read_xlsx("../../coldata.xlsx")

coldata$Trat_01 <- as.factor(coldata$Trat_01) 
coldata$Trat_02 <- as.factor(coldata$Trat_02)
coldata$Trat_03 <- as.factor(coldata$Trat_03)

mat <- mat_gse$counts
mat <- round(mat)

dds.mat <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = coldata,
                                  design = ~ Trat_01 + Trat_02 + Trat_03)

# normalizacao vst
vsd <- vst(dds.mat, blind = FALSE)

# samples distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$Trat_01, vsd$Trat_02, sep = " - " )
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

table(vsd$Trat_01)
table(vsd$Trat_02)
table(vsd$Trat_03)

## method Poisson
poisd <- PoissonDistance(t(counts(dds.mat)))

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- colnames(counts(dds.mat))
colnames(samplePoisDistMatrix) <- colnames(counts(dds.mat))
#rownames(samplePoisDistMatrix) <- paste(dds.mat$Trat_01, dds.mat$Trat_02, sep=" - " )
#colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# -- PCA
## -- 01
plotPCA(vsd, intgroup = c("Trat_01", "Trat_02", "Trat_03"))

## -- 02
pcaData <- plotPCA(vsd, intgroup = c("Trat_01", "Trat_02", "Trat_03"),return=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Trat_03, shape = Trat_02)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


library(pcaExplorer)
pcaExplorer()
?pcaExplorer

pcaExplorer(dds = dds.mat)


mat_gse_filt <- mat_gse
idx <- match(c("CT_Cut_100_3","CT_Cut_5_3","AH_Cut_0_2"),
             colnames(mat_gse_filt$abundance))
mat_gse_filt$abundance <- mat_gse_filt$abundance[,-idx]
mat_gse_filt$counts <- mat_gse_filt$counts[,-idx]
mat_gse_filt$length <- mat_gse_filt$length[,-idx]

ncol(mat_gse_filt$abundance)

save(mat_gse_filt, file="mat_gse_filt_20230519.RData")
