################################################################################
## main components
## jean resende
## RNAseq_zebrafish
################################################################################
library(magrittr)
library(DESeq2)
library(SummarizedExperiment)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(pcaExplorer)

# -- coldata --
colData <- readxl::read_xlsx("coldata.xlsx")

colData$Trat_01 <- as.factor(colData$Trat_01)
colData$Trat_02 <- as.factor(colData$Trat_02)
colData$Trat_03 <- as.factor(colData$Trat_03)

# -- matrix of counts --
load("matrix_salmon_tximport_20230519.RData")
mat <- mat_gse$counts %>% round()
head(mat)

# -- dds object --
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = colData,
                              design = ~ Trat_01 + Trat_02 + Trat_03)

# -- vst normalization
vsd <- vst(dds, blind = FALSE)

# -- samples distances --
## -- euclidean method
sampleDists <- dist(t(assay(vsd)))
sampleDistsMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(255)
pdf("sampleDistsMatrix_euclidean.pdf", height = 7, width = 7)
pheatmap(sampleDistsMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()

## -- poisson method
sampleDistPoisson <- PoissonDistance(t(counts(dds)))
sampleDistPoissonMatrix <- as.matrix(sampleDistPoisson$dd)
rownames(sampleDistPoissonMatrix) <- colnames(counts(dds))
colnames(sampleDistPoissonMatrix) <- colnames(counts(dds))
pheatmap(sampleDistPoissonMatrix,
         clustering_distance_rows = sampleDistPoisson$dd,
         clustering_distance_cols = sampleDistPoisson$dd,
         col = colors)
pcaExplorer(dds)
