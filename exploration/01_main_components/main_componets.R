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

dds <- DESeqDataSetFromTximport(mat_gse,
                                colData = colData,
                                design = ~ Trat_01 + Trat_02 + Trat_03)

# -- dds object --
#dds <- DESeqDataSetFromMatrix(countData = mat,
#                              colData = colData,
#                              design = ~ Trat_01 + Trat_02 + Trat_03)

# -- vst normalization
vsd <- vst(dds, blind = FALSE)

# -- samples distances --
## -- euclidean method
sampleDists <- dist(t(assay(vsd)))
sampleDistsMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(255)

pdf("sampleDistsMatrix_euclidean_all.pdf", height = 7, width = 7)
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

pdf("sampleDistsMatrix_poisson_all.pdf", height = 7, width = 7)
pheatmap(sampleDistPoissonMatrix,
         clustering_distance_rows = sampleDistPoisson$dd,
         clustering_distance_cols = sampleDistPoisson$dd,
         col = colors)
dev.off()
################################################################################
library(factoextra)
library(FactoMineR)
library(readxl)

# -- matrix of counts --
load("matrix_salmon_tximport_20230519.RData")
coldata <- read_xlsx("coldata.xlsx")

coldata$Trat_01 <- as.factor(coldata$Trat_01)
coldata$Trat_02 <- as.factor(coldata$Trat_02)
coldata$Trat_03 <- as.factor(coldata$Trat_03)

mat <- mat_gse$counts
head(mat)
mat.t <- t(mat)

res.pca <- PCA(mat.t, graph = FALSE)

#fviz_eig(res.pca, addlabels=TRUE, ylim = c(0,50))

all(rownames(mat.t)==coldata$names)

pdf("pca_all_AHCT.pdf", width = 8, height = 6)
fviz_pca_ind(res.pca,
             label = "none",
             addEllipses = TRUE,
             habillage = coldata$Trat_01,
             ellipse.type = "confidence",
             ellipse.level=0.95,
             repel = TRUE)
dev.off()

mat.AH <- mat[,1:16]
mat.AH.t <- t(mat.AH)

res.pca <- PCA(mat.AH.t, graph = FALSE)

#fviz_eig(res.pca, addlabels=TRUE, ylim = c(0,50))

coldata.AH <- coldata[1:16,]
all(rownames(mat.AH.t)==coldata.AH$names)

pdf("pca_AH_cutUncut.pdf", width = 8, height = 6)
fviz_pca_ind(res.pca,
             label = "none",
             addEllipses = TRUE,
             habillage = coldata.AH$Trat_03,
             ellipse.type = "confidence",
             ellipse.level=0.95,
             repel = TRUE)
)
dev.off()

pdf("pca_AH_05100.pdf", width = 8, height = 6)
fviz_pca_ind(res.pca,
             label = "none",
             addEllipses = TRUE,
             habillage = coldata.AH$Trat_03,
             ellipse.type = "confidence",
             ellipse.level=0.95,
             repel = TRUE)
)
dev.off()

################################################################################


plot(res.pca,choix="ind",habillage=2)

pcaExplorer(dds)
