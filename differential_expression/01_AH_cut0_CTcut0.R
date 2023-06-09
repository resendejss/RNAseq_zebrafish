################################################################################
## differential expression
## jean resende
## RNAseq_zebrafish
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

mat.gse.ahctCut0 <- mat_gse

colnames(mat.gse.ahctCut0$counts)

mat.gse.ahctCut0$abundance <- mat.gse.ahctCut0$abundance[idx,c(1:4,17:20)]
mat.gse.ahctCut0$counts <- mat.gse.ahctCut0$counts[idx,c(1:4,17:20)]
mat.gse.ahctCut0$length <- mat.gse.ahctCut0$length[idx,c(1:4,17:20)]

coldata.ahct <- coldata[c(1:4,17:20),1:2]

ddsTxi <- DESeqDataSetFromTximport(mat.gse.ahctCut0,
                                   colData = coldata.ahct,
                                   design = ~ Trat_01)

# pre-filtration
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

# setting the reference
dds$Trat_01 <- relevel(ddsTxi$Trat_01, ref = "AH")

# differential expression
dds <- DESeq(dds)
res <- results(dds, contrast = c("Trat_01","AH","CT"))
res.padj05.lfc0 <- results(dds, contrast = c("Trat_01","AH","CT"), alpha = 0.05)

#summary(res)
summary(res.padj05.lfc0)

write.csv(res.padj05.lfc0, file = "AHcut0_CTcut0_res05_sig_fc0.csv")
