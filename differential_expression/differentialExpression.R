################################################################################
# Name: differentialExpression                                                 #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/06/28                                                    #
# Last update date: 2023/06/28                                                 #
# Description: Differential expression                                         #
################################################################################
library(DESeq2)
library(readxl)

# -- building object
load("matrix_salmon_tximport_20230519.RData")
load("mat_gse_filt_20230519.RData")
load("EnsDbAnnotation_20230519_atual.RData")
coldata <- read_xlsx("coldata.xlsx")

# -- arrumando a estrutura do objeto coldata
str(coldata)
coldata <- as.data.frame(coldata)

for (i in 2:ncol(coldata)){
  coldata[,i] <- as.factor(coldata[,i])
}

str(coldata)

# -- matrizes com allSamples e filtSamples
data_allSamples <- mat_gse; rm(mat_gse)
data_filtSamples  <- mat_gse_filt; rm(mat_gse_filt)

# -- somente genes codificantes de proteinas
data_protCod <- EnsDbAnnotation[EnsDbAnnotation$gene_biotype=="protein_coding",]
head(data_protCod)
table(data_protCod$gene_biotype) # 30153

head(rownames(data_allSamples$counts)) # sem num da versao
head(rownames(data_filtSamples$counts)) # sem num da versao
head(rownames(data_protCod)) # sem num da versao

nrow(data_allSamples$counts)
nrow(data_filtSamples$counts)
identical(rownames(data_allSamples$counts), rownames(data_filtSamples$counts))

idx <- match(rownames(data_protCod), rownames(data_allSamples$counts))
table(is.na(idx))
idx <- idx[!is.na(idx)]

colnames(data_allSamples$counts)
colnames(data_filtSamples$counts)

ahctCut0_allSamples <- grep(".._Cut_0_.", colnames(data_allSamples$counts)) # ahctCut0_allSamples
ahctCut0_filtSamples <- grep(".._Cut_0_.", colnames(data_filtSamples$counts)) # ahctCut0_filtSamples

ahctUncut0_allSamples <- grep(".._Uncut_0_.", colnames(data_allSamples$counts)) # ahctUncut0_allSamples
ahcUncut0_filtSamples <- grep(".._Uncut_0_.", colnames(data_allSamples$counts)) # ahcUncut0_filtSamples

ctCutUncut0_allSamples <- grep("^CT.*_0_.", colnames(data_allSamples$counts)) # ctCutUncut0_allSamples
ctCutUncut0_filtSamples <- grep("^CT.*_0_.", colnames(data_filtSamples$counts)) # ctCutUncut0_filtSamples

ahctCut5_allSamples <- grep(".._Cut_5_.", colnames(data_allSamples$counts)) # ahctCut5_allSamples
ahctCut5_filtSamples <- grep(".._Cut_5_.", colnames(data_filtSamples$counts)) # ahctCut5_filtSamples

ahctCut100_allSamples <- grep(".._Cut_100_.", colnames(data_allSamples$counts)) # ahctCut100_allSamples
ahctCut100_filtSamples <- grep(".._Cut_100_.", colnames(data_filtSamples$counts)) # ahctCut100_filtSamples

ctCut0100_allSamples <- grep("CT_Cut_.*0_.", colnames(data_allSamples$counts)) # ctCut0100_allSamples
ctCut0100_filtSamples <- grep("CT_Cut_.*0_.", colnames(data_filtSamples$counts)) # ctCut0100_filtSamples

################################################################################
data <- data_allSamples
data$abundance <- data$abundance[idx,ahctCut0_allSamples]
data$counts <- data$counts[idx,ahctCut0_allSamples]
data$length <- data$length[idx,ahctCut0_allSamples]

coldata_samples <- coldata[ahctCut0_allSamples, c("names","Trat_01")]
################################################################################

ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                   colData = coldata_samples,
                                   design = ~ Trat_01)

# -- pre-filtragem
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]

# -- setting the reference
dds$Trat_01 < relevel(dds$Trat_01, ref = "AH")

# -- differential expression
dds <- DESeq(dds)
res <- results(dds, contrast = c("Trat_01","AH","CT"))
res_padj05_lfc0 <- results(dds, contrast = c("Trat_01","AH","CT"), alpha = 0.05)

summary(res)
summary(res_padj05_lfc0)

################################################################################
write.csv(res_padj05_lfc0, file = "allSamples/AHcut0_CTcut0_res05_sig_fc0.csv")
################################################################################