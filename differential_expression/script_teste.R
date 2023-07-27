################################################################################
# Name: differentialExpression                                                 #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/06/28                                                    #
# Last update date: 2023/07/27                                                 #
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
ahcUncut0_filtSamples <- grep(".._Uncut_0_.", colnames(data_filtSamples$counts)) # ahcUncut0_filtSamples

ctCutUncut0_allSamples <- grep("^CT.*_0_.", colnames(data_allSamples$counts)) # ctCutUncut0_allSamples
ctCutUncut0_filtSamples <- grep("^CT.*_0_.", colnames(data_filtSamples$counts)) # ctCutUncut0_filtSamples

ahctCut5_allSamples <- grep(".._Cut_5_.", colnames(data_allSamples$counts)) # ahctCut5_allSamples
ahctCut5_filtSamples <- grep(".._Cut_5_.", colnames(data_filtSamples$counts)) # ahctCut5_filtSamples

ahctCut100_allSamples <- grep(".._Cut_100_.", colnames(data_allSamples$counts)) # ahctCut100_allSamples
ahctCut100_filtSamples <- grep(".._Cut_100_.", colnames(data_filtSamples$counts)) # ahctCut100_filtSamples

ctCut0100_allSamples <- grep("CT_Cut_.*0_.", colnames(data_allSamples$counts)) # ctCut0100_allSamples
ctCut0100_filtSamples <- grep("CT_Cut_.*0_.", colnames(data_filtSamples$counts)) # ctCut0100_filtSamples

ctCut05_allSamples <- grep("CT_Cut_._.", colnames(data_allSamples$counts)) # ctCut05_allSamples
ctCut05_filtSamples <- grep("CT_Cut_._.", colnames(data_filtSamples$counts)) # ctCut05_filtSamples

names_allSamples <- ls()[grep(".allSamples", ls())]
names_allSamples <- names_allSamples["data_allSamples"!= names_allSamples]
#names_allSamples <- as.factor(names_allSamples)

names_filtSamples <- ls()[grep(".filtSamples", ls())]
names_filtSamples <- names_filtSamples["data_filtSamples"!= names_filtSamples]
#names_filtSamples <- as.factor(names_filtSamples)

for (i in names_allSamples) {
  data <- data_allSamples
  data$abundance <- data$abundance[idx,ls()[i]]
  data$counts <- data$counts[idx,ls()[i]]
  data$length <- data$length[idx,ls()[i]]
  
  # -- coldata
  if(i == ahctCut0_allSamples | i == ahctCut100_allSamples |
     i == ahctCut5_allSamples | i == ahctUncut0_allSamples){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_01")]
    
  } else if(i == ctCut0100_allSamples | i == ctCut05_allSamples){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_03")]
  
  } else if(i == ctCutUncut0_allSamples){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_02")]
  }else{
  cat("verifique a comparação 1")
  }
 
  # -- preparacao do objeto
  if(i == ahctCut0_allSamples | i == ahctCut100_allSamples |
     i == ahctCut5_allSamples | i == ahctUncut0_allSamples){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_01)
    
  } else if(i == ctCut0100_allSamples | i == ctCut05_allSamples){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_03)
    
  } else if(i == ctCutUncut0_allSamples){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_02)
  }else{
    cat("verifique a comparação 2")
  }
  
  
  # -- pre-filtragem
  keep <- rowSums(counts(ddsTxi)) >= 10
  dds <- ddsTxi[keep,]
  
  
  # -- setting the reference
  if(i == ahctCut0_allSamples | i == ahctCut100_allSamples |
     i == ahctCut5_allSamples | i == ahctUncut0_allSamples){
    dds$Trat_01 <- relevel(dds$Trat_01, ref = "CT")

  } else if(i == ctCut0100_allSamples | i == ctCut05_allSamples){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = 0)
    
  } else if(i == ctCutUncut0_allSamples){
    dds$Trat_02 <- relevel(dds$Trat_02, ref = "uncut")
  }else{
    cat("verifique a comparação 3")
  }

  # -- differential expression
  if(i == ahctCut0_allSamples | i == ahctCut100_allSamples |
     i == ahctCut5_allSamples | i == ahctUncut0_allSamples){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_01","AH","CT"), alpha = 0.05)
    
  } else if(i == ctCut0100_allSamples | i == ctCut05_allSamples){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03",100,0), alpha = 0.05)
    
  } else if(i == ctCutUncut0_allSamples){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_02","cut","uncut"), alpha = 0.05)
  }else{
    cat("verifique a comparação 4")
  }
  
  write.csv(res_padj05_lfc0,
            file = paste("allSamples/",i,"_res05_sig_fc0.csv", sep = ""))
}

