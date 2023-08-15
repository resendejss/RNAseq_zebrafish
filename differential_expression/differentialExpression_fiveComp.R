################################################################################
# Name: differentialExpression_fiveComp                                        #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/08/15                                                    #
# Last update date: 2023/08/15                                                 #
# Description: Differential expression                                         #
################################################################################

# 1) AH_cut_0 vs AH_uncut_0 = Trat 02
# 2) CT_cut_100 vs CT_cut_5 = Trat_03
# 3) AH_cut_5 vs AH_cut_0 = Trat_03
# 4) AH_cut_100 vs AH_cut_0 = Trat_03
# 5) AH_cut_100 vs AH_cut_5 = Trat_03

######################## -- all samples -- #####################################

# -- funcoes -- ################################################################
DEGs_allSamples <- function(coldata, i, data){
  require(DESeq2)
  # -- coldata
  if(i == "ahCutUncut0_allSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_02")]
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut50_allSamples" |
            i == "ahCut1000_allSamples" | i == "ahCut1005_allSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_03")]
  }else{
    cat("verifique a comparação 1")
  }
  
  # -- preparacao do objeto
  if(i == "ahCutUncut0_allSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_02)
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut50_allSamples" |
            i == "ahCut1000_allSamples" | i == "ahCut1005_allSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_03)
  }else{
    cat("verifique a comparação 2")
  }
  
  
  # -- pre-filtragem
  keep <- rowSums(counts(ddsTxi)) >= 10
  dds <- ddsTxi[keep,]
  
    # -- setting the reference
  if(i == "ahCutUncut0_allSamples"){
    dds$Trat_02 <- relevel(dds$Trat_02, ref = "uncut")
    
  } else if(i == "ahCut50_allSamples" | i == "ahCut1000_allSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "0")
    
  } else if(i == "ctCut1005_allSamples" | i == "ahCut1005_allSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "5")
  }else{
    cat("verifique a comparação 3")
  }
  
  # -- differential expression
  if(i == "ahCutUncut0_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_02","cut","uncut"), alpha = 0.05)
    
  } else if(i == "ahCut50_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","5","0"), alpha = 0.05)
    
  }else if(i == "ahCut1000_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","0"), alpha = 0.05)
    
  } else if(i == "ctCut1005_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)
    
  }else if(i == "ahCut1005_allSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)

  }else{
    cat("verifique a comparação 4")
  }
  
  
  if(i == "ahCutUncut0_allSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/allSamples/ahCutUncut0_allSamples_res05_sig_fc0.csv")
  } else if(i == "ahCut50_allSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/allSamples/ahCut50_allSamples_res05_sig_fc0.csv")
  } else if(i == "ahCut1000_allSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/allSamples/ahCut1000_allSamples_res05_sig_fc0.csv")
  } else if(i == "ctCut1005_allSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/allSamples/ctCut1005_allSamples_res05_sig_fc0.csv")
  }else if(i == "ahCut1005_allSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/allSamples/ahCut1005_allSamples_res05_sig_fc0.csv")
  } else{
    cat("verifique a comparacao 5")
  }
}

################################################################################
# -- Anotacao -- 
load("EnsDbAnnotation_20230519_atual.RData")
data_protCod <- EnsDbAnnotation[EnsDbAnnotation$gene_biotype=="protein_coding",]
head(data_protCod)
table(data_protCod$gene_biotype) # 30153

# -- colldata -- 
coldata <- readxl::read_xlsx("coldata.xlsx")
str(coldata)
coldata <- as.data.frame(coldata)

for (i in 2:ncol(coldata)){
  coldata[,i] <- as.factor(coldata[,i])
}

str(coldata)

# -- preparando objeto -- ######################################################
load("matrix_salmon_tximport_20230519.RData")

## -- matrizes com allSamples 
data_allSamples <- mat_gse; rm(mat_gse)

head(rownames(data_allSamples$counts)) # sem num da versao
head(rownames(data_protCod)) # sem num da versao

nrow(data_allSamples$counts)

idx <- match(rownames(data_protCod), rownames(data_allSamples$counts))
table(is.na(idx))
idx <- idx[!is.na(idx)]

colnames(data_allSamples$counts)

ahCutUncut0_allSamples <- grep("^AH.*_0_.", colnames(data_allSamples$counts)) # ahCutUncut0_allSamples

ahCut50_allSamples <- grep("AH_Cut_._.", colnames(data_allSamples$counts)) # ahCut50_allSamples

ahCut1000_allSamples <- grep("AH_Cut_.*0_.", colnames(data_allSamples$counts)) # ahCut1000_allSamples

ctCut1005_allSamples <- c(grep("CT_Cut_100_.", colnames(data_allSamples$counts)),
                          grep("CT_Cut_5_.", colnames(data_allSamples$counts)))# ctCut1005_allSamples

ahCut1005_allSamples <- c(grep("AH_Cut_100_.", colnames(data_allSamples$counts)),
                          grep("AH_Cut_5_.", colnames(data_allSamples$counts)))# ahCut1005_allSamples

# -- expressao diferencial: all Samples -- #####################################

## -- ahCutUncut0_allSamples -- 
data <- data_allSamples
i="ahCutUncut0_allSamples"
data$abundance <- data$abundance[idx,ahCutUncut0_allSamples]
data$counts <- data$counts[idx,ahCutUncut0_allSamples]
data$length <- data$length[idx,ahCutUncut0_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut50_allSamples
data <- data_allSamples
i="ahCut50_allSamples"
data$abundance <- data$abundance[idx,ahCut50_allSamples]
data$counts <- data$counts[idx,ahCut50_allSamples]
data$length <- data$length[idx,ahCut50_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut1000_allSamples
data <- data_allSamples
i="ahCut1000_allSamples"
data$abundance <- data$abundance[idx,ahCut1000_allSamples]
data$counts <- data$counts[idx,ahCut1000_allSamples]
data$length <- data$length[idx,ahCut1000_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ctCut1005_allSamples
data <- data_allSamples
i="ctCut1005_allSamples"
data$abundance <- data$abundance[idx,ctCut1005_allSamples]
data$counts <- data$counts[idx,ctCut1005_allSamples]
data$length <- data$length[idx,ctCut1005_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

## -- ahCut1005_allSamples
data <- data_allSamples
i="ahCut1005_allSamples"
data$abundance <- data$abundance[idx,ahCut1005_allSamples]
data$counts <- data$counts[idx,ahCut1005_allSamples]
data$length <- data$length[idx,ahCut1005_allSamples]

DEGs_allSamples(coldata=coldata, i= i, data= data)

################################################################################

############################ -- filt samples -- ################################

# -- funcoes -- ################################################################
DEGs_filtSamples <- function(coldata, i, data){
  require(DESeq2)
  # -- coldata
  if(i == "ahCutUncut0_filtSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_02")]
    
  } else if(i == "ctCut1005_filtSamples" | i == "ahCut50_filtSamples" |
            i == "ahCut1000_filtSamples" | i == "ahCut1005_filtSamples"){
    coldata_samples <- coldata[coldata$names %in% colnames(data$counts),
                               c("names","Trat_03")]
  }else{
    cat("verifique a comparação 1")
  }
  
  # -- preparacao do objeto
  if(i == "ahCutUncut0_filtSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_02)
    
  } else if(i == "ctCut1005_filtSamples" | i == "ahCut50_filtSamples" |
            i == "ahCut1000_filtSamples" | i == "ahCut1005_filtSamples"){
    ddsTxi <- DESeqDataSetFromTximport(txi = data,
                                       colData = coldata_samples,
                                       design = ~ Trat_03)
  }else{
    cat("verifique a comparação 2")
  }
  
  
  # -- pre-filtragem
  keep <- rowSums(counts(ddsTxi)) >= 10
  dds <- ddsTxi[keep,]
  
  # -- setting the reference
  if(i == "ahCutUncut0_filtSamples"){
    dds$Trat_02 <- relevel(dds$Trat_02, ref = "uncut")
    
  } else if(i == "ahCut50_filtSamples" | i == "ahCut1000_filtSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "0")
    
  } else if(i == "ctCut1005_filtSamples" | i == "ahCut1005_filtSamples"){
    dds$Trat_03 <- relevel(dds$Trat_03, ref = "5")
  }else{
    cat("verifique a comparação 3")
  }
  
  # -- differential expression
  if(i == "ahCutUncut0_filtSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_02","cut","uncut"), alpha = 0.05)
    
  } else if(i == "ahCut50_filtSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","5","0"), alpha = 0.05)
    
  }else if(i == "ahCut1000_filtSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","0"), alpha = 0.05)
    
  } else if(i == "ctCut1005_filtSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)
    
  }else if(i == "ahCut1005_filtSamples"){
    dds <- DESeq(dds)
    res_padj05_lfc0 <- results(dds, contrast = c("Trat_03","100","5"), alpha = 0.05)
    
  }else{
    cat("verifique a comparação 4")
  }
  
  
  if(i == "ahCutUncut0_filtSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/filtSamples/ahCutUncut0_filtSamples_res05_sig_fc0.csv")
  } else if(i == "ahCut50_filtSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/filtSamples/ahCut50_filtSamples_res05_sig_fc0.csv")
  } else if(i == "ahCut1000_filtSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/filtSamples/ahCut1000_filtSamples_res05_sig_fc0.csv")
  } else if(i == "ctCut1005_filtSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/filtSamples/ctCut1005_filtSamples_res05_sig_fc0.csv")
  }else if(i == "ahCut1005_filtSamples"){
    write.csv(res_padj05_lfc0,
              file = "five_comp/filtSamples/ahCut1005_filtSamples_res05_sig_fc0.csv")
  } else{
    cat("verifique a comparacao 5")
  }
}

################################################################################
# -- Anotacao -- 
load("EnsDbAnnotation_20230519_atual.RData")
data_protCod <- EnsDbAnnotation[EnsDbAnnotation$gene_biotype=="protein_coding",]
head(data_protCod)
table(data_protCod$gene_biotype) # 30153

# -- colldata -- 
coldata <- readxl::read_xlsx("coldata.xlsx")
str(coldata)
coldata <- as.data.frame(coldata)

for (i in 2:ncol(coldata)){
  coldata[,i] <- as.factor(coldata[,i])
}

str(coldata)

# -- preparando objeto -- ######################################################
load("matrix_salmon_tximport_20230519.RData")

## -- matrizes com allSamples 
data_filtSamples <- mat_gse; rm(mat_gse)

head(rownames(data_filtSamples$counts)) # sem num da versao
head(rownames(data_protCod)) # sem num da versao

nrow(data_filtSamples$counts)

idx <- match(rownames(data_protCod), rownames(data_filtSamples$counts))
table(is.na(idx))
idx <- idx[!is.na(idx)]

colnames(data_filtSamples$counts)

ahCutUncut0_filtSamples <- grep("^AH.*_0_.", colnames(data_filtSamples$counts)) # ahCutUncut0_allSamples

ahCut50_filtSamples <- grep("AH_Cut_._.", colnames(data_filtSamples$counts)) # ahCut50_allSamples

ahCut1000_filtSamples <- grep("AH_Cut_.*0_.", colnames(data_filtSamples$counts)) # ahCut1000_allSamples

ctCut1005_filtSamples <- c(grep("CT_Cut_100_.", colnames(data_filtSamples$counts)),
                          grep("CT_Cut_5_.", colnames(data_filtSamples$counts)))# ctCut1005_allSamples

ahCut1005_filtSamples <- c(grep("AH_Cut_100_.", colnames(data_filtSamples$counts)),
                          grep("AH_Cut_5_.", colnames(data_filtSamples$counts)))# ahCut1005_allSamples

# -- expressao diferencial: all Samples -- #####################################

## -- ahCutUncut0_filtSamples -- 
data <- data_filtSamples
i="ahCutUncut0_filtSamples"
data$abundance <- data$abundance[idx,ahCutUncut0_filtSamples]
data$counts <- data$counts[idx,ahCutUncut0_filtSamples]
data$length <- data$length[idx,ahCutUncut0_filtSamples]

DEGs_filtSamples(coldata=coldata, i= i, data= data)

## -- ahCut50_filtSamples
data <- data_filtSamples
i="ahCut50_filtSamples"
data$abundance <- data$abundance[idx,ahCut50_filtSamples]
data$counts <- data$counts[idx,ahCut50_filtSamples]
data$length <- data$length[idx,ahCut50_filtSamples]

DEGs_filtSamples(coldata=coldata, i= i, data= data)

## -- ahCut1000_filtSamples
data <- data_filtSamples
i="ahCut1000_filtSamples"
data$abundance <- data$abundance[idx,ahCut1000_filtSamples]
data$counts <- data$counts[idx,ahCut1000_filtSamples]
data$length <- data$length[idx,ahCut1000_filtSamples]

DEGs_filtSamples(coldata=coldata, i= i, data= data)

## -- ctCut1005_filtSamples
data <- data_filtSamples
i="ctCut1005_filtSamples"
data$abundance <- data$abundance[idx,ctCut1005_filtSamples]
data$counts <- data$counts[idx,ctCut1005_filtSamples]
data$length <- data$length[idx,ctCut1005_filtSamples]

DEGs_filtSamples(coldata=coldata, i= i, data= data)

## -- ahCut1005_allSamples
data <- data_filtSamples
i="ahCut1005_filtSamples"
data$abundance <- data$abundance[idx,ahCut1005_filtSamples]
data$counts <- data$counts[idx,ahCut1005_filtSamples]
data$length <- data$length[idx,ahCut1005_filtSamples]

DEGs_filtSamples(coldata=coldata, i= i, data= data)

################################ -- fim -- #####################################
