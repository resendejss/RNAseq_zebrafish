################################################################################
# Name: listDEGS                                                               #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/06/28                                                    #
# Last update date: 2023/08/01                                                 #
# Description: List of differential genes                                      #
################################################################################

# -- all samples -- ############################################################
list.files("allSamples/")[1:7]

ahctCut0_allSamples <- read.csv("allSamples/ahctCut0_allSamples_res05_sig_fc0.csv")
ahctCut100_allSamples <- read.csv("allSamples/ahctCut100_allSamples_res05_sig_fc0.csv")
ahctCut5_allSamples <- read.csv("allSamples/ahctCut5_allSamples_res05_sig_fc0.csv")
ahctUncut0_allSamples <- read.csv("allSamples/ahctUncut0_allSamples_res05_sig_fc0.csv")
ctCut1000_allSamples <- read.csv("allSamples/ctCut1000_allSamples_res05_sig_fc0.csv")
ctCut50_allSamples <- read.csv("allSamples/ctCut50_allSamples_res05_sig_fc0.csv")
ctCutUncut0_allSamples <- read.csv("allSamples/ctCutUncut0_allSamples_res05_sig_fc0.csv")

colnames(ahctCut0_allSamples)
head(rownames(ahctCut0_allSamples))

lista <- mget(ls())

for (i in 1:length(lista)) {
  #rownames(lista[[i]]) <- lista[[i]][["X"]]
  lista[[i]]$expDiff <- rep(names(lista[i]), nrow(lista[[i]]))
  lista[[i]] <- lista[[i]][c(1,3,6:8)]
  colnames(lista[[i]]) <- paste(colnames(lista[[i]]),
                                lista[[i]]$expDiff[1],
                                sep = "_")
  colnames(lista[[i]])[1] <- "geneid"
}


data <- Reduce(function(x,y) merge(x,y,by="geneid",all=TRUE),lista)
rownames(data) <- data$geneid

load("EnsDbAnnotation_20230519_atual.RData")
idx <- match(rownames(data), EnsDbAnnotation$ensemblid) 
table(is.na(idx))

data$symbol <- EnsDbAnnotation$symbol[idx]
data <- data[,c(1,30,2:29)]

data <- data[,-c(6,10,14,18,22,26,30)]
colnames(data)[1] <- "ensemblid"

write.csv(data, file="listDEGs/LIST_differentialExpression_allSamples.csv")
rm(list = ls())

# -- filt samples -- ###########################################################
list.files("filtSamples/")[1:7]

ahctCut0_filtSamples <- read.csv("filtSamples/ahctCut0_filtSamples_res05_sig_fc0.csv")
ahctCut100_filtSamples <- read.csv("filtSamples/ahctCut100_filtSamples_res05_sig_fc0.csv")
ahctCut5_filtSamples <- read.csv("filtSamples/ahctCut5_filtSamples_res05_sig_fc0.csv")
ahctUncut0_filtSamples <- read.csv("filtSamples/ahctUncut0_filtSamples_res05_sig_fc0.csv")
ctCut1000_filtSamples <- read.csv("filtSamples/ctCut1000_filtSamples_res05_sig_fc0.csv")
ctCut50_filtSamples <- read.csv("filtSamples/ctCut50_filtSamples_res05_sig_fc0.csv")
ctCutUncut0_filtSamples <- read.csv("filtSamples/ctCutUncut0_filtSamples_res05_sig_fc0.csv")

colnames(ahctCut0_filtSamples)
head(rownames(ahctCut0_filtSamples))

lista <- mget(ls())

for (i in 1:length(lista)) {
  #rownames(lista[[i]]) <- lista[[i]][["X"]]
  lista[[i]]$expDiff <- rep(names(lista[i]), nrow(lista[[i]]))
  lista[[i]] <- lista[[i]][c(1,3,6:8)]
  colnames(lista[[i]]) <- paste(colnames(lista[[i]]),
                                lista[[i]]$expDiff[1],
                                sep = "_")
  colnames(lista[[i]])[1] <- "ensemblid"
}


data <- Reduce(function(x,y) merge(x,y,by="ensemblid",all=TRUE),lista)
rownames(data) <- data$ensemblid

load("EnsDbAnnotation_20230519_atual.RData")
idx <- match(rownames(data), EnsDbAnnotation$ensemblid) 
table(is.na(idx))

data$symbol <- EnsDbAnnotation$symbol[idx]
data <- data[,c(1,30,2:29)]

data <- data[,-c(6,10,14,18,22,26,30)]
#colnames(data)[1] <- "ensemblid"

write.csv(data, file="listDEGs/LIST_differentialExpression_filtSamples.csv")
rm(list = ls())
################################################################################