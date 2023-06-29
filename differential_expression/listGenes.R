################################################################################
# Name: listGenes                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/06/28                                                    #
# Last update date: 2023/06/28                                                 #
# Description: List of differential genes                                      #
################################################################################
# -- all samples
ahctCut0_allSamples <- read.csv("allSamples/AHcut0_CTcut0_res05_sig_fc0.csv")
ahctUncut0_allSamples <- read.csv("allSamples/AHuncut0_CTuncut0_res05_sig_fc0.csv")
ctCutUncut0_allSamples <- read.csv("allSamples/CTuncut0_CTcut0_res05_sig_fc0.csv")
ahctCut5_allSamples <- read.csv("allSamples/AHcut5_CTcut5_res05_sig_fc0.csv")
ahctCut100_allSamples <- read.csv("allSamples/AHcut100_CTcut100_res05_sig_fc0.csv")
ctCut0100_allSamples <- read.csv("allSamples/CTcut100_CTcut0_res05_sig_fc0.csv")
ctCut05_allSamples <- read.csv("allSamples/CTcut5_CTcut0_res05_sig_fc0.csv")

colnames(ahctCut0_allSamples)
rownames(ahctCut0_allSamples)

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
