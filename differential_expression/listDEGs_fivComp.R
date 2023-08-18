################################################################################
# Name: listDEGS_fivComp                                                       #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/08/18                                                    #
# Last update date: 2023/08/18                                                 #
# Description: List of differential genes                                      #
################################################################################

# -- all samples -- ############################################################
list.files("five_comp/allSamples/")[1:5]

ahCut1000_allSamples <- read.csv("five_comp/allSamples/ahCut1000_allSamples_res05_sig_fc0.csv")
ahCut1005_allSamples <- read.csv("five_comp/allSamples/ahCut1005_allSamples_res05_sig_fc0.csv")
ahCut50_allSamples <- read.csv("five_comp/allSamples/ahCut50_allSamples_res05_sig_fc0.csv")
ahCutUncut0_allSamples <- read.csv("five_comp/allSamples/ahCutUncut0_allSamples_res05_sig_fc0.csv")
ctCut1005_allSamples <- read.csv("five_comp/allSamples/ctCut1005_allSamples_res05_sig_fc0.csv")

colnames(ahCut1000_allSamples)
head(rownames(ahCut1000_allSamples))

rownames(ahCut1000_allSamples) <- ahCut1000_allSamples$X
rownames(ahCut1005_allSamples) <- ahCut1005_allSamples$X
rownames(ahCut50_allSamples) <- ahCut50_allSamples$X
rownames(ahCutUncut0_allSamples) <- ahCutUncut0_allSamples$X
rownames(ctCut1005_allSamples) <- ctCut1005_allSamples$X

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

colnames(data)
data <- data[,c(1,22,2:4,6:8,10:12,14:16,18:20)]

colnames(data)[1] <- "ensemblid"

write.csv(data, file="listDEGs/LIST_differentialExpression_allSamples_fiveComp.csv")
rm(list = ls())

# -- filt samples -- ###########################################################

list.files("five_comp/filtSamples/")[1:5]

ahCut1000_filtSamples <- read.csv("five_comp/filtSamples/ahCut1000_filtSamples_res05_sig_fc0.csv")
ahCut1005_filtSamples <- read.csv("five_comp/filtSamples/ahCut1005_filtSamples_res05_sig_fc0.csv")
ahCut50_filtSamples <- read.csv("five_comp/filtSamples/ahCut50_filtSamples_res05_sig_fc0.csv")
ahCutUncut0_filtSamples <- read.csv("five_comp/filtSamples/ahCutUncut0_filtSamples_res05_sig_fc0.csv")
ctCut1005_filtSamples <- read.csv("five_comp/filtSamples/ctCut1005_filtSamples_res05_sig_fc0.csv")

colnames(ahCut1000_filtSamples)
head(rownames(ahCut1000_filtSamples))

rownames(ahCut1000_filtSamples) <- ahCut1000_filtSamples$X
rownames(ahCut1005_filtSamples) <- ahCut1005_filtSamples$X
rownames(ahCut50_filtSamples) <- ahCut50_filtSamples$X
rownames(ahCutUncut0_filtSamples) <- ahCutUncut0_filtSamples$X
rownames(ctCut1005_filtSamples) <- ctCut1005_filtSamples$X

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

colnames(data)
data <- data[,c(1,22,2:4,6:8,10:12,14:16,18:20)]

colnames(data)[1] <- "ensemblid"

write.csv(data, file="listDEGs/LIST_differentialExpression_filtSamples_fiveComp.csv")
rm(list = ls())
################################################################################

list_allSamples_all <- read.csv("listDEGs/LIST_differentialExpression_allSamples.csv")
list_allSamples_fivComp <- read.csv("listDEGs/LIST_differentialExpression_allSamples_fiveComp.csv")

list_allSamples_fivComp <- list_allSamples_fivComp[,-c(1,3)]


list_all_allSamples <- merge(list_allSamples_all,list_allSamples_fivComp,
                             by="ensemblid", all = T)
list_all_allSamples <- list_all_allSamples[,-2]

head(list_all_allSamples$ensemblid)
head(list_allSamples_all$ensemblid)

write.csv(list_all_allSamples, file = "listDEGs/LIST_all_allSamples.csv")

rm(list = ls())

################################################################################

list_filtSamples_all <- read.csv("listDEGs/LIST_differentialExpression_filtSamples.csv")
list_filtSamples_fivComp <- read.csv("listDEGs/LIST_differentialExpression_filtSamples_fiveComp.csv")

list_filtSamples_fivComp <- list_filtSamples_fivComp[,-c(1,3)]


list_all_filtSamples <- merge(list_filtSamples_all,list_filtSamples_fivComp,
                             by="ensemblid", all = T)
list_all_filtSamples <- list_all_filtSamples[,-2]

head(list_all_filtSamples$ensemblid)
head(list_filtSamples_all$ensemblid)

write.csv(list_all_filtSamples, file = "listDEGs/LIST_all_filtSamples.csv")

rm(list = ls())
