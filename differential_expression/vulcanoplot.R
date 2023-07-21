################################################################################
# Name: vulcanoplot                                                 #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/07/03                                                    #
# Last update date: 2023/07/03                                                 #
# Description: graphic differential expression                                 #
################################################################################
library(ggplot2)
library(ggrepel)

data <- read.csv("filtSamples/CTcut5_CTcut0_res05_sig_fc0.csv")
is.data.frame(data)
rownames(data) <- data$X

data <- data[data$padj < 0.05,]
data <- na.omit(data)
View(data)

load("EnsDbAnnotation_20230519_atual.RData")

symbol_idx <- match(rownames(data), EnsDbAnnotation$ensemblid)
symbol_idx <- EnsDbAnnotation$symbol[symbol_idx]
table(is.na(symbol_idx))
data$symbol <- symbol_idx

data$diffexpressed <- "NO"
data$diffexpressed[data$log2FoldChange < (-0.5)] <- "DOWN"
data$diffexpressed[data$log2FoldChange > (0.5)] <- "UP"

#length(data$diffexpressed == "UP")
#length(data$diffexpressed == "DOWN")

data$delabel <- NA
data$delabel[data$diffexpressed != "NO"] <- data$symbol[data$diffexpressed != "NO"]

mycolors <- c("#00AFBB","#bb0c00","grey")
names(mycolors) <- c("DOWN","UP","NO")

ggplot(data, aes(x=log2FoldChange, y= -log10(padj), col=diffexpressed, label=delabel))+
  geom_point(size=0.5)+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = mycolors)
