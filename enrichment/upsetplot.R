################################################################################
## RNA-Seq - zebrafish
## UP set plot
## date of creation: 2023/08/14
## date of last update: 2023/08/14
## jean resende
################################################################################
library(VennDetail)
library(ggplot2)

filtGene <- function(x,y){
  data = read.csv(paste("../differential_expression/", y, "/", x,
                        "_", y, "_res05_sig_fc0.csv", sep = ""))
  datafilt = data[data$log2FoldChange < (-0.5) | data$log2FoldChange > 0.5,]
  datafilt = datafilt[datafilt$padj < 0.05,]
  datafilt = datafilt[!is.na(datafilt$padj),]
  return(datafilt$X)
}

# -- all samples -- ############################################################
ahctCut0_allSamples <- filtGene("ahctCut0","allSamples")
ahctCut100_allSamples <- filtGene("ahctCut100","allSamples")
ahctCut5_allSamples <- filtGene("ahctCut5","allSamples")
ahctUncut0_allSamples <- filtGene("ahctUncut0","allSamples")
ctCut1000_allSamples <- filtGene("ctCut1000","allSamples")
ctCut50_allSamples <- filtGene("ctCut50","allSamples")
ctCutUncut0_allSamples <- filtGene("ctCutUncut0","allSamples")

venn <- venndetail(list(ahctCut0 = ahctCut0_allSamples,
                        ahctCut100 = ahctCut100_allSamples,
                        ahctCut5 = ahctCut5_allSamples,
                        ahctUncut0 = ahctUncut0_allSamples,
                        ctCut1000 = ctCut1000_allSamples,
                        ctCut50 = ctCut50_allSamples,
                        ctCutUncut0 = ctCutUncut0_allSamples))
plot(venn, type = "upset")

pdf(file = "upsetplot_allsamples.pdf", width = 9, height = 5)
plot(venn, type = "upset", cex = 0.5, text.scale = c(1,1,1,1,1,1))
dev.off()

venn_select <- venndetail(list(ahctUncut0 = ahctUncut0_allSamples,
                               ahctCut0 = ahctCut0_allSamples,
                               ctCutUncut0 = ctCutUncut0_allSamples,
                               ahctCut5 = ahctCut5_allSamples,
                               ahctCut100 = ahctCut100_allSamples))
plot(venn_select, type = "upset")

pdf(file = "upsetplot_allsamples_select.pdf", width = 9, height = 5)
plot(venn_select, type = "upset", cex = 0.5, text.scale = c(1,1,1,1,1,1))
dev.off()

# -- filt samples -- ###########################################################
ahctCut0_filtSamples <- filtGene("ahctCut0","filtSamples")
ahctCut100_filtSamples <- filtGene("ahctCut100","filtSamples")
ahctCut5_filtSamples <- filtGene("ahctCut5","filtSamples")
ahctUncut0_filtSamples <- filtGene("ahctUncut0","filtSamples")
ctCut1000_filtSamples <- filtGene("ctCut1000","filtSamples")
ctCut50_filtSamples <- filtGene("ctCut50","filtSamples")
ctCutUncut0_filtSamples <- filtGene("ctCutUncut0","filtSamples")

venn <- venndetail(list(ahctCut0 = ahctCut0_filtSamples,
                        ahctCut100 = ahctCut100_filtSamples,
                        ahctCut5 = ahctCut5_filtSamples,
                        ahctUncut0 = ahctUncut0_filtSamples,
                        ctCut1000 = ctCut1000_filtSamples,
                        ctCut50 = ctCut50_filtSamples,
                        ctCutUncut0 = ctCutUncut0_filtSamples))
plot(venn, type = "upset")

pdf(file = "upsetplot_filtsamples.pdf", width = 9, height = 5)
plot(venn, type = "upset", cex = 0.5, text.scale = c(1,1,1,1,1,1))
dev.off()

venn_select <- venndetail(list(ahctUncut0 = ahctUncut0_filtSamples,
                               ahctCut0 = ahctCut0_filtSamples,
                               ctCutUncut0 = ctCutUncut0_filtSamples,
                               ahctCut5 = ahctCut5_filtSamples,
                               ahctCut100 = ahctCut100_filtSamples))
plot(venn_select, type = "upset")

pdf(file = "upsetplot_filtsamples_select.pdf", width = 9, height = 5)
plot(venn_select, type = "upset", cex = 0.5, text.scale = c(1,1,1,1,1,1))
dev.off()
