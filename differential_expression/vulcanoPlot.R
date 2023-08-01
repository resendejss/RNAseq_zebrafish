################################################################################
# Name: vulcanoplot                                                            #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Project: RNAseq_zebrafish                                                    #
#                                                                              #
# Creation date: 2023/07/03                                                    #
# Last update date: 2023/07/31                                                 #
# Description: graphic differential expression                                 #
################################################################################
vulcanoPlot <- function(allFilt){
  require(ggplot2)
  require(ggrepel)
  require(magrittr)
  
  files <- list.files(allFilt)
  
  for (i in files) {
    data <- read.csv(paste(allFilt, i, sep="/"))
    rownames(data) <- data$X
    
    
    data <- data[data$padj < 0.05,] %>% na.omit()
    
    load("EnsDbAnnotation_20230519_atual.RData")
    
    symbol_idx <- match(rownames(data), EnsDbAnnotation$ensemblid)
    data$symbol <- EnsDbAnnotation$symbol[symbol_idx]
    
    data$diffexpressed <- "NO"
    data$diffexpressed[data$log2FoldChange < (-0.5)] <- "DOWN"
    data$diffexpressed[data$log2FoldChange > (0.5)] <- "UP"
    
    data$delabel <- NA
    data$delabel[data$diffexpressed != "NO"] <- data$symbol[data$diffexpressed != "NO"]
    
    # -- salvar a tabela com os UP DOWN
    data_upDown <- data[data$diffexpressed != "NO",]
    write.csv(data_upDown, paste(allFilt,"/tableVulcanoplot_",i, sep = ""))
    
    mycolors <- c("#00AFBB","#bb0c00","grey")
    names(mycolors) <- c("DOWN","UP","NO")
    
    graph <- ggplot(data, aes(x=log2FoldChange, y= -log10(padj),
                              col=diffexpressed, label=delabel))+
      geom_point(size=0.5)+
      labs(color = 'DEGs', x=expression("log"[2]*"FC"), y=expression("-log"[10]*"padj"))+
      coord_cartesian(ylim = c(0,30), xlim = c(-5, 5)) +
      #theme_minimal()+
      geom_text_repel()+
      scale_color_manual(values = mycolors,
                         labels=c(paste("DOWN\n (",
                                        length(data$diffexpressed[data$diffexpressed=="DOWN"]),
                                        ")", sep = ""),
                                  paste("NO\n (", length(data$diffexpressed[data$diffexpressed=="NO"]),
                                        ")", sep = ""),
                                  paste("UP\n (",
                                        length(data$diffexpressed[data$diffexpressed=="UP"]),
                                        ")", sep = "")))+
      ggtitle(gsub(".csv","", i))+
      guides(color = guide_legend(label.position = "bottom",
                                  title.position="top",
                                  title.hjust=0.5,
                                  override.aes=list(size=3,
                                                    shape=21,
                                                    fill=c("#00AFBB","grey","#bb0c00"))))+
      theme(
        plot.title=element_text(hjust=0.5),
        legend.position="bottom",
        legend.margin = margin(-5, -5, 5, -5),
        legend.spacing.y = unit(0, 'cm'), 
        title = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.y=element_text(size=8), 
        axis.text.x=element_text(size=8),
        axis.ticks = element_line(linetype=1, color="grey"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10 ),
        legend.key.size = unit(0,"mm"),
        legend.background = element_rect(fill=NULL, color=NULL),
        axis.line = element_blank(),
        panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
        panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
        panel.background = element_rect(fill = "grey98", colour = "grey50"),
        panel.border = element_rect(colour = "grey", fill=NA, size=1))
    
    ggsave(paste(allFilt,"/vulcanoPlot_",gsub(".csv",".pdf",i),sep = ""),
           graph, width = 5, height = 5)
    
  }
}

################################################################################

vulcanoPlot(allFilt = "allSamples")
vulcanoPlot(allFilt = "filtSamples")

