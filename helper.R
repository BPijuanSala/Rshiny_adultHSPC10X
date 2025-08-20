########################################################################################
## Title: helper.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for HSCs (NK Wilson and FK Hamey)
## Date: 22 November 2017
########################################################################################
#setwd("/Users/blancap/Documents/PhD_CAMBRIDGE/PhD/Experiments/PhD_BPS39/scripts/website3/")
# helper.R
library(ggplot2)
library(Matrix)
library(data.table)
load("./www/data/data.RData")

#Plot v1
#redRamp <- colorRampPalette(c("#ffe6e6", "#e80700"))

#Plot v2 - 19.01.2018
redRamp <- colorRampPalette(c("#ffe6e6", "#e80700"))

#xl <- -15; yb <- -45; xr <- 5; yt <- -40


#Create palette for clusters in LK WT and W41
colorPal <- c("#5b113c", #0
             "#0505ce", #1
             "#ffff33", #2
             "#b2df8a", #3
             "#b599d7", #4
             "#fb9a99", #5
             "#608eff", #6
             "#c9b1a9", #7
             "#ff7f00", #8
             "#e704c4", #9
             "#e31a1c", #10
             "#33a02c", #11
             "#6a3d9a") #12
  names(colorPal)<- seq(0,12,1)
  
  #This is the order I would like to obtain
  vioOrder <- c(3,  9,  5,  4, 12,  1,  2,  8, 10,  6, 11,  7,  0)
  vioOrderPos <- vioOrder+1
  colOrder <- colorPal[vioOrderPos]
  
  
###------------------
##   geneExpression Levels
###------------------


plotGeneLevels <- function(gene="Procr",coordDataset,genes,nameDataset="WT_LK",barcodes) {
  if((substr(gene,1,7)== "ENSMUSG")==TRUE){
    id <- gene
    geneCorrected <- as.character(geneTable[geneTable[,"Gene.ID"] == gene,"Associated.Gene.Name"])
  } else {
    geneCorrected <- paste0(toupper(substr(gene,1,1)), substr(gene,2,nchar(gene)))
    id <- as.character(geneTable[toupper(geneTable[,"Associated.Gene.Name"]) == toupper(gene),"Gene.ID"])
  }
  load(file=paste0("./www/dataSubsets/",nameDataset,"_indices.RData"))
  idxGene <- indices[toupper(names(indices))==toupper(geneCorrected)]
  colGenes <- names(indices[which(indices==unname(idxGene))])
  idxGene2 <- which(toupper(colGenes)%in%toupper(geneCorrected))
  
  load(file=paste0("./www/dataSubsets/",nameDataset,"_",idxGene,".RData"))
  
  dataGene <- as.vector(countsSubset[,idxGene2])
  rm(countsSubset)
  rm(indices)
  
  names(dataGene) <- barcodes
  dataGeneOrder <- dataGene[rownames(coordDataset)]

  #==================================
  #Plot version v2 - Created 19.01.2018
  #==================================
  df <- data.frame(x = coordDataset$x, y = coordDataset$y, exp = log10(dataGeneOrder+1))
  df <- df[order(df$exp,decreasing=F),]
  dfsub <- df[df$exp>0,]
  
  
  interval <- findInterval(dfsub$exp, seq(min(dfsub$exp), 
                                          max(dfsub$exp), 
                                          (max(dfsub$exp)-min(dfsub$exp))/10))
  
  interval[interval==0]<-1
  colorsPlot <- redRamp(11)[interval]
  
  
  par(mar=c(8,4,8,4),xpd=NA)
  plot(df$x, df$y, col="lightgrey", pch=16, cex=0.75, xlab="", 
       ylab="", main=gene, axes=F, cex.main=1.5)
       
  box(bty="l")
  points(dfsub$x, dfsub$y, pch=20, cex=0.5, 
         col=colorsPlot)
  
  xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4))); yb <- (min(df$y))-(-0.34*min(df$y)); xr <- max(df$x); yt <- (min(df$y))-(min(df$y)*(-0.17))
  rect(xleft=head(seq(xl,xr,(xr-xl)/10),-1), ybottom = yb, 
       xright = tail(seq(xl,xr,(xr-xl)/10),-1), ytop = yt,
       col=redRamp(10), border=redRamp(10), lwd=0.5, xpd=NA)
  rect(xl,yb,xr,yt, xpd=NA)
  segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.0381), xpd=NA)
  text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.089), labels = round(seq(min(dfsub$exp), max(dfsub$exp), (max(dfsub$exp)-min(dfsub$exp))/5),1), cex=1, xpd=NA)
  text(stats::median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.216), labels = expression('log'[10]*' normalized counts + 1'),cex=1.2, xpd=NA)
  
  
  
  ###Version v1
  #par(mar=c(8,4,8,4),xpd=NA)
  #df <- data.frame(x = coordDataset$x, y = coordDataset$y, exp = log10(dataGeneOrder+1))
  #dfsub <- df[df$exp>0,]
  #plot(df$x, df$y, col="lightgrey", pch=16, cex=0.75,
  #     xlab="", ylab="", main=gene, axes=F, cex.main=1.5)
  #box(bty="l")
  #points(dfsub$x, dfsub$y, pch=20, cex=0.5, 
  #       col=redRamp(10)[findInterval(df$exp, seq(min(dfsub$exp), 
  #                                                max(dfsub$exp), 
  #                                                (max(dfsub$exp)-min(dfsub$exp))/10))] )
  
  #xl <- min(df$x)- (min(df$x)*(-5.57*10^(-4))); yb <- (min(df$y))-(-0.34*min(df$y)); xr <- max(df$x); yt <- (min(df$y))-(min(df$y)*(-0.17))
  #rect(xleft=head(seq(xl,xr,(xr-xl)/10),-1), ybottom = yb, 
  #     xright = tail(seq(xl,xr,(xr-xl)/10),-1), ytop = yt,
  #     col=redRamp(10), border=redRamp(10), lwd=0.5, xpd=NA)
  #rect(xl,yb,xr,yt, xpd=NA)
  #segments(seq(xl,xr,((xr-xl)/5)),yb,seq(xl,xr,((xr-xl)/5)),yb-(-yb*0.0381), xpd=NA)
  #text(seq(xl,xr,((xr-xl)/5)), yb-(-yb*0.089), labels = round(seq(min(dfsub$exp), max(dfsub$exp), (max(dfsub$exp)-min(dfsub$exp))/5),1), cex=1, xpd=NA)
  #text(median(seq(xl,xr,((xr-xl)/5))), yb-(-yb*0.216), labels = expression('log'[10]*' normalized counts + 1'),cex=1.2, xpd=NA)
  
  
  
  
  }




###------------------
#     Violinplot
###------------------
vioplotPerCluster <- function(gene="Procr",clusters,type="LK",genes,barcodes,nameDataset){
  if((substr(gene,1,7)== "ENSMUSG")==TRUE){
    id <- gene
    geneCorrected <- as.character(geneTable[geneTable[,"Gene.ID"] == gene,"Associated.Gene.Name"])
  } else {
    geneCorrected <- paste0(toupper(substr(gene,1,1)), substr(gene,2,nchar(gene)))
    id <- as.character(geneTable[toupper(geneTable[,"Associated.Gene.Name"]) == toupper(gene),"Gene.ID"])
  }
  load(file=paste0("./www/dataSubsets/",nameDataset,"_indices.RData"))
  idxGene <- indices[toupper(names(indices))==toupper(geneCorrected)]
  
  colGenes <- names(indices[which(indices==unname(idxGene))])
  idxGene2 <- which(toupper(colGenes)%in%toupper(geneCorrected))
  load(file=paste0("./www/dataSubsets/",nameDataset,"_",idxGene,".RData"))
  dataGene <- as.vector(countsSubset[,idxGene2])
  rm(countsSubset)
  rm(indices)
  
  names(dataGene) <- barcodes
  dataGeneOrder <- dataGene[names(clusters)]
  df <- data.frame(cluster=clusters, exp = log10(dataGeneOrder+1))
  rownames(df) <- names(clusters)
  par(mar=c(8,4,8,4),xpd=NA)
  if (type == "LSK"){
    setOrder <- df[order(match(as.character(df$cluster),c("red2","purple4"))),]
    df$cluster <- factor(df$cluster, levels = unique(df[rownames(setOrder),"cluster"]))
    
    ggplot(df, aes(x=cluster, y=exp,fill=cluster)) + 
      geom_violin(scale="width") + 
      scale_fill_manual(values=c(as.character(toupper(unique(df$cluster))))) +
      theme_classic() +
      guides(fill=FALSE) +
      scale_x_discrete(limit = as.character((unique(df$cluster))),
                       labels = c("LSK","LK")) +
      theme(
        #axis.text.x = element_blank(),
        #axis.ticks = element_blank(),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=60,hjust=1),
        axis.title=element_text(size=14)) +
      labs(x = "Cluster", y=expression('log'[10]*' normalized counts + 1')) +
      theme(plot.margin = unit(c(2,1,1,1), "cm"))
    
  } else {
    setOrder <- df[order(match(as.character(df$cluster),colOrder)),]
    #Here I am reordering the factor according the order I wish to obtain, 
    #defined at the beginning of the script.
    df$cluster <- factor(df$cluster, levels = unique(df[rownames(setOrder),"cluster"]))
    ggplot(df, aes(x=cluster, y=exp,fill=cluster)) + 
      geom_violin(scale="width") + 
      scale_fill_manual(values=c(as.character(toupper((colOrder))))) +
      theme_classic() +
      guides(fill=FALSE) +
      scale_x_discrete(limit = as.character(colOrder),
                       labels = seq(1,13,1)) +
      theme(
        #axis.text.x = element_blank(),
        #axis.ticks = element_blank(),
        axis.text=element_text(size=12),
        axis.text.x=element_text(angle=60,hjust=1),
        axis.title=element_text(size=14)) +
      labs(x = "Cluster", y=expression('log'[10]*' normalized counts + 1')) +
      theme(plot.margin = unit(c(2,1,1,1), "cm"))
    
  }
  
}

