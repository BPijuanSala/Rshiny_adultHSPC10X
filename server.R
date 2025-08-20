########################################################################################
## Title: server.R
## Author: Blanca Pijuan-Sala
## Description: Shiny app for HSCs (NK Wilson and FK Hamey)
## Date: 22 November 2017 
########################################################################################




library(shiny)
#library(DT)
library(png)
library(Matrix)
library(data.table)

source("helper.R")




server <- function(input, output) {

  #Set images
  output$LSKLK_colLSK <- renderImage({
    outfile <- tempfile(fileext='.png')
    p <- readPNG("./www/images/LSKLK_colLSK.png")
    writePNG(p, target=outfile)
    list(src = outfile, contentType = 'image/png', width = 370, height = 450)
  }, deleteFile = TRUE)
  
  output$LSKLK_colLK <- renderImage({
    outfile <- tempfile(fileext='.png')
    p <- readPNG("./www/images/LSKLK_colLK.png")
    writePNG(p, target=outfile)
    list(src = outfile, contentType = 'image/png', width = 370, height = 450)
  }, deleteFile = TRUE)
  
  
  output$LK_WT <- renderImage({
    outfile <- tempfile(fileext='.png')
    p <- readPNG("./www/images/LK_WT.png")
    writePNG(p, target=outfile)
    list(src = outfile, contentType = 'image/png', width = 370, height = 450)
  }, deleteFile = TRUE)
  
  
  output$LK_W41 <- renderImage({
    outfile <- tempfile(fileext='.png')
    p <- readPNG("./www/images/LK_W41.png")
    writePNG(p, target=outfile)
    list(src = outfile, contentType = 'image/png', width = 370, height = 450)
  }, deleteFile = TRUE)
  
  
  
  
  
  
  ####----------------------------------
  ##   WT LSK LK GeneLevels + download
  ####----------------------------------
  
  output$WTLK_LSK_plot <- renderPlot({ 
    validate(
      if((substr(input$gene,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene)),"Associated.Gene.Name"]%in%geneNames.WT.LSK.LK, 'This gene does not exist or is not expressed. Please choose another one.')
        
      } else {
        geneCorrected <- toupper(input$gene)
        need(geneCorrected %in% toupper(geneNames.WT.LSK.LK), 'This gene does not exist or is not expressed. Please choose another one.')
        
      })
    plotGeneLevels(gene=input$gene,
                   coordDataset=coords.WT.LSK.LK,genes=geneNames.WT.LSK.LK,nameDataset="WT_LSK",
                   barcodes=barcode.WT.LSK.LK)
  }, height=500, width=400)
  
  
  
  
  output$download_WTLK_LSK_plot <- downloadHandler(
    filename = function() { paste0(input$gene, '_expression_WTLSKLK_plot.pdf') },
    content = function(file) {
      pdf(file, width = 7, height = 8)
      plotGeneLevels(gene=input$gene,
                     coordDataset=coords.WT.LSK.LK,genes=geneNames.WT.LSK.LK,
                     nameDataset="WT_LSK",barcodes=barcode.WT.LSK.LK)
      dev.off()
    }
  )
  
  
  ####----------------------------------
  ##   WT LK GeneLevels + download
  ####----------------------------------
  
  
  output$WTLK_plot <- renderPlot({ 
    validate(
      if((substr(input$gene2,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene2)),"Associated.Gene.Name"]%in%geneNames.WT.LK, 'This gene does not exist or is not expressed. Please choose another one.')
        
      } else {
        geneCorrected <- toupper(input$gene2)
        need(geneCorrected %in% toupper(geneNames.WT.LK), 'This gene does not exist or is not expressed. Please choose another one.')
        
      })
    plotGeneLevels(gene=input$gene2,
                   coordDataset=coords.WT.LK,genes=geneNames.WT.LK,
                   nameDataset="WT_LK",barcodes=barcode.WT.LK)
  }, height=460, width=380)
  
  
  
  output$download_WTLK_plot <- downloadHandler(
    filename = function() { paste0(input$gene2, '_expression_WTLK_plot.pdf') },
    content = function(file) {
      pdf(file, width = 7, height = 8)
      plotGeneLevels(gene=input$gene2,
                     coordDataset=coords.WT.LK,genes=geneNames.WT.LK,
                     nameDataset="WT_LK",barcodes=barcode.WT.LK)
      dev.off()
    }
  )
  
  
  ####----------------------------------
  ##   W41 LK GeneLevels + download
  ####----------------------------------
  
  output$W41LK_plot <- renderPlot({ 
    validate(
      if((substr(input$gene2,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene2)),"Associated.Gene.Name"]%in%geneNames.W41.LK, '')
        
      } else {
        geneCorrected <- toupper(input$gene2)
        need(geneCorrected %in% toupper(geneNames.W41.LK), '')
        
      })
    plotGeneLevels(gene=input$gene2,
                   coordDataset=coords.W41.LK,genes=geneNames.W41.LK,
                   nameDataset="W41_LK",barcodes=barcode.W41.LK
                   )
  }, height=460, width=380)
  
  
  
  output$downloadPlot_W41LK_plot <- downloadHandler(
    filename = function() { paste0(input$gene2, '_expression_W41LK_plot.pdf') },
    content = function(file) {
      pdf(file, width = 7, height = 8)
      plotGeneLevels(gene=input$gene2,coordDataset=coords.W41.LK,
                     genes=geneNames.W41.LK, nameDataset="W41_LK",barcodes=barcode.W41.LK)
      dev.off()
    }
  )
  
  
  
  ####----------------------------------
  ##   WT LSK LK BOXPLOTS + download
  ####----------------------------------  
  
  
  output$WTLK_LSK_boxplot <- renderPlot({
    validate(
      if((substr(input$gene,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene)),"Associated.Gene.Name"]%in%geneNames.WT.LSK.LK, '')
        
      } else {
        geneCorrected <- toupper(input$gene)
        need(geneCorrected %in% toupper(geneNames.WT.LSK.LK), '')
        
      })
    vioplotPerCluster(gene=input$gene,clusters=clusters.WT.LSK.LK.vec,
                      type="LSK",genes=geneNames.WT.LSK.LK,
                      barcodes=barcode.WT.LSK.LK,nameDataset="WT_LSK")
  }, height=400, width=420)
  
  
  
  output$download_WTLK_LSK_boxplot <- downloadHandler(
    filename = function() { paste0(input$gene, '_expression_WTLSKLK_boxplot.pdf') },
    content = function(file) {
      pdf(file, width = 8, height = 8)
      vioplotPerCluster(gene=input$gene,clusters=clusters.WT.LSK.LK.vec,
                        type="LSK",
                        genes=geneNames.WT.LSK.LK,barcodes=barcode.WT.LSK.LK,
                        nameDataset="WT_LSK")
      dev.off()
    }
  )
  
  
  ####----------------------------------
  ##   WT LK BOXPLOTS + download
  ####----------------------------------
  
  output$WTLK_boxplot <- renderPlot({
    validate(
      if((substr(input$gene2,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene2)),"Associated.Gene.Name"]%in%geneNames.WT.LK, '')
        
      } else {
        geneCorrected <- toupper(input$gene2)
        need(geneCorrected %in% toupper(geneNames.WT.LK), '')
        
      })
    vioplotPerCluster(gene=input$gene2,clusters=clusters.WT.LK.vec,
                      genes=geneNames.WT.LK,barcodes=barcode.WT.LK,
                      nameDataset="WT_LK")
  }, height=380, width=400)
  
  
  output$download_WTLK_boxplot <- downloadHandler(
    filename = function() { paste0(input$gene2, '_expression_WTLK_boxplot.pdf') },
    content = function(file) {
      pdf(file, width = 8, height = 8)
      vioplotPerCluster(gene=input$gene2,clusters=clusters.WT.LK.vec,
                       genes=geneNames.WT.LK,barcodes=barcode.WT.LK,
                        nameDataset="WT_LK")
      dev.off()
    }
  )
  
  
  ####----------------------------------
  ##   W41 LK BOXPLOTS + download
  ####----------------------------------
  
  output$W41LK_boxplot <- renderPlot({
    validate(
      if((substr(input$gene2,1,7)== "ENSMUSG")==TRUE){
        need(geneTable[geneTable[,"Gene.ID"]==toupper(as.character(input$gene2)),"Associated.Gene.Name"]%in%geneNames.W41.LK, '')
        
      } else {
        geneCorrected <- toupper(input$gene2)
        need(geneCorrected %in% toupper(geneNames.W41.LK), '')
        
      })
    vioplotPerCluster(gene=input$gene2,clusters=clusters.W41.LK.vec,
                      genes=geneNames.W41.LK,barcodes=barcode.W41.LK,
                      nameDataset="W41_LK")
  }, height=380, width=400)
  
  output$download_W41LK_boxplot <- downloadHandler(
    filename = function() { paste0(input$gene2, '_expression_W41LK_boxplot.pdf') },
    content = function(file) {
      pdf(file, width = 8, height = 8)
      vioplotPerCluster(gene=input$gene2,clusters=clusters.W41.LK.vec,
                        genes=geneNames.W41.LK,barcodes=barcode.W41.LK,
                        nameDataset="W41_LK")
      dev.off()
    }
  )
  

  
  
  ###-----------------------------
  ##Download data
  ###----------------------
  
  output$download_WT_LSK_LK_data <- downloadHandler(
    
    filename <- function() {
      "LSK_data.tar.gz"
    },
    
    content <- function(file) {
      file.copy("./www/data/LSK_data.tar.gz", file)
    },
    contentType = "application/zip"
  )
  
  output$download_WT_LK_data <- downloadHandler(
    
    filename <- function() {
      "LK_WT_data.tar.gz"
    },
    
    content <- function(file) {
      file.copy("./www/data/LK_WT_data.tar.gz", file)
    },
    contentType = "application/zip"
  )  

  output$download_W41_LK_data <- downloadHandler(
    
    filename <- function() {
      "LK_W41_data.tar.gz"
    },
    
    content <- function(file) {
      file.copy("./www/data/LK_W41_data.tar.gz", file)
    },
    contentType = "application/zip"
  )  
  
  
  
  }







